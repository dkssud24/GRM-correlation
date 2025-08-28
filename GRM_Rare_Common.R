############################################################
## 0) 필수 패키지
############################################################
suppressWarnings(library(data.table))

############################################################
## 1) 헬퍼 함수들 (공통)
############################################################

.read_grm_ids <- function(prefix) {
  read.table(paste0(prefix, ".grm.id"), stringsAsFactors = FALSE)
}

# y를 공변량으로 잔차화 (없으면 평균중심)
.residualize_y <- function(y, covar = NULL) {
  y <- as.numeric(y)
  if (is.null(covar) || ncol(covar) == 0) {
    return(y - mean(y, na.rm = TRUE))
  } else {
    covar <- as.data.frame(covar)
    fit <- lm(y ~ ., data = covar, na.action = na.exclude)
    return(as.numeric(residuals(fit)))  # 결측행은 NA residual
  }
}

# 큰 n에서도 안전한 오프대각 (i>j) 샘플러
.sample_offdiag_idx <- function(n, K = 5e5, seed = 1, verbose = TRUE) {
  total_offdiag <- as.numeric(n) * (as.numeric(n) - 1) / 2
  K <- as.integer(min(K, total_offdiag))
  if (verbose) cat(sprintf("[INFO] n=%d; total off-diag=%s; K=%s\n",
                           n, format(total_offdiag, big.mark=","), format(K, big.mark=",")))
  set.seed(seed)
  want <- K; buf <- numeric(0)
  while (length(buf) < K) {
    batch <- as.integer(max(want * 1.3, 1e6))
    cand  <- floor(runif(batch, min = 0, max = total_offdiag))
    buf   <- unique(c(buf, cand))
    want  <- K - length(buf)
  }
  k0 <- sort(buf[seq_len(K)])  # 0-based

  # 0-based k -> (i,j)
  i_floor <- floor((1 + sqrt(1 + 8 * k0)) / 2)
  i <- i_floor + 1
  Tij <- i_floor * (i_floor - 1) / 2
  j <- as.numeric(k0 - Tij + 1)

  # 하삼각(대각 포함) 1-based 선형 인덱스
  lin_idx <- as.numeric(i) * (as.numeric(i) - 1) / 2 + as.numeric(j)
  max_L   <- as.numeric(n) * (as.numeric(n) + 1) / 2
  if (!all(is.finite(lin_idx))) stop("lin_idx non-finite")
  if (any(lin_idx < 1) || any(lin_idx > max_L)) stop("lin_idx out of range")

  list(i = as.integer(i), j = as.integer(j), lin_idx = lin_idx, K = K, n = n)
}

# 선형 인덱스 위치에서 GRM 값 읽기
.read_grm_at <- function(prefix, lin_idx, endian = "little", verbose = TRUE) {
  con <- file(paste0(prefix, ".grm.bin"), "rb")
  on.exit(close(con), add = TRUE)
  bytes <- 8L
  out <- numeric(length(lin_idx))
  if (verbose) cat("[INFO] Reading GRM values...\n")
  for (t in seq_along(lin_idx)) {
    where <- (as.numeric(lin_idx[t]) - 1) * bytes
    seek(con, where = where, origin = "start")
    val <- readBin(con, what = "numeric", n = 1L, size = bytes, endian = endian)
    out[t] <- if (length(val) == 0L) NA_real_ else val
    if (verbose && (t %% max(1L, floor(length(lin_idx)/10))) == 0L) {
      cat(sprintf("  ... %d / %d\n", t, length(lin_idx)))
    }
  }
  out
}

# (i,j) -> 하삼각 선형 인덱스 (permutation용)
.lt_index <- function(i, j) {
  ii <- pmax(i, j); jj <- pmin(i, j)
  as.numeric(ii) * (as.numeric(ii) - 1) / 2 + as.numeric(jj)
}

############################################################
## 2) 분석 A: GRM 태깅 R^2 (rare ~ common) + 간단 CI
############################################################

grm_tagging_R2 <- function(prefix_common, prefix_rare, K = 5e5,
                           seed = 1, endian = "little", verbose = TRUE) {
  ids <- .read_grm_ids(prefix_common)
  n <- nrow(ids)
  samp <- .sample_offdiag_idx(n, K = K, seed = seed, verbose = verbose)
  x <- .read_grm_at(prefix_common, samp$lin_idx, endian = endian, verbose = verbose)
  y <- .read_grm_at(prefix_rare,   samp$lin_idx, endian = endian, verbose = FALSE)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  R2 <- summary(lm(y ~ x))$r.squared
  list(R2 = R2, unique_frac = 1 - R2, K = length(x))
}

# 반복 리샘플로 CI (간단)
grm_tagging_R2_ci <- function(prefix_common, prefix_rare, K = 5e5,
                              B = 10, seed = 1, endian = "little") {
  set.seed(seed)
  vals <- numeric(B)
  for (b in 1:B) {
    res <- grm_tagging_R2(prefix_common, prefix_rare, K = K,
                          seed = seed + b, endian = endian, verbose = FALSE)
    vals[b] <- res$R2
  }
  ci <- quantile(vals, c(0.025, 0.5, 0.975))
  list(R2_mean = mean(vals), R2_median = unname(ci[2]),
       R2_CI = unname(ci[c(1,3)]), unique_frac_mean = 1 - mean(vals), B = B, K = K)
}

############################################################
## 3) 분석 B+C: 정직교(rare⊥) + HE 회귀 + permutation 대조
############################################################

he_ortho_and_perm <- function(prefix_common, prefix_rare,
                              y_vec, covar_df = NULL,
                              K = 1e6, seed = 42,
                              endian = "little", verbose = TRUE) {
  # 0) 기본 세팅/잔차화
  ids <- .read_grm_ids(prefix_common); n <- nrow(ids)
  y_star <- .residualize_y(y_vec, covar = covar_df)
  s2 <- var(y_star, na.rm = TRUE)

  # 1) (i,j) 샘플 분할 (학습/평가)
  samp <- .sample_offdiag_idx(n, K = K, seed = seed, verbose = verbose)
  idx <- seq_len(samp$K)
  set.seed(seed)
  idxA <- sample(idx, floor(length(idx)/2))
  idxB <- setdiff(idx, idxA)

  # 2) A세트에서 Gr ~ Gc 학습
  GcA <- .read_grm_at(prefix_common, samp$lin_idx[idxA], endian = endian, verbose = verbose)
  GrA <- .read_grm_at(prefix_rare,   samp$lin_idx[idxA], endian = endian, verbose = FALSE)
  fitA <- lm(GrA ~ GcA)
  a0 <- coef(fitA)[1]; a1 <- coef(fitA)[2]

  # 3) B세트에 적용 → rare⊥
  GcB <- .read_grm_at(prefix_common, samp$lin_idx[idxB], endian = endian, verbose = FALSE)
  GrB <- .read_grm_at(prefix_rare,   samp$lin_idx[idxB], endian = endian, verbose = FALSE)
  GrB_perp <- GrB - (a0 + a1 * GcB)
  tag_cor  <- suppressWarnings(cor(GcB, GrB_perp, use = "complete.obs"))

  # 4) HE 회귀 on B: z ~ Gc + rare⊥
  iB <- samp$i[idxB]; jB <- samp$j[idxB]
  zB <- y_star[iB] * y_star[jB]
  ok <- is.finite(zB) & is.finite(GcB) & is.finite(GrB_perp)
  df <- data.frame(z = zB[ok], Gc = GcB[ok], Grp = GrB_perp[ok])

  fit_r <- lm(z ~ Gc,       data = df)
  fit_f <- lm(z ~ Gc + Grp, data = df)

  co <- coef(summary(fit_f))
  beta_Gr_perp <- unname(co["Grp","Estimate"])
  p_Gr_perp    <- unname(co["Grp","Pr(>|t|)"])
  R2_r <- summary(fit_r)$r.squared
  R2_f <- summary(fit_f)$r.squared
  part_R2 <- (R2_f - R2_r) / max(1 - R2_r, .Machine$double.eps)
  h2_rare_perp <- beta_Gr_perp / s2

  # 5) permutation control: rare ID 섞기
  set.seed(seed + 123)
  perm <- sample.int(n, n, replace = FALSE)
  lin_perm_B <- .lt_index(perm[iB], perm[jB])
  GrB_perm <- .read_grm_at(prefix_rare, lin_perm_B, endian = endian, verbose = FALSE)
  okp <- is.finite(zB) & is.finite(GcB) & is.finite(GrB_perm)
  dfp <- data.frame(z = zB[okp], Gc = GcB[okp], Grp = GrB_perm[okp])
  fit_p <- lm(z ~ Gc + Grp, data = dfp)
  beta_perm <- coef(summary(fit_p))["Grp","Estimate"]
  p_perm    <- coef(summary(fit_p))["Grp","Pr(>|t|)"]

  list(
    kept_pairs = nrow(df),
    corr_Gc_Gr_perp = tag_cor,
    beta_Gr_perp = beta_Gr_perp,
    h2_rare_perp = h2_rare_perp,
    p_Gr_perp = p_Gr_perp,
    partial_R2 = part_R2,
    beta_perm = beta_perm,
    p_perm = p_perm
  )
}

############################################################
## 4) 실제 데이터 로드/정렬 (여기만 경로 맞춰서 실행)
############################################################

# ---- 경로 세팅
prefix_common <- "/data/BiO/hae/WGS_rev/0_common/v2.common.grm"
prefix_ld005  <- "/data/BiO/hae/WGS_rev/Z_LD005/v6.ABC.RARE.Aset"   # 실제 prefix 확인
prefix_ld010  <- "/data/BiO/hae/WGS_rev/2_LD01/v4.grm.LD01"         # 실제 prefix 확인

# ---- 파일 읽기
pheno_dt <- fread("X50", header = TRUE)                 # cols: FID IID norm_pheno
cov_dt   <- fread("cov.txt", header = FALSE)            # cols: FID IID <binary...> (없으면 자동 NULL)
qcov_dt  <- fread("qcov.Aset.Rare.pc", header = FALSE)  # cols: FID IID <PC...>    (없으면 자동 NULL)

# ---- 컬럼 이름 정리
if (!all(c("FID","IID") %in% names(pheno_dt))) setnames(pheno_dt, 1:2, c("FID","IID"))
y_col <- setdiff(names(pheno_dt), c("FID","IID"))[1]
setnames(pheno_dt, y_col, "PHENO")

setnames(cov_dt,  1:2, c("FID","IID"))
if (ncol(cov_dt)  > 2) setnames(cov_dt,  3:ncol(cov_dt),  paste0("COV",  seq_len(ncol(cov_dt)-2)))
setnames(qcov_dt, 1:2, c("FID","IID"))
if (ncol(qcov_dt) > 2) setnames(qcov_dt, 3:ncol(qcov_dt), paste0("QCOV", seq_len(ncol(qcov_dt)-2)))

# ---- .grm.id 순서로 정렬 (common 기준)
ids <- as.data.table(.read_grm_ids(prefix_common)); setnames(ids, 1:2, c("FID","IID"))
ph  <- merge(ids, pheno_dt, by = c("FID","IID"), all.x = TRUE, sort = FALSE)

cov_merged <- ids
if (ncol(cov_dt)  > 2) cov_merged <- merge(cov_merged, cov_dt,  by = c("FID","IID"), all.x = TRUE, sort = FALSE)
if (ncol(qcov_dt) > 2) cov_merged <- merge(cov_merged, qcov_dt, by = c("FID","IID"), all.x = TRUE, sort = FALSE)

# ---- covar_df 만들기 (없으면 NULL)
covar_df <- NULL
if (ncol(cov_merged) > 2) {
  covar_df <- as.data.frame(cov_merged[, -(1:2)])
  names(covar_df) <- make.names(names(covar_df), unique = TRUE)
}
y_vec <- as.numeric(ph$PHENO)

cat(sprintf("[CHECK] y non-finite count: %d / %d\n", sum(!is.finite(y_vec)), length(y_vec)))
if (!is.null(covar_df)) {
  cat(sprintf("[CHECK] covar non-finite count: %d\n", sum(!is.finite(as.matrix(covar_df)))))
}

############################################################
## 5) 실행 예시
############################################################

## A) 태깅 R^2 (컷별 비교 + 간단 CI)
r2_ld005 <- grm_tagging_R2_ci(prefix_common, prefix_ld005, K = 5e5, B = 10, seed = 42)
r2_ld010 <- grm_tagging_R2_ci(prefix_common, prefix_ld010, K = 5e5, B = 10, seed = 42)
print(r2_ld005)
print(r2_ld010)
cat(sprintf("\n[Tagging] LD<=0.05: R2_mean=%.4f (unique=%.4f)\n", r2_ld005$R2_mean, 1 - r2_ld005$R2_mean))
cat(sprintf("[Tagging] LD<=0.10: R2_mean=%.4f (unique=%.4f)\n\n", r2_ld010$R2_mean, 1 - r2_ld010$R2_mean))

## B+C) 정직교 rare⊥ + permutation 대조 (컷별)
res_ortho_005 <- he_ortho_and_perm(prefix_common, prefix_ld005, y_vec, covar_df, K = 1e6, seed = 42)
res_ortho_010 <- he_ortho_and_perm(prefix_common, prefix_ld010, y_vec, covar_df, K = 1e6, seed = 43)

print(res_ortho_005)
print(res_ortho_010)

cat(sprintf("\n[Ortho] LD<=0.05: corr(Gc,Gr⊥)=%.3g, h2_rare⊥=%.3g, p=%.3g, perm p=%.3g\n",
            res_ortho_005$corr_Gc_Gr_perp, res_ortho_005$h2_rare_perp,
            res_ortho_005$p_Gr_perp, res_ortho_005$p_perm))
cat(sprintf("[Ortho] LD<=0.10: corr(Gc,Gr⊥)=%.3g, h2_rare⊥=%.3g, p=%.3g, perm p=%.3g\n",
            res_ortho_010$corr_Gc_Gr_perp, res_ortho_010$h2_rare_perp,
            res_ortho_010$p_Gr_perp, res_ortho_010$p_perm))
