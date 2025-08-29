## ============================================================
## Off-diagonal pair sampler (64-bit safe, low RAM) + R^2 only
## - 두 GRM에서 같은 (i,j) 오프대각 값을 샘플링해서 lm(y ~ x)의 R^2 계산
## - GRM 전체를 메모리에 올리지 않음 (파일에서 필요한 쌍만 읽음)
## - total_offdiag > 2^31-1 환경에서도 동작
## ============================================================

grm_offdiag_sample_pairs <- function(prefix_x,
                                     prefix_y,
                                     K = 5e5,                 # 샘플링할 (i,j) 쌍 개수
                                     endian = .Platform$endian,
                                     check_ids = TRUE,
                                     skip_nonfinite = TRUE,
                                     seed = 1,
                                     verbose = TRUE) {
  ## 1) ID 읽고 동일성 확인
  id_x <- read.table(paste0(prefix_x, ".grm.id"), stringsAsFactors = FALSE)
  id_y <- read.table(paste0(prefix_y, ".grm.id"), stringsAsFactors = FALSE)
  if (nrow(id_x) != nrow(id_y)) stop("ID length mismatch")
  if (check_ids) {
    same <- all(id_x[[1]] == id_y[[1]] & id_x[[2]] == id_y[[2]])
    if (!same) stop(".grm.id mismatch in IDs/order between X and Y")
  }
  n  <- nrow(id_x)
  total_offdiag <- as.numeric(n) * (as.numeric(n) - 1) / 2  # double
  K  <- as.integer(min(K, total_offdiag))

  if (verbose) cat(sprintf("[INFO] n=%d; total off-diag=%s; sampling K=%s\n",
                           n, format(total_offdiag, big.mark=","), format(K, big.mark=",")))
  set.seed(seed)

  ## 2) 0-based 오프대각 인덱스 샘플링 (sample.int 한계 회피)
  want <- K
  buf  <- numeric(0)   # double 벡터
  while (length(buf) < K) {
    batch <- as.integer(max(want * 1.3, 1e6))     # 넉넉히 뽑아 unique
    cand  <- floor(runif(batch, min = 0, max = total_offdiag))  # [0, total_offdiag)
    buf   <- unique(c(buf, cand))
    want  <- K - length(buf)
  }
  samp0 <- sort(buf[seq_len(K)])  # 0-based 인덱스 K개

  ## 3) 0-based k -> (i,j) 역변환 (double 산술: overflow 방지)
  i_floor <- floor((1 + sqrt(1 + 8 * samp0)) / 2)   # numeric(double)
  row_i   <- i_floor + 1                             # 1..n (double)
  Tij     <- i_floor * (i_floor - 1) / 2            # T_i (double)
  j_col   <- as.numeric(samp0 - Tij + 1)            # 1..(row_i-1), double

  ## 4) 하삼각(대각 포함) 1-based 선형 인덱스: L = i*(i-1)/2 + j (double 산술)
  lin_idx <- as.numeric(row_i) * (as.numeric(row_i) - 1) / 2 + as.numeric(j_col)
  max_L   <- as.numeric(n) * (as.numeric(n) + 1) / 2
  if (!all(is.finite(lin_idx))) stop("lin_idx contains non-finite values (overflow?)")
  if (any(lin_idx < 1) || any(lin_idx > max_L)) stop("lin_idx out of range")

  ## 5) 같은 위치에서 X,Y 값을 파일에서 직접 읽기 (메모리 절약)
  conx <- file(paste0(prefix_x, ".grm.bin"), "rb")
  cony <- file(paste0(prefix_y, ".grm.bin"), "rb")
  on.exit({ close(conx); close(cony) }, add = TRUE)

  bytes_per <- 8
  read_one <- function(con, index) {
    where <- (as.numeric(index) - 1) * bytes_per   # double byte offset
    seek(con, where = where, origin = "start")
    val <- readBin(con, what = "numeric", n = 1L, size = bytes_per, endian = endian)
    if (length(val) == 0L) NA_real_ else val
  }

  if (verbose) cat("[INFO] Reading sampled pairs...\n")
  Kint <- length(lin_idx)
  x <- numeric(Kint); y <- numeric(Kint)
  for (t in seq_len(Kint)) {
    x[t] <- read_one(conx, lin_idx[t])
    y[t] <- read_one(cony, lin_idx[t])
    if (verbose && (t %% max(1L, floor(Kint/10))) == 0L) {
      cat(sprintf("  ... %d / %d sampled\n", t, Kint))
    }
  }

  if (skip_nonfinite) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
  }

  list(n = n, K = length(x), x = x, y = y)
}

## ============================================================
## R^2 만 계산 (rare ~ common)
## ============================================================
grm_R2 <- function(prefix_common,
                   prefix_rare,
                   K = 5e5,        # 빠르게 확인: 5e5, 최종 보고: 2e6 권장
                   seed = 1,
                   endian = "little",
                   verbose = TRUE) {
  samp <- grm_offdiag_sample_pairs(prefix_x = prefix_common,
                                   prefix_y = prefix_rare,
                                   K = K, seed = seed,
                                   endian = endian,
                                   verbose = verbose)
  if (length(samp$x) < 10L) stop("Too few sampled pairs after filtering")
  fit <- lm(samp$y ~ samp$x)
  summary(fit)$r.squared
}


R2_005 <- grm_R2(
  "/data/BiO/hae/WGS_rev/0_common/v2.common.grm",
  "/data/BiO/hae/WGS_rev/Z_LD005/v4.grm.LD005",
  K = 5e5, seed = 42, endian = "little"
)

R2_010 <- grm_R2(
  "/data/BiO/hae/WGS_rev/0_common/v2.common.grm",
  "/data/BiO/hae/WGS_rev/2_LD01/v4.grm.LD01",
  K = 5e5, seed = 42, endian = "little"
)

cat("Cutoff 0.05 R^2 =", R2_005, "\n")
cat("Cutoff 0.10 R^2 =", R2_010, "\n")
