## ===========================================
## GRM off-diagonal correlation (low-memory, fast-combine)
## - GCTA .grm.bin (lower-tri, double) 스트리밍 행단위 읽기
## - 메모리 O(1), 원소별 루프 제거(배치 결합)로 속도 개선
## - n이 아주 클 때는 sampled 버전 사용 권장
## ===========================================

.read_grm_ids <- function(prefix) {
  read.table(paste0(prefix, ".grm.id"), stringsAsFactors = FALSE)
}

.check_same_ids <- function(id1, id2) {
  if (nrow(id1) != nrow(id2)) return(FALSE)
  all(id1[[1]] == id2[[1]] & id1[[2]] == id2[[2]])
}

## 한 배치(x, y 벡터)에 대한 통계 요약(평균, M2, 공분산 합) 계산
.batch_stats <- function(x, y, skip_nonfinite = TRUE) {
  if (skip_nonfinite) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
  }
  m <- length(x)
  if (m == 0) {
    return(list(n = 0.0, meanx = 0.0, meany = 0.0, M2x = 0.0, M2y = 0.0, Cxy = 0.0))
  }
  sx  <- sum(x);  sy  <- sum(y)
  mx  <- sx / m;  my  <- sy / m
  dx  <- x - mx;  dy  <- y - my
  M2x <- sum(dx * dx)
  M2y <- sum(dy * dy)
  Cxy <- sum(dx * dy)
  list(n = as.numeric(m), meanx = mx, meany = my, M2x = M2x, M2y = M2y, Cxy = Cxy)
}

## 두 통계 요약을 Chan/Welford 방식으로 결합
.combine_stats <- function(acc, b) {
  ## acc: list(n, meanx, meany, M2x, M2y, Cxy)
  ## b  : list(n, meanx, meany, M2x, M2y, Cxy)
  n1 <- acc$n; n2 <- b$n
  if (n2 == 0) return(acc)
  if (n1 == 0) return(b)

  n  <- n1 + n2
  dx <- b$meanx - acc$meanx
  dy <- b$meany - acc$meany

  meanx <- acc$meanx + dx * (n2 / n)
  meany <- acc$meany + dy * (n2 / n)

  M2x <- acc$M2x + b$M2x + (dx * dx) * (n1 * n2 / n)
  M2y <- acc$M2y + b$M2y + (dy * dy) * (n1 * n2 / n)
  Cxy <- acc$Cxy + b$Cxy + (dx * dy) * (n1 * n2 / n)

  list(n = n, meanx = meanx, meany = meany, M2x = M2x, M2y = M2y, Cxy = Cxy)
}

## 온라인 피어슨 상관(오프대각만) – 배치 결합 버전
grm_offdiag_cor <- function(prefix1,
                            prefix2 = NULL,
                            check_ids = TRUE,
                            endian = .Platform$endian,  # 보통 "little"
                            skip_nonfinite = TRUE,
                            verbose = TRUE) {
  id1 <- .read_grm_ids(prefix1)
  n   <- nrow(id1)

  if (is.null(prefix2)) {
    prefix2 <- prefix1
    id2 <- id1
  } else {
    id2 <- .read_grm_ids(prefix2)
    if (check_ids && !.check_same_ids(id1, id2)) {
      stop("IDs or order differ between GRMs (.grm.id mismatch).")
    }
  }

  con1 <- file(paste0(prefix1, ".grm.bin"), "rb")
  con2 <- file(paste0(prefix2, ".grm.bin"), "rb")
  on.exit({ close(con1); close(con2) }, add = TRUE)

  ## 누적 통계 (모두 double로 유지: 오버플로우 방지)
  acc <- list(n = 0.0, meanx = 0.0, meany = 0.0, M2x = 0.0, M2y = 0.0, Cxy = 0.0)
  offdiag_pairs <- 0.0

  if (verbose) cat(sprintf("[INFO] n=%d; streaming rows (batch combine)...\n", n))

  for (i in 1:n) {
    seg1 <- readBin(con1, what = "numeric", n = i, size = 8, endian = endian)
    seg2 <- readBin(con2, what = "numeric", n = i, size = 8, endian = endian)

    if (i > 1) {
      x <- seg1[1:(i-1)]  # off-diagonal
      y <- seg2[1:(i-1)]
      if (skip_nonfinite) {
        ok <- is.finite(x) & is.finite(y)
        x <- x[ok]; y <- y[ok]
      }
      b <- .batch_stats(x, y, skip_nonfinite = FALSE)  # 이미 ok 필터함
      acc <- .combine_stats(acc, b)
      offdiag_pairs <- offdiag_pairs + b$n
    }

    if (verbose && (i %% max(1L, floor(n/10))) == 0L) {
      cat(sprintf("  ... %d / %d rows\n", i, n))
    }
  }

  corr <- if (is.finite(acc$M2x) && is.finite(acc$M2y) && acc$M2x > 0 && acc$M2y > 0) {
    acc$Cxy / sqrt(acc$M2x * acc$M2y)
  } else NA_real_

  list(n = n, offdiag_pairs = as.numeric(offdiag_pairs), correlation = corr)
}

## (옵션) 빠른 미리보기: off-diagonal 일부 샘플만 읽기
## - 균등 샘플링된 오프대각 K개만 읽어 상관 추정
## - 인덱싱 버그 수정 반영
grm_offdiag_cor_sampled <- function(prefix1,
                                    prefix2,
                                    K = 2e6,
                                    check_ids = TRUE,
                                    endian = .Platform$endian,  # 보통 "little"
                                    skip_nonfinite = TRUE,
                                    verbose = TRUE) {
  id1 <- .read_grm_ids(prefix1)
  id2 <- .read_grm_ids(prefix2)
  if (check_ids && !.check_same_ids(id1, id2)) {
    stop("IDs or order differ between GRMs (.grm.id mismatch).")
  }
  n <- nrow(id1)
  total_offdiag <- as.numeric(n) * (as.numeric(n) - 1) / 2
  K <- as.integer(min(K, total_offdiag))

  if (verbose) {
    cat(sprintf("[INFO] n=%d, total off-diag=%s, sampling K=%s\n",
                n, format(total_offdiag, big.mark=","), format(K, big.mark=",")))
  }

  ## 1) 오프대각(대각 제외)에서 0-기반 선형 인덱스 샘플
  samp <- sort(sample.int(as.integer(total_offdiag), K, replace = FALSE) - 1L)

  ## 2) 0-기반 오프대각 인덱스를 (i, j)로 역변환: i=1..n, j=1..(i-1)
  i_vals <- integer(K); j_vals <- integer(K)
  cum <- 0L
  i <- 1L
  ptr <- 1L
  while (ptr <= K) {
    need <- samp[ptr]
    while (cum + (i - 1L) <= need) {
      cum <- cum + (i - 1L)
      i <- i + 1L
    }
    j0 <- need - cum
    i_vals[ptr] <- i
    j_vals[ptr] <- j0 + 1L
    ptr <- ptr + 1L
  }

  ## 3) 하삼각(대각 포함) 1-기반 선형 인덱스로 변환: L = i*(i-1)/2 + j
  lin_idx <- i_vals*(i_vals - 1L)/2L + j_vals

  con1 <- file(paste0(prefix1, ".grm.bin"), "rb")
  con2 <- file(paste0(prefix2, ".grm.bin"), "rb")
  on.exit({ close(con1); close(con2) }, add = TRUE)

  bytes_per <- 8L
  read_one <- function(con, index) {
    seek(con, where = (index - 1L) * bytes_per, origin = "start")
    readBin(con, what = "numeric", n = 1L, size = bytes_per, endian = endian)
  }

  if (verbose) cat("[INFO] Sampling values...\n")
  x <- numeric(K); y <- numeric(K)
  for (t in seq_len(K)) {
    x[t] <- read_one(con1, lin_idx[t])
    y[t] <- read_one(con2, lin_idx[t])
    if (verbose && (t %% max(1L, floor(K/10))) == 0L) {
      cat(sprintf("  ... %d / %d sampled\n", t, K))
    }
  }

  if (skip_nonfinite) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
  }

  list(n = n, sampled_offdiag = length(x), correlation = if (length(x) > 1) cor(x, y) else NA_real_)
}


fisher_ci_p <- function(r, m, conf = 0.95) {
  if (!is.finite(r) || m <= 3) return(list(se = NA_real_, lwr = NA_real_, upr = NA_real_, p = NA_real_))
  z <- atanh(r)
  se.z <- 1 / sqrt(m - 3)
  alpha <- 1 - conf
  z.l <- z - qnorm(1 - alpha/2) * se.z
  z.u <- z + qnorm(1 - alpha/2) * se.z
  ci <- tanh(c(z.l, z.u))
  # p-value for H0: rho=0 (equivalently t-test)
  tval <- r * sqrt((m - 2) / (1 - r^2))
  pval <- 2 * pt(-abs(tval), df = m - 2)
  list(se = se.z * (1 - r^2),   # delta method on r (optional; often report SE on z-scale)
       se_z = se.z,             # SE on z-scale(권장 표기)
       lwr = ci[1], upr = ci[2],
       p = pval)
}

      
## ========== 사용 예시 ==========
## 0) self-check (같은 파일로 비교: corr ≈ 1)
## grm_offdiag_cor("v6.ABC.RARE.Aset", endian = "little")

## 1) 서로 다른 GRM 비교 (ID 순서 동일해야 함)
## grm_offdiag_cor("v6.ABC.RARE.Aset", "v6.ABC.COMMON.Aset", endian = "little")

## 2) 빠른 미리보기(샘플링)
## grm_offdiag_cor_sampled("v6.ABC.RARE.Aset", "v6.ABC.COMMON.Aset", K = 2e6, endian = "little")

#res <- grm_offdiag_cor("v6.ABC.RARE.Aset", "v6.ABC.COMMON.Aset", endian="little")
#m <- res$offdiag_pairs
#r <- res$correlation
#out <- fisher_ci_p(r, m)
#out
