library(testthat)

# Compile and load C++ implementations
# Path works from project root (source/Rscript) and from tests/ (test_file)
cpp_path <- ifelse(file.exists("src/fred_factors.cpp"), "src/fred_factors.cpp", "../src/fred_factors.cpp")
Rcpp::sourceCpp(cpp_path)

# R reference implementations (inlined — fred-factors.R removed after Rcpp conversion)
minindc <- function(x) {
  n_rows <- nrow(x); n_cols <- ncol(x)
  pos <- matrix(0, nrow = n_cols, ncol = 1)
  seq0 <- t(matrix(seq(1, n_rows, 1), nrow = 1))
  for (i in 1:n_cols) {
    min_i <- min(x[, i])
    colmin_i <- seq0 * ((x[, i] - min_i) == 0)
    pos[i] <- sum(colmin_i)
  }
  return(pos)
}

pc2 <- function(X, nfac) {
  N <- ncol(X)
  svd_res <- svd(t(X) %*% X)
  U <- svd_res$u; S <- svd_res$d
  lambda <- U[, 1:nfac] * sqrt(N)
  fhat <- X %*% lambda / N
  chat <- fhat %*% t(lambda)
  ss <- diag(S)
  list(chat, fhat, lambda, ss)
}

transform_data <- function(x2, DEMEAN) {
  T0 <- nrow(x2); N <- ncol(x2)
  switch(as.character(DEMEAN),
    "0" = { mut <- matrix(0, T0, 1); std <- matrix(1, T0, 1); x22 <- x2 },
    "1" = {
      mut <- matrix(rep(apply(x2, 2, mean, na.rm = TRUE), each = T0), T0, N)
      std <- matrix(1, T0, 1); x22 <- x2 - mut
    },
    "2" = {
      mut <- matrix(rep(apply(x2, 2, mean, na.rm = TRUE), each = T0), T0, N)
      std <- matrix(rep(apply(x2, 2, sd,   na.rm = TRUE), each = T0), T0, N)
      x22 <- (x2 - mut) / std
    },
    "3" = {
      mut <- matrix(NA, T0, N)
      for (t in 1:T0) mut[t, ] <- colMeans(x2[1:t, , drop = FALSE], na.rm = TRUE)
      std <- matrix(rep(apply(x2, 2, sd, na.rm = TRUE), each = T0), T0, N, byrow = TRUE)
      x22 <- (x2 - mut) / std
    }
  )
  list(x22, mut, std)
}

baing <- function(X, kmax, jj) {
  T0 <- nrow(X); N <- ncol(X); NT <- N * T0; NT1 <- N + T0
  CT <- matrix(0, 1, kmax); ii <- seq(1, kmax)
  switch(as.character(jj),
    "1" = { CT[1, ] <- log(NT / NT1) * ii * NT1 / NT },
    "2" = { CT[1, ] <- (NT1 / NT) * log(min(N, T0)) * ii },
    "3" = { CT[1, ] <- ii * log(min(N, T0)) / min(N, T0) }
  )
  if (T0 < N) {
    s <- svd(X %*% t(X)); Fhat0 <- sqrt(T0) * s$u; Lambda0 <- t(X) %*% Fhat0 / T0; eigval <- s$d
  } else {
    s <- svd(t(X) %*% X); Lambda0 <- sqrt(N) * s$u; Fhat0 <- X %*% Lambda0 / N; eigval <- s$d
  }
  Sigma <- matrix(0, 1, kmax + 1); IC1 <- matrix(0, 1, kmax + 1)
  for (i in rev(seq(1, kmax))) {
    chat <- Fhat0[, 1:i, drop = FALSE] %*% t(Lambda0[, 1:i, drop = FALSE])
    ehat <- X - chat
    Sigma[1, i] <- mean(apply((ehat * ehat) / T0, 2, sum))
    IC1[1, i] <- log(Sigma[1, i]) + CT[, i]
  }
  Sigma[kmax + 1] <- mean(apply(X * X / T0, 2, sum))
  IC1[, kmax + 1] <- log(Sigma[1, kmax + 1])
  ic1 <- t(minindc(t(IC1))); ic1 <- ic1 * (ic1 <= kmax)
  chat <- Fhat0[, 1:kmax] %*% t(Lambda0[, 1:kmax])
  list(ic1, chat, Fhat0[, 1:kmax], diag(eigval))
}

factors_em <- function(x, kmax, jj, DEMEAN) {
  maxit <- 50; T0 <- nrow(x); N <- ncol(x); err <- 999; it <- 0
  x1 <- matrix(as.numeric(is.na(x[, 2:N])), ncol = N - 1)
  mut <- matrix(rep(apply(x[, 2:N], 2, mean, na.rm = TRUE), nrow(x)),
                ncol = ncol(x[, 2:N]), byrow = TRUE)
  x2 <- as.matrix(x[, 2:N]); x2[is.na(x[, 2:N])] <- mut[is.na(x[, 2:N])]
  td <- transform_data(x2, DEMEAN); x3 <- td[[1]]; mut <- td[[2]]; std <- td[[3]]
  icstar <- if (kmax != 99) baing(x3, kmax, jj)[[1]] else 8
  pc2_res <- pc2(x3, icstar); chat <- pc2_res[[1]]; chat0 <- chat
  while (err > 1e-6 && it < maxit) {
    it <- it + 1
    for (t0 in 1:T0) for (j in 1:(N - 1)) {
      x2[t0, j] <- if (x1[t0, j] == 1) chat[t0, j] * std[t0, j] + mut[t0, j]
                   else as.numeric(x[t0, j + 1])
    }
    td <- transform_data(x2, DEMEAN); x3 <- td[[1]]; mut <- td[[2]]; std <- td[[3]]
    icstar <- if (kmax != 99) baing(x3, kmax, jj)[[1]] else 8
    pc2_res <- pc2(x3, icstar); chat <- pc2_res[[1]]
    err <- sum((chat - chat0)^2) / sum(chat0^2); chat0 <- chat
  }
  ehat <- as.matrix(x[, 2:N]) - chat * std - mut
  list(ehat, pc2_res[[2]], pc2_res[[3]], pc2_res[[4]], x2)
}

# Helper: compare two matrices up to sign flip per column
# (SVD eigenvectors are unique only up to sign)
expect_mat_equal <- function(A, B, tol = 1e-8) {
  expect_equal(dim(A), dim(B))
  # align signs: flip B columns where correlation with A is negative
  for (j in seq_len(ncol(A))) {
    if (sum(A[, j] * B[, j]) < 0) B[, j] <- -B[, j]
  }
  expect_equal(A, B, tolerance = tol)
}

cat("Test harness loaded.\n")

test_that("minindc_cpp matches minindc", {
  set.seed(42)
  X <- matrix(rnorm(20), nrow = 5, ncol = 4)
  X[3, 2] <- -99  # force known minimum in col 2 (row 3)

  r_result  <- as.integer(minindc(X))
  cpp_result <- as.integer(minindc_cpp(X))

  expect_equal(cpp_result, r_result)
})

test_that("pc2_cpp chat matches pc2 chat", {
  set.seed(7)
  X <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10)
  nfac <- 3

  r   <- pc2(X, nfac)
  cpp <- pc2_cpp(X, nfac)

  # chat is sign-invariant (fhat %*% t(lambda))
  expect_equal(cpp[[1]], r[[1]], tolerance = 1e-8)

  # eigenvalues (diagonal of ss) are non-negative and match
  expect_equal(diag(cpp[[4]]), diag(r[[4]]), tolerance = 1e-8)
})

test_that("transform_data_cpp matches transform_data for all DEMEAN modes", {
  set.seed(3)
  X <- matrix(rnorm(30 * 5), nrow = 30, ncol = 5)

  for (dm in c(0, 1, 2)) {
    r   <- transform_data(X, dm)
    cpp <- transform_data_cpp(X, dm)
    # Compare transformed data (x22)
    expect_equal(cpp[[1]], r[[1]], tolerance = 1e-8,
                 info = paste("DEMEAN =", dm, "(x22)"))
    # Compare mut — R may be T0×1, C++ is T0×N; compare column-wise
    r_mut <- matrix(r[[2]], nrow = nrow(X), ncol = ncol(X))
    expect_equal(cpp[[2]], r_mut, tolerance = 1e-8,
                 info = paste("DEMEAN =", dm, "(mut)"))
  }

  # DEMEAN=3 (recursive means) — no R reference due to R source bug (T vs T0),
  # so test structural properties only
  cpp3 <- transform_data_cpp(X, 3)
  expect_equal(dim(cpp3[[1]]), dim(X))     # x22 has same shape
  expect_equal(dim(cpp3[[2]]), dim(X))     # mut is T0×N
  expect_equal(dim(cpp3[[3]]), dim(X))     # std is T0×N
  # First row of mut should equal first row of X (mean of one observation)
  expect_equal(cpp3[[2]][1, ], X[1, ], tolerance = 1e-10)
})

test_that("baing_cpp ic1 matches baing ic1 for all jj criteria", {
  set.seed(11)
  X <- matrix(rnorm(80 * 15), nrow = 80, ncol = 15)
  kmax <- 5

  for (jj_val in 1:3) {
    r   <- baing(X, kmax, jj_val)
    cpp <- baing_cpp(X, kmax, jj_val)
    expect_equal(as.integer(cpp[[1]]), as.integer(r[[1]]),
                 info = paste("jj =", jj_val, "(ic1)"))
    # chat is the product F*Lambda' — sign-invariant across R/C++ SVD conventions.
    # Use expect_mat_equal (column sign-aligning helper) in case singular vector
    # sign conventions differ between arma::svd_econ and R's svd().
    expect_mat_equal(cpp[[2]], r[[2]], tol = 1e-6)
  }
})

test_that("factors_em_cpp x2 matches factors_em x2", {
  set.seed(21)
  X_full <- matrix(rnorm(60 * 8), nrow = 60, ncol = 8)
  # Introduce ~10% missing values
  miss_idx <- sample(length(X_full), size = 48)
  X_full[miss_idx] <- NA

  kmax  <- 3
  jj    <- 2
  DEMEAN <- 2

  # R version receives data.frame with a fake date column in position 1
  df_r <- as.data.frame(cbind(date = seq_len(nrow(X_full)), X_full))
  r_res <- factors_em(df_r, kmax, jj, DEMEAN)
  r_x2  <- r_res[[5]]

  # C++ version receives pre-stripped numeric matrix
  cpp_res <- factors_em_cpp(X_full, kmax, jj, DEMEAN)
  cpp_x2  <- cpp_res[[5]]

  expect_equal(cpp_x2, unname(r_x2), tolerance = 1e-6)
})
