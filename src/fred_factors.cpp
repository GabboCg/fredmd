// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// minindc_cpp: for each column, return 1-based index of the minimum value
// [[Rcpp::export]]
uvec minindc_cpp(mat X) {
  int n_cols = X.n_cols;
  uvec pos(n_cols);
  for (int i = 0; i < n_cols; i++) {
    pos(i) = X.col(i).index_min() + 1;  // 1-based
  }
  return pos;
}

// pc2_cpp: PCA via economy SVD of X'X
// Returns: List{chat, fhat, lambda, diagmat(s)}
// [[Rcpp::export]]
List pc2_cpp(mat X, int nfac) {
  int N = X.n_cols;
  mat U, V;
  vec s;
  svd_econ(U, s, V, X.t() * X);

  mat lambda = U.cols(0, nfac - 1) * std::sqrt((double)N);
  mat fhat   = X * lambda / N;
  mat chat   = fhat * lambda.t();
  mat ss     = diagmat(s);

  return List::create(chat, fhat, lambda, ss);
}

// transform_data_cpp: demean/standardize x2
// DEMEAN: 0=none, 1=demean, 2=demean+standardize, 3=recursive demean+standardize
// Returns: List{x22, mut, std}  — all T0×N matrices
// [[Rcpp::export]]
List transform_data_cpp(mat x2, int DEMEAN) {
  int T0 = x2.n_rows;
  int N  = x2.n_cols;
  mat mut, std_mat, x22;

  if (DEMEAN == 0) {
    mut     = zeros(T0, N);
    std_mat = ones(T0, N);
    x22     = x2;

  } else if (DEMEAN == 1) {
    rowvec col_means = mean(x2, 0);          // 1×N
    mut     = repmat(col_means, T0, 1);
    std_mat = ones(T0, N);
    x22     = x2 - mut;

  } else if (DEMEAN == 2) {
    rowvec col_means = mean(x2, 0);
    rowvec col_sds   = stddev(x2, 0, 0);    // N-1 normalization, column-wise
    mut     = repmat(col_means, T0, 1);
    std_mat = repmat(col_sds, T0, 1);
    x22     = (x2 - mut) / std_mat;

  } else {  // DEMEAN == 3
    mut = mat(T0, N);
    for (int t = 0; t < T0; t++) {
      mut.row(t) = mean(x2.rows(0, t), 0);  // recursive column means
    }
    rowvec col_sds = stddev(x2, 0, 0);
    std_mat = repmat(col_sds, T0, 1);
    x22 = (x2 - mut) / std_mat;
  }

  return List::create(x22, mut, std_mat);
}

// baing_cpp: Bai & Ng (2002) information criterion for factor selection
// Returns: List{ic1 (int), chat, Fhat, diagmat(eigval)}
// [[Rcpp::export]]
List baing_cpp(mat X, int kmax, int jj) {
  int T0  = X.n_rows;
  int N   = X.n_cols;
  double NT  = (double)(N * T0);
  double NT1 = (double)(N + T0);
  int GCT = std::min(N, T0);

  // Penalty vector CT (length kmax, 0-indexed; CT(i) = penalty for i+1 factors)
  rowvec ii  = linspace<rowvec>(1, kmax, kmax);
  rowvec CT(kmax);
  if (jj == 1) {
    CT = std::log(NT / NT1) * ii * NT1 / NT;
  } else if (jj == 2) {
    CT = (NT1 / NT) * std::log((double)std::min(N, T0)) * ii;
  } else {
    CT = ii * std::log((double)GCT) / GCT;
  }

  // SVD branch
  mat Fhat0, Lambda0;
  vec eigval;
  if (T0 < N) {
    mat U, V; vec s;
    svd_econ(U, s, V, X * X.t());
    eigval  = s;
    Fhat0   = std::sqrt((double)T0) * U;
    Lambda0 = X.t() * Fhat0 / T0;
  } else {
    mat U, V; vec s;
    svd_econ(U, s, V, X.t() * X);
    eigval  = s;
    Lambda0 = std::sqrt((double)N) * U;
    Fhat0   = X * Lambda0 / N;
  }

  // Information criterion loop
  vec Sigma(kmax + 1, fill::zeros);
  vec IC1(kmax + 1, fill::zeros);

  for (int i = kmax - 1; i >= 0; i--) {
    mat Fhat   = Fhat0.cols(0, i);
    mat lambda = Lambda0.cols(0, i);
    mat chat   = Fhat * lambda.t();
    mat ehat   = X - chat;
    Sigma(i)   = mean(sum(square(ehat) / T0, 0));  // mean of column sums
    IC1(i)     = std::log(Sigma(i)) + CT(i);
  }
  // Zero-factor baseline
  Sigma(kmax) = mean(sum(square(X) / T0, 0));
  IC1(kmax)   = std::log(Sigma(kmax));

  // Optimal factor count (1-based; 0 if zero-factor wins)
  int ic1 = (int)index_min(IC1) + 1;
  if (ic1 > kmax) ic1 = 0;

  // Return chat and Fhat using kmax factors
  mat Fhat_out   = Fhat0.cols(0, kmax - 1);
  mat Lambda_out = Lambda0.cols(0, kmax - 1);
  mat chat_out   = Fhat_out * Lambda_out.t();

  return List::create(ic1, chat_out, Fhat_out, diagmat(eigval));
}

// factors_em_cpp: EM algorithm for missing value imputation + factor extraction
// x: T0×N numeric matrix (date column already stripped by caller)
// Returns: List{ehat, Fhat, lamhat, ve2, x2}
// [[Rcpp::export]]
List factors_em_cpp(mat x, int kmax, int jj, int DEMEAN) {
  int maxit = 50;
  int T0 = x.n_rows;
  int N  = x.n_cols;
  double err = 999.0;
  int it = 0;
  int icstar;

  // Missing value indicator (1 = missing/NA → NaN in C++)
  mat x1(T0, N, fill::zeros);
  for (int t = 0; t < T0; t++)
    for (int j = 0; j < N; j++)
      if (!std::isfinite(x(t, j))) x1(t, j) = 1.0;

  // Column means ignoring NAs (NaN in C++)
  rowvec col_means(N, fill::zeros);
  for (int j = 0; j < N; j++) {
    vec col_j = x.col(j);
    uvec finite_idx = find_finite(col_j);
    if (finite_idx.n_elem > 0)
      col_means(j) = mean(col_j.elem(finite_idx));
  }

  // Initialize x2: replace NAs with column means
  mat x2 = x;
  for (int j = 0; j < N; j++)
    for (int t = 0; t < T0; t++)
      if (x1(t, j) == 1.0) x2(t, j) = col_means(j);

  // Initial transform + factor selection + PCA
  List td = transform_data_cpp(x2, DEMEAN);
  mat x3      = as<mat>(td[0]);
  mat mut     = as<mat>(td[1]);
  mat std_mat = as<mat>(td[2]);

  if (kmax != 99) {
    icstar = as<int>(baing_cpp(x3, kmax, jj)[0]);
    if (icstar <= 0) icstar = 1;  // guard: zero-factor solution → use 1 factor
  } else {
    icstar = 8;
  }

  List pc2_res = pc2_cpp(x3, icstar);
  mat chat  = as<mat>(pc2_res[0]);
  mat chat0 = chat;

  // EM loop
  while (err > 1e-6 && it < maxit) {
    it++;
    Rcpp::Rcout << "Iteration " << it << ": obj " << err
                << " IC " << icstar << "\n";

    // Hot inner loop: impute missing, keep observed
    for (int t = 0; t < T0; t++) {
      for (int j = 0; j < N; j++) {
        if (x1(t, j) == 1.0) {
          x2(t, j) = chat(t, j) * std_mat(t, j) + mut(t, j);
        } else {
          x2(t, j) = x(t, j);
        }
      }
    }

    td      = transform_data_cpp(x2, DEMEAN);
    x3      = as<mat>(td[0]);
    mut     = as<mat>(td[1]);
    std_mat = as<mat>(td[2]);

    if (kmax != 99) {
      icstar = as<int>(baing_cpp(x3, kmax, jj)[0]);
      if (icstar <= 0) icstar = 1;  // guard: zero-factor solution → use 1 factor
    } else {
      icstar = 8;
    }

    pc2_res = pc2_cpp(x3, icstar);
    chat = as<mat>(pc2_res[0]);

    mat diff_ = chat - chat0;
    err   = accu(pow(diff_, 2)) / accu(pow(chat0, 2));
    chat0 = chat;
  }

  if (it == maxit)
    Rcpp::warning("Maximum number of iterations reached in EM algorithm");

  mat ehat   = x - chat % std_mat - mut;  // % = element-wise multiply
  mat Fhat   = as<mat>(pc2_res[1]);
  mat lamhat = as<mat>(pc2_res[2]);
  mat ve2    = as<mat>(pc2_res[3]);

  return List::create(ehat, Fhat, lamhat, ve2, x2);
}
