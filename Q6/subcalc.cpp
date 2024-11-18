#include <RcppArmadillo.h>

#include <iostream>

#define USING_OPENMP 1

#if USING_OPENMP
// [[Rcpp::plugins(openmp)]]

#pragma omp declare reduction(+ : arma::mat : omp_out += omp_in) \
    initializer(omp_priv = omp_orig)
#endif

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List find_svd(arma::mat X) {
    Rcpp::List answer;
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U, s, V, X);
    answer = Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("s") = s,
                                Rcpp::Named("V") = V);
    return answer;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat calc_reduced(arma::mat X, arma::uword k, bool resid) {
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U, s, V, X);
    arma::mat U_k = U.submat(arma::span::all, arma::span(0, k - 1));
    arma::mat s_k = arma::diagmat(s(arma::span(0, k - 1)));
    arma::mat V_k = V.submat(arma::span::all, arma::span(0, k - 1)).t();
    arma::mat result;
    if (resid) {
        result = X - ((U_k * s_k) * V_k);  // E_k = X - X_k
    } else {
        result = ((U_k * s_k) * V_k);  // X_k
    }
    return result;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat calc_sigma2(arma::mat E) {
    arma::mat ET = E.t();
    arma::mat sigma = arma::cov(ET, 1);
    return sigma;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat calc_sigma(arma::mat E) {
    arma::uword N = E.n_cols;
    arma::uword M = E.n_rows;
    arma::mat sigma(M, M);
    sigma.fill(0.0);

    arma::colvec mu = arma::mean(E, 1);
    // std::cout << E.n_rows << " " << E.n_cols << " " << mu.size() << std::endl;
    // this is calculating Sigma
#if USING_OPENMP
#pragma omp parallel for reduction(+ : sigma)
#endif
    for (arma::uword i = 0; i < N; i++) {
        arma::mat tmp = (E.col(i) - mu) * (E.col(i) - mu).t();
        tmp /= N;
        sigma += tmp;
    }
    return sigma;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat calc_gamma(arma::mat E) {
    arma::uword N = E.n_rows;
    arma::uword M = E.n_cols;
    arma::mat gamma(M, M);
    gamma.fill(0.0);

    // this is calculating Gamma
    arma::rowvec mu = arma::mean(E, 0);
    // std::cout << E.n_rows << " " << E.n_cols << " " << mu.size() << std::endl;
#if USING_OPENMP
#pragma omp parallel for reduction(+ : gamma)
#endif
    for (arma::uword i = 0; i < N; i++) {
        arma::mat tmp = (E.row(i) - mu).t() * (E.row(i) - mu);
        tmp /= N;
        gamma += tmp;
    }
    return gamma;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat calc_gamma2(arma::mat E) {
    arma::mat gamma = arma::cov(E, 1);
    return gamma;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat calc_invchol(arma::mat m) {
    return arma::inv(arma::chol(m, "lower"));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double calc_teststat(arma::mat E) {
    arma::mat z = E.t() * E;
    return arma::trace(z);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double getfast_statistic(arma::mat X, arma::uword k) {
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U, s, V, X);
    arma::mat U_k = U.submat(arma::span::all, arma::span(0, k - 1));
    arma::mat s_k = arma::diagmat(s(arma::span(0, k - 1)));
    arma::mat V_k = V.submat(arma::span::all, arma::span(0, k - 1)).t();
    arma::mat E_k;
    E_k = X - ((U_k * s_k) * V_k);  // E_k = X - X_k
    arma::mat z = E_k.t() * E_k;
    return arma::trace(z);
}
