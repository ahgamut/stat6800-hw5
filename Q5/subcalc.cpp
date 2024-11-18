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
Rcpp::List calc_neibSame(arma::umat vals) {
    arma::uword M = vals.n_rows;
    arma::uword N = vals.n_cols;

    arma::mat dvals(M, N);
    arma::mat numpos(M, N);
    arma::mat same_neib(M, N);
    arma::mat diff_neib(M, N);
#if USING_OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (arma::uword i = 0; i < M; i++) {
        for (arma::uword j = 0; j < N; j++) {
            // std::cout << i << " " << j << std::endl;
            unsigned int same = 0;
            unsigned int diff = 0;
            unsigned int total = 0;
            bool toprow = i == 0;
            bool botrow = i == M - 1;
            bool leftcol = j == 0;
            bool rightcol = j == N - 1;
            same += toprow ? 0 : (vals(i, j) == vals(i - 1, j));
            same += botrow ? 0 : (vals(i, j) == vals(i + 1, j));
            same += leftcol ? 0 : (vals(i, j) == vals(i, j - 1));
            same += rightcol ? 0 : (vals(i, j) == vals(i, j + 1));
            total = (!toprow + !botrow + !leftcol + !rightcol);
            diff = total - same;
            //
            dvals(i, j) = vals(i, j);
            same_neib(i, j) = same;
            numpos(i, j) = total;
            diff_neib(i, j) = diff;
        }
    }

    auto answer = Rcpp::List::create(Rcpp::Named("value") = dvals,
                                     Rcpp::Named("pos_neib") = numpos,
                                     Rcpp::Named("same_neib") = same_neib,
                                     Rcpp::Named("diff_neib") = diff_neib);
    return answer;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double fastcalc_plhd(double alpha, double beta, const arma::mat vals,
                     const arma::mat same_neib, const arma::mat diff_neib) {
    /* this is slow because of the memory allocations */
    arma::mat num = arma::exp(alpha + beta * same_neib);
    arma::mat den = arma::exp(beta * diff_neib);
    arma::mat p = num / (num + den);
    arma::mat q = 1 - p;
    arma::mat plhd = (vals % arma::log(p)) + ((1 - vals) % arma::log(q));
    return arma::accu(plhd);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double f2stcalc_plhd(const double alpha, const double beta,
                     const arma::mat vals, const arma::mat same_neib,
                     const arma::mat diff_neib) {
    arma::uword M = vals.n_rows;
    arma::uword N = vals.n_cols;
    // arma::mat plhd(M, N);
    double plsum = 0.0;
#if USING_OPENMP
#pragma omp parallel for collapse(2) shared(vals, same_neib, diff_neib) reduction(+:plsum)
#endif
    for (arma::uword i = 0; i < M; i++) {
        for (arma::uword j = 0; j < N; j++) {
            double num = std::exp(alpha + beta * same_neib(i, j));
            double den = num + std::exp(beta * diff_neib(i, j));
            double p = num / den;
            double q = 1 - p;
            double plval = vals(i, j) * std::log(p) + (1 - vals(i, j)) * std::log(q);
            plsum += plval;
            // plhd(i, j) = vals(i, j) * std::log(p) + (1 - vals(i, j)) * std::log(q);
        }
    }
    return plsum; // arma::accu(plhd);
}
