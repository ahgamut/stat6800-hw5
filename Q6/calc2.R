library(optparse)
library(Rcpp)
sourceCpp('./subcalc.cpp')

# Z_star selected as bootstrap from Z_k

load_data <- function() {
    X <- matrix(scan("./jsm.dat"), nrow=2107, byrow=T)
    X# [1:20, 1:100]
}

get_statistic <- function(X, k) {
    X_k <- calc_reduced(X, k, F)
    E_k <- X - X_k
    calc_teststat(E_k)
}

runner <- function(k, bsim = 10) {
    X <- load_data()
    true_statistic <- getfast_statistic(X, k)
    # print(true_statistic)

    # setup for simulation
    X_k <- calc_reduced(X, k, F)
    E_k <- X - X_k
    lambda1 <- 1e-12
    lambda2 <- 1e-10

    sigma <- calc_sigma2(E_k)
    sigma <- sigma + lambda1 * diag(nrow(sigma))
    A <- t(chol(sigma))
    # print(dim(A))

    gamma <- calc_gamma2(E_k)
    gamma <- gamma + lambda2 * diag(nrow(gamma))
    B <- t(chol(gamma))
    # print(dim(B))

    Ainv <- solve(A)
    Binv <- solve(B)

    Z <- c(Ainv %*% E_k %*% Binv)

    simulated_stat <- c()
    for (i in 1:bsim) {
        print(i)
        Z_star <- matrix(sample(Z, size=length(Z), replace=T), nrow=nrow(sigma))
        X_star <- X_k + A %*% (Z_star %*% t(B))
        simulated_stat[i] <- getfast_statistic(X_star, k)
    }
    bpvalue <- mean(true_statistic < simulated_stat)

    print(simulated_stat)
    print(true_statistic)
    print(bpvalue)
    res <- list()
    res[["k"]] <- k
    res[["B"]] <- bsim
    res[["true_stat"]] <- true_statistic
    res[["simul_stats"]] <- simulated_stat
    res[["p_value"]] <- bpvalue
    res
}

submain <- function(k, b) {
    info <- runner(k, b)
    save(info, file=sprintf("q6-k%d-b%d.Rdata", k ,b))
}

main <- function() {
    option_list <- list(
        make_option(c("-k", "--num-svd"), type="integer", default=1,
                    dest="num_svd", help="number of singular values to use"),
        make_option(c("-n", "--num-samples"), type="integer", default=10,
                    dest="num_samples", help="num of bootstrap samples")
    )
    parser <- OptionParser(option_list=option_list)
    args <- parse_args(parser)
    submain(k=args$num_svd, b=args$num_samples)
}

main()
