library(optparse)
library(jpeg)
library(Rcpp)
sourceCpp('./subcalc.cpp')

# this is for H0: beta = 0
# run this after calc.R

imread <- function(x, LSB=TRUE) {
    img_float <- readJPEG(x, native=F)
    img_8bit <- img_float * 255
    storage.mode(img_8bit) <- "integer"
    if (LSB) {
        img_8bit %% 2
    } else {
        img_8bit
    }
}

calc2_plhd <- function(params, sdat) {
    alpha <- params[1]
    beta <- params[2]
    f2stcalc_plhd(alpha, beta, sdat$value, sdat$same_neib, sdat$diff_neib)
}

get_stat2 <- function(samp, a0) {
    sdat <- calc_neibSame(samp)
    s <- optim(c(a0, 0), calc2_plhd, control=list(fnscale=-1), sdat=sdat)
    # print(s)
    abs(s$par[2])
}

do_bootstrap2 <- function(data, N=10, p0=0.5) {
    a0 <- log(p0 / (1-p0))
    true_statistic <- get_stat2(data, a0)
    simulated_stat <- c()
    for (i in 1:N) {
        print(i)
        sim_data <- matrix(rbinom(prod(dim(data)), size=1, prob=p0), nrow=nrow(data))
        simulated_stat[i] <- get_stat2(sim_data, a0)
    }
    bpvalue <- (1 + sum(true_statistic < simulated_stat)) / (1 + length(simulated_stat))
    print(sprintf("true stat %f", true_statistic))
    print(simulated_stat)
    print(sprintf("p-value (True < simul) = %f", bpvalue))
    bpvalue
    res <- list(true_stat=true_statistic, simul_stat=simulated_stat, p_value=bpvalue)
    res
}

submain <- function(i=1, N=10) {
    files <- sort(Sys.glob("./*.Rdata"))[i]
    load(files)
    image <- imread(res$name)
    parts <- c(1,2,3)
    colors <- c("R", "G", "B")
    newres <- list()
    for (i in parts) {
            r <- res[[colors[i]]]
            # print(r)
            check1 <- r[[1]]
            p_value <- (1 + sum(check1$simul_stat > check1$true_stat)) / (1 + length(check1$simul_stat))
            if (p_value < 0.05) {
                subimg <- image[,,i]
                p_hat <- mean(subimg)
                check2 <- do_bootstrap2(subimg, N, p_hat) 
            } else {
                check2 <- NULL
            }
            newpart <- list(check1=check1, check2=check2)
            newres[[colors[i]]] <- newpart
    }
    newres$name <- res$name
    newres
}

main <- function() {
    option_list <- list(
        make_option(c("-i", "--image-id"), type="integer", default=1,
                    dest="image_id", help="image id"),
        make_option(c("-n", "--num-samples"), type="integer", default=10,
                    dest="num_samples", help="num of bootstrap samples")
    )
    parser <- OptionParser(option_list=option_list)
    args <- parse_args(parser)
    # print(args)
    res <- submain(i=args$image_id, N=args$num_samples)
    name <- gsub("\\.jpg", "", basename(res$name))
    save(res, file=sprintf("%s-%d.Rdata", name, args$num_samples))
}

main()
