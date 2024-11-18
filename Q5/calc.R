library(optparse)
library(jpeg)
library(Rcpp)
sourceCpp('./subcalc.cpp')

# this is for H0: (alpha, beta) = (0, 0)
# run this first

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

load_data <- function() {
    files <- sort(Sys.glob("UnitedColorsOfFall/*.jpg"))
    images <- lapply(files, imread)
    images
}

load_single <- function(i=1) {
    files <- sort(Sys.glob("UnitedColorsOfFall/*.jpg"))[i]
    images <- imread(files)
    res <- list(name=files, img=images)
    res
}

calc_loglike <- function(params, sdat) {
    alpha <- params[1]
    beta <- params[2]
    if (alpha < 0) {
        return(-Inf)
    }
    same_neib <- sdat$same_neib
    diff_neib <- sdat$pos_neib - sdat$same_neib
    xij <- sdat$value

    num <- exp(alpha + beta * sdat$same_neib) 
    den <- exp(beta * diff_neib) + num

    pij <- num / den
    qij <- 1 - pij

    plhd <- xij * log(pij) + (1 - xij) * log(qij)
    res <- sum(plhd)
    # print(c(alpha, beta, res))
    res
}

calc2_plhd <- function(params, sdat) {
    alpha <- params[1]
    beta <- params[2]
    if (alpha < 0) {
        return(-Inf)
    }
    f2stcalc_plhd(alpha, beta, sdat$value, sdat$same_neib, sdat$diff_neib)
}

get_stat <- function(samp) {
    sdat <- calc_neibSame(samp)
    alpha <- 0
    beta <- 0
    s <- optim(c(0, 0), calc2_plhd, control=list(fnscale=-1), sdat=sdat)
    # s <- optim(c(0, 0), calc_loglike, control=list(fnscale=-1), sdat=sdat)
    s$value
}

check_speeds <- function() {
    x <- load_single(7)
    sdat <- calc_neibSame(x[,,3])
    alpha <- 1.2
    beta <- 0.7

    f1_times <- sapply(1:100, function(i) {
        a <- Sys.time()
        p <- fastcalc_plhd(alpha, beta, sdat$value, sdat$same_neib, sdat$diff_neib)
        b <- Sys.time()
        b - a
    })
    print(sprintf("f1 %f", mean(f1_times)))
    f2_times <- sapply(1:100, function(i) {
        a <- Sys.time()
        p <- f2stcalc_plhd(alpha, beta, sdat$value, sdat$same_neib, sdat$diff_neib)
        b <- Sys.time()
        b - a
    })
    print(sprintf("f2 %f", mean(f2_times)))
}

do_bootstrap <- function(data, N=10, p0=0.5) {
    true_statistic <- get_stat(data)
    simulated_stat <- c()
    for (i in 1:N) {
        print(i)
        sim_data <- matrix(rbinom(prod(dim(data)), size=1, prob=p0), nrow=nrow(data))
        simulated_stat[i] <- get_stat(sim_data)
    }
    bpvalue <- mean(true_statistic > simulated_stat)
    print(sprintf("true stat %f", true_statistic))
    print(simulated_stat)
    print(sprintf("p-value (True > simul) = %f", bpvalue))
    bpvalue
    res <- list(true_stat=true_statistic, simul_stat=simulated_stat, p_value=bpvalue)
    res
}

run0 <- function(data, N) {
    check1 <- do_bootstrap(data, N, p0=0.5)
    check2 <- NULL
    if (check1$p_value < 0.05) {
        print(sprintf("check1 has p-value of %f", check1$p_value))
        p_hat <- mean(data)
        alpha_hat <- log(p_hat / (1-p_hat))
        print(sprintf("testing for beta=0, alpha = %f (i.e p = %f)", alpha_hat, p_hat))
        check2 <- do_bootstrap(data, N, p0=p_hat)
        print(sprintf("check2 has p-value of %f", check2$p_value))
    }
    res <- list(check1, check2)
    res
}

submain <- function(imnum, N) {
    res <- list()
    imdata <- load_single(imnum)
    name <- gsub("\\.jpg", "", basename(imdata$name))
    image <- imdata$img
    res[["name"]] <- imdata$name
    res[["R"]] <- run0(image[,,1], N)
    res[["G"]] <- run0(image[,,2], N)
    res[["B"]] <- run0(image[,,3], N)
    print(res)
    save(res, file=sprintf("%s-%d.Rdata", name, N))
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
    submain(imnum=args$image_id, N=args$num_samples)
}

main()
