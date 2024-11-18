
# post-processing
# run after calc.R and calc2.R

get_data <- function() {
    files <- Sys.glob("./*.Rdata")
    data <- do.call(rbind, lapply(files, function(x) {
        load(x)
        # print(res$name)
        iname <- gsub("\\.jpg", "", basename(res$name))
        colors <- c("R", "G", "B")
        d1 <- do.call(rbind, lapply(colors, function(col) {
            check1 <- res[[col]][[1]]
            check2 <- res[[col]][[2]]
            c1.pv <- check1$p_value
            c2.pv <- ifelse(is.null(check2), NA, check2$p_value)
            dc <- data.frame(name=iname, color=col, h0_ab=c1.pv, h0_beta=c2.pv)
            dc
        }))
        d1
    }))
    data
}

get_data()
