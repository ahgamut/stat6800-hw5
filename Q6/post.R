

get_data <- function() {
    files <- Sys.glob("./*.Rdata")
    data <- do.call(rbind, lapply(files, function(x) {
        load(x)
        # print(info)
        print(info$k)
        d1 <- data.frame(k=info$k, B=info$B,
                         true_stat=info$true_stat, 
                         avg_boot=mean(info$simul_stat), 
                         p_value=info$p_value)
        d1
    }))
    data <- data[order(data$k),]
}

main <- function() {
}
