
doubleGammaPoisson <- function(params, data) {
  
    t <- params[1];
    alpha <- params[2];
    beta <- params[3]

    x1 <- data[[1]]
    x2 <- data[[2]]
    N <- length(x1)

    ML <- sum(x1) * log(t)
    ML <- ML - sum(x1 + x2 + rep(alpha, N)) * log(t + beta + 1)
    ML <- ML + N * alpha * log(beta)
    ML <- ML + sum(vapply(x1 + x2 + rep(alpha, N), lgamma, numeric(1)))
    ML <- ML - N * lgamma(alpha) 
 
    return(-ML);  
  
}


# alpha and beta are fixed in this function
doubleGammaPoisson2 <- function(params, data) {

    t <- params;
    alpha <- data[[3]];
    beta <- data[[4]];

    x1 <- data[[1]]
    x2 <- data[[2]]
    N <- length(x1)

    ML <- sum(x1) * log(t)
    ML <- ML - sum(x1 + x2 + rep(alpha, N)) * log(t + beta + 1)
    ML <- ML + N * alpha * log(beta)
    ML <- ML + sum(vapply(x1 + x2 + rep(alpha, N), lgamma, numeric(1)))
    ML <- ML - N * lgamma(alpha)

    return(-ML);

}


bdata <- read.table(commandArgs()[5], sep="\t",header=F);
targetDepth_T <- bdata[, 1]
targetDepth_N <- bdata[, 2]
flankingDepth_T <- bdata[, 3]
flankingDepth_N <- bdata[, 4]
infoData <- bdata[, 5:ncol(bdata)]
PVs <- matrix(0, nrow(bdata), 1)
tumorRateCheck <- matrix(0, nrow(bdata), 1)

for (i in 1:nrow(bdata)) {

    y1 <- as.integer(unlist(strsplit(as.character(targetDepth_T[i]), ",")));
    y2 <- as.integer(unlist(strsplit(as.character(targetDepth_N[i]), ",")));
    removeInd <- (y2 < 1000)
    y1 <- y1[removeInd];
    y2 <- y2[removeInd];

    x1 <- as.integer(unlist(strsplit(as.character(flankingDepth_T[i]), ",")));
    x2 <- as.integer(unlist(strsplit(as.character(flankingDepth_N[i]), ",")));
    removeInd <- (x2 < 1000)
    x1 <- x1[removeInd];
    x2 <- x2[removeInd];


    if (length(y1) >= 1 & length(x1) >= 1) {

        res <- constrOptim(c(10, 10, 10), doubleGammaPoisson, grad = NULL, ui = diag(3), ci = c(0.01, 0.01, 0.01), data = list(x1, x2));
        param_est <- res$par

        L1 <- doubleGammaPoisson(param_est, list(y1, y2)) 
        res2 <- optim(c(param_est[1]), doubleGammaPoisson2, gr = NULL, method = "Brent", lower = 0.01, upper = 100, data = list(y1, y2, param_est[2], param_est[3]));

        PVs[i, 1] <- -log10(1 - pchisq(L1 - res2$value, df = 1));
        if (PVs[i, 1] == Inf) PVs[i, 1] <- 60;
        if (res2$par[1] > param_est[1]) {
            tumorRateCheck[i, 1] <- 1;
        }

    }
    if (PVs[i, 1] > 5 & tumorRateCheck[i, 1] == 0) {
        print(PVs[i, 1]);
        print(infoData[i,]);
    }
}

write.table(cbind(infoData[PVs >= 3 & tumorRateCheck == 0,,drop=FALSE], as.matrix(PVs)[PVs >= 3 & tumorRateCheck == 0,,drop=FALSE]), commandArgs()[6], sep="\t", row.names= FALSE, col.names = FALSE, quote = FALSE);

