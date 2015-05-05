# script for outlier test assuming negative binomial distribution
# after estimating the parameter of negative binomial distribution, 
# the maximum number is checked whether the estimated parameter can resonably explain the value,
# and if the p-value is below the threshold, the junction is filtered out.

thresPValue <- 5;

# function for calculating the part of likelihood of negative binomial distribution (mainly for estimating the "r")
partNBLiklihood <- function(param, data) {

    N <- length(data);
    return( sum(lgamma(data + rep(param - 1, N))) - N * lgamma(param) + N * param * log(param / (param + mean(data))) );
}

# function for getting the p-value of the maximum number
getPout <- function(x) {

    res <- optim(1, partNBLiklihood, gr = NULL, method = "Brent", lower = 0.001, upper = 100, data = as.numeric(x));
    r_est <- res$par;
    p_est <- mean(x) / (mean(x) + r_est);
    return(-log10( 1 - pnbinom(max(x) - 1, r_est, 1 - p_est)))

}

inputFile <- commandArgs(trailingOnly=TRUE)[1]
outputFile <- commandArgs(trailingOnly=TRUE)[2]

bdata <- read.table(inputFile, sep="\t", header=FALSE, stringsAsFactors = FALSE);
juncInfo <- bdata[, 1:3];
countInfo <- bdata[, 4:ncol(bdata)];

Pout <- apply(countInfo, 1, getPout)
Pout[Pout == Inf] <- 60

write.table(cbind(juncInfo[Pout >= thresPValue, ], Pout[Pout >= thresPValue], countInfo[Pout >= thresPValue, ]), outputFile, sep="\t", row.names= FALSE, col.names = FALSE, quote = FALSE)

