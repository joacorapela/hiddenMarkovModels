
source("../src/stats/hmm/emEstimationHMM.R")
source("../src/stats/hmm/viterbi.R")
require(signal)
require(mvtnorm)

getInitPoissonParams <- function(x, K) {
    kMeansRes <- kmeans(x=x[1,], centers=K)
    assignation <- kMeansRes$cluster
    lambda <- array(NA, dim=K)
    for(k in 1:K) {
        indices <- which(assignation==k)
        lambda[k] <- mean(x[1,indices])
    }
    phi <- list(lambda=lambda)
    return(phi)
}
getUpdatedPoissonParams <- function(x, gamma, phi) {
    gammaTilde <- gamma%*%diag(1/colSums(gamma))
    lambda <- x%*%gammaTilde
    phi <- list(lambda=lambda)
    return(phi)
}
getPoissonProbabilities <- function(x, phi) {
    N <- ncol(x)
    K <- length(phi$lambda)
    p <- matrix(NA, nrow=N, ncol=K)
    for(n in 1:N) {
        for(k in 1:K) {
            p[n,k] <- dpois(x=x[,n], lambda=phi$lambda[k])
        }
    }
    return(p)
}

processAll <- function() {
    # generate data from two different Poisson distributions
    nSamples <- 300
    lambda1 <- 1
    lambda2 <- 3.5
    K <- 2

    y1 <- rpois(n=nSamples/2, lambda=lambda1)
    y2 <- rpois(n=nSamples/2, lambda=lambda2)
    # this creates data with a single change point
    y <- matrix(c(y1,y2), nrow=1)
    
    estimationRes <- emEstimationHMM(x=y,
                                      K=K, 
                                      getInitEmissionModelParams=
                                       getInitPoissonParams,
                                      getEmissionProbabilities=
                                       getPoissonProbabilities,
                                      getUpdatedEmissionModelParams=
                                       getUpdatedPoissonParams,
                                      convergenceTol=1e-4)
    vpath <- viterbi(Pi=estimationRes$Pi, p=estimationRes$p, A=estimationRes$A)

    png("figures/testPoisson.png")
    plot(vpath, xlab="Sample", ylab="State", yaxt="n")
    abline(v=length(y1), col="red")
    axis(side=2, at=c(1,2))
    dev.off()

    browser()
}

processAll()

rm(processAll)
