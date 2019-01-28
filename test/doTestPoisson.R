
source("../emEstimationHMM.R")
source("../viterbi.R")
source("../computeHMMCrossValidatedLogLikelihood.R")
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
    nSamplesTrain <- 300
    nSamplesTest <- 100
    lambda1 <- 1
    lambda2 <- 3.5
    K <- 2

    y1Train <- rpois(n=nSamplesTrain/2, lambda=lambda1)
    y2Train <- rpois(n=nSamplesTrain/2, lambda=lambda2)
    # this creates data with a single change point
    yTrain <- matrix(c(y1Train,y2Train), nrow=1)
    
    y1Test <- rpois(n=nSamplesTest/2, lambda=lambda1)
    y2Test <- rpois(n=nSamplesTest/2, lambda=lambda2)
    # this creates data with a single change point
    yTest <- matrix(c(y1Test,y2Test), nrow=1)
    
    estimationRes <- emEstimationHMM(x=yTrain,
                                      K=K, 
                                      getInitEmissionModelParams=
                                       getInitPoissonParams,
                                      getEmissionProbabilities=
                                       getPoissonProbabilities,
                                      getUpdatedEmissionModelParams=
                                       getUpdatedPoissonParams,
                                      convergenceTol=1e-4)
    vpath <- viterbi(Pi=estimationRes$Pi, p=estimationRes$p, A=estimationRes$A)

    gamma <- estimationRes$gamma
    PiTest <- gamma[nrow(gamma),]
    crossValidatedLL <- 
     computeHMMCrossValidatedLogLikelihood(x=yTest, 
                                            Pi=PiTest, 
                                            A=estimationRes$A,
                                            phi=estimationRes$phi,
                                            getEmissionProbabilities=
                                             getPoissonProbabilities)
    png("figures/testPoisson.png")
    plot(vpath, xlab="Sample", ylab="State", yaxt="n",
                main=sprintf("LL=%f", crossValidatedLL))
    abline(v=length(y1Train), col="red")
    axis(side=2, at=c(1,2))
    dev.off()

    browser()
}

processAll()

rm(processAll)
