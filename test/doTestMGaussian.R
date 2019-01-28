
source("../emEstimationHMM.R")
source("../getAlphaHat.R")
source("../viterbi.R")
source("../computeHMMCrossValidatedLogLikelihood.R")
source("../../../math/weightedCov.R")
require(signal)
require(mvtnorm)
require(MASS)

getInitGaussianParams <- function(x, K) {
    D <- nrow(x)
    kMeansRes <- kmeans(x=t(x), centers=K)
    assignation <- kMeansRes$cluster
    mu <- matrix(NA, nrow=D, ncol=K)
    sigma <- array(NA, dim=c(D, D, K))
    for(k in 1:K) {
        indices <- which(assignation==k)
        mu[,k] <- rowMeans(x[,indices])
        sigma[,,k] <- cov(t(x[,indices]))
    }
    phi <- list(mu=mu, sigma=sigma)
    return(phi)
}
getUpdatedGaussianParams <- function(x, gamma, phi) {
    mu <- getMu(x=x, gamma=gamma)
    sigma <- getSigma(x=x, gamma=gamma, mu=mu)
    phi <- list(mu=mu, sigma=sigma)
    return(phi)
}
getMu <- function(x, gamma) {
    gammaTilde <- gamma%*%diag(1/colSums(gamma))
    mu <- x%*%gammaTilde
    return(mu)
}
getSigma <-  function(x, gamma, mu) {
    D <- nrow(x)
    K <- ncol(gamma)
    sigma <- array(NA, dim=c(D, D, K))
    gammaTilde <- gamma%*%diag(1/colSums(gamma))
    for(k in 1:K) {
        sigma[,,k] <- weightedCov(x=x, mu=mu[,k], weight=gammaTilde[,k])
    }
    return(sigma)
}
getGaussianProbabilities <- function(x, phi) {
    N <- ncol(x)
    K <- ncol(phi$mu)
    p <- matrix(NA, nrow=N, ncol=K)
    for(n in 1:N) {
        for(k in 1:K) {
            p[n,k] <- dmvnorm(x=x[,n], mean=phi$mu[,k], sigma=phi$sigma[,,k])
        }
    }
    return(p)
}

processAll <- function() {
    nSamplesTrain <- 150
    nSamplesTest <- 50
    m1 <- c(0,1)
    sd1 <- matrix(c(1,0.7,.7,1),2,2)
    m2 <- c(1,0)
    sd2 <- matrix(c(2,.1,.1,1),2,2)
    set.seed(2)

    y1Train <- mvrnorm(nSamplesTrain/2,m1,sd1)
    y2Train <- mvrnorm(nSamplesTrain/2,m2,sd2)
    # this creates data with a single change point
    yTrain <- rbind(y1Train,y2Train)
    
    y1Test <- mvrnorm(nSamplesTest/2,m1,sd1)
    y2Test <- mvrnorm(nSamplesTest/2,m2,sd2)
    # this creates data with a single change point
    yTest <- rbind(y1Test,y2Test)
    
    K <- 3
    estimationRes <- emEstimationHMM(x=t(yTrain),
                                      K=K, 
                                      getInitEmissionModelParams=
                                       getInitGaussianParams,
                                      getEmissionProbabilities=
                                       getGaussianProbabilities,
                                      getUpdatedEmissionModelParams=
                                       getUpdatedGaussianParams,
                                      convergenceTol=1e-4)
    vpath <- viterbi(Pi=estimationRes$Pi, p=estimationRes$p, A=estimationRes$A)
    gamma <- estimationRes$gamma
    PiTest <- gamma[nrow(gamma),]
    crossValidatedLL <- 
     computeHMMCrossValidatedLogLikelihood(x=t(yTest), 
                                            Pi=PiTest, 
                                            A=estimationRes$A,
                                            phi=estimationRes$phi,
                                            getEmissionProbabilities=
                                             getGaussianProbabilities)
    png("figures/testMGaussian.png")
    plot(vpath, xlab="Sample", ylab="State", yaxt="n", 
                main=sprintf("LL=%f", crossValidatedLL))
    abline(v=nrow(y1Train), col="red")
    axis(side=2, at=c(1,2))
    dev.off()

    browser()
}

processAll()

rm(processAll)
