
source("../src/stats/hmm/initProbabilityVector.R")
source("../src/stats/hmm/computeEStep.R")
source("../src/stats/hmm/emEstimationHMM.R")
source("../src/stats/hmm/getAlphaHat.R")
source("../src/stats/hmm/getBetaHat.R")
source("../src/stats/hmm/computeHMMCrossValidatedLogLikelihood.R")
source("../src/stats/hmm/getARProbabilities.R")
source("../src/stats/hmm/buildXHist.R")

getInitARParamsFunc <- function(L, nSamplesKMean=1e4, sigmaSqMean=.1, sigmaSqSD=.01, aSD=.1) {
    answer <- function(x, K) {
        x <- matrix(sample(x, size=nSamplesKMean), nrow=nrow(x))
        kMeansRes <- kmeans(x[1,], centers=K)
        a <- matrix(0, nrow=L+1, ncol=K)
        sigmaSq <- rep(0, times=K)
        for(k in 1:K) {
            a[1,k] <- mean(x[kMeansRes$cluster==k])
            sigmaSq[k] <- var(x[kMeansRes$cluster==k])
        }
        phi <- list(a=a, sigmaSq=sigmaSq)
        return(phi)
    }
    return(answer)
}
iterDisplayFn <- function(phi) {
    show("a:")
    show(phi$a)
    show("sigmaSq:")
    show(phi$sigmaSq)
}
processAll <- function() {
    sigmaSqMean <- .1
    sigmaSqSD <- .01
    aSD <- .1
    percentageTrain <- .5
    Ls <- 10
    Ks <- 2
    convergenceTol <- 1e-6
    # prefix <- "results/simulationResAR_delta1N10000a10.100_0.200_0.200_0.150_0.100_0.050_0.000_-0.050_-0.025_-0.010MuOne0.00SigmaSqOne0.2a2-0.30_-0.20_-0.10_0.00_0.05_0.10_0.15_0.20_0.15_0.05MuTwo0.00SigmaSqTwo0.1_train"
    prefix <- "simulationResAR_delta1N10000a10.100_0.200_0.200_0.150_0.100_0.050_0.000_-0.050_-0.025_-0.010MuOne1.00SigmaSqOne0.2a2-0.30_-0.20_-0.10_0.00_0.05_0.10_0.15_0.20_0.15_0.05MuTwo2.00SigmaSqTwo0.1_train"
    dataFilenamePattern <- "data/%s.RData"
    resultsFilenamePattern <- "results/%s_%dDStateL%d_hmmEstimatedParams.RData"

    getARProbabilitiesWrapper <- function(x, phi) {
        return(getARProbabilities(x=x, phi=phi, xHist=xHist))
    }
    getUpdatedARParams <- function(x, gamma, phi, ...) {
        # x \in \Re^{N \times 1}
        # gamma \in Re^{N\times K}
        # phi=list(sigma, a): sigma\in Re^K, a\in Re^{L\times K}

        df <- data.frame(y=x[1,], X=I(xHist[1,,]))
        a <- matrix(NA, nrow=L+1, ncol=K)
        sigmaSq <- array(NA, dim=K)
        for(k in 1:K) {
            # first compute
            lmRes <- lm(y~X, data=df, weights=gamma[,k])
            coefs <- coefficients(lmRes)
            a[,k] <- coefs

            # next compute sigma^2
            auxVector <- x[1,]-cbind(1, xHist[1,,])%*%coefs
            numerator <- sum(gamma[,k]*auxVector^2)
            denominator <- sum(gamma[,k])
            sigmaSq[k] <- numerator/denominator
        }
        phi <- list(a=a, sigmaSq=sigmaSq)
        return(phi)
    }
    dataFilename <- sprintf(dataFilenamePattern, prefix)
    loadRes <- get(load(dataFilename))
    x <- matrix(loadRes$x, nrow=1)

    for(L in Ls) {
        xHist <- buildXHist(x=x, L=L)
        for(K in Ks) {
            getInitARParams <- getInitARParamsFunc(L=L, sigmaSqMean=sigmaSqMean,
                                                        sigmaSqSD=sigmaSqSD,
                                                        aSD=aSD)
            estimationRes <- emEstimationHMM(x=x, 
                                              K=K, 
                                              getInitEmissionModelParams=
                                               getInitARParams,
                                              getEmissionProbabilities=
                                               getARProbabilitiesWrapper,
                                              getUpdatedEmissionModelParams=
                                               getUpdatedARParams,
                                              convergenceTol=convergenceTol,
                                              iterDisplayFn=iterDisplayFn)
            resultsFilename <- sprintf(resultsFilenamePattern, prefix, K, L)
            save(estimationRes, file=resultsFilename)
        }
    }
    for(k in 1:K) {
        stateMean <- estimationRes$phi$a[1,k]/(1-sum(estimationRes$phi$a[-1,k]))
        show(sprintf("mu[%d]=%f", k, stateMean))
    }

    browser()
}

processAll()

rm(processAll)
