
source("../src/stats/hmm/emEstimationHMM.R")
source("../src/stats/hmm/fEscola.R")
source("../src/stats/hmm/fEscolaDer.R")
require(signal)

getInitPoissonGLMParams <- function(x, K, nb=2, noiseSD=.5) {
#     kMeansRes <- kmeans(x=x[1,], centers=K)
#     assignation <- kMeansRes$cluster
#     lambda <- array(NA, dim=K)
#     for(k in 1:K) {
#         indices <- which(assignation==k)
#         lambda[k] <- mean(x[1,indices])
#     }
    # hTrue <- cbind(c(.4, -.3), c(0, 0))
    # h <- hTrue + rnorm(length(hTrue), mean=0, sd=.1)
    h <- matrix(rnorm((nb+1)*K, mean=0, sd=noiseSD), ncol=K)
    phi <- list(h=h)
    return(phi)
}
getLogProb <- function(xn, xnHist, h, f) { 
    # get log(p(xn|\gamma_n, \theta))
    # h \in \Re^L
    g <- function(h, gamma) {
        return(f(sum(h*gamma)))
    }

    gValue <- g(h=h, gamma=c(1, xnHist))
    lnProb <- xn*log(gValue)-gValue-lfactorial(xn)
    return(lnProb)
}
getLogProbXn <- function(xn, xnHist, h, f) {# h \in \Re^L
    sum <- 0
    for(i in 1:length(xn)) {
        sum <- sum + getLogProb(xn=xn[i], xnHist=xnHist[i,], h=h, f=f)
    }
    return(sum)
}
getGradLogProb <- function(xn, xnHist, h, f, fDer) { 
    # get \nabla log(p(xn|\gamma_n, \theta))
    # h \in \Re^L
    g <- function(h, gamma) {
        return(f(sum(h*gamma)))
    }
    gGrad <- function(h, gamma) {
        return(fDer(sum(h*gamma))*gamma)
    } 

    gValue <- g(h=h, gamma=c(1, xnHist))
    gGradValue <- gGrad(h=h, gamma=c(1, xnHist))
    gradLogProb <- (xn/gValue-1)*gGradValue
    return(gradLogProb)
}
getGradLogProbXn <- function(xn, xnHist, h, f, fDer) { # h \in \Re^L
    sum <- rep(0, times=length(h))
    for(i in 1:length(xn)) {
        sum <- sum + getGradLogProb(xn=xn[i], xnHist=xnHist[i,], h=h, 
                                              f=f, fDer=fDer)
    }
    return(sum)
}
iterDisplayFn <- function(phi) {
    show("h:")
    show(phi$h)
}
buildXHist <- function(x, L) {
    D <- nrow(x)
    N <- ncol(x)
    xHist <- array(NA, dim=c(D, N, L))
    zeroPad <- matrix(0, nrow=D, ncol=L)
    xZeroPadded <- cbind(zeroPad, x)
    for(n in 1:N) {
        xHist[, n,] <- xZeroPadded[,n+L-(1:L)]
    }
    return(xHist)
}
reduceDimensionalityXHist <- function(xHist, Phi) {
    D <- dim(xHist)[1]
    N <- dim(xHist)[2]
    nb <- ncol(Phi)
    xHistReduced <- array(NA, dim=c(D, N, nb))
    for(n in 1:N) {
        xHistReduced[, n,] <- xHist[,n,]%*%Phi
    }
    return(xHistReduced)
}
processAll <- function() {
    getUpdatedPoissonGLMParams <- function(x, gamma, phi,
                                              f=fEscola, fDer=fEscolaDer,
                                              method="BFGS", reltol=1e-6, ...) {
        # h \in \Re^{L \times K}
        h <- phi$h
        L <- dim(xHist)[3]
        N <- ncol(x)
        D <- nrow(x)
        K <- ncol(phi$h)
        for(k in 1:K) {
            qProp <- function(h) { # h \in \Re^L
                sum <- 0
                for(n in 1:N) {
                    lnProbXn <- getLogProbXn(xn=x[,n], xnHist=xHist[,n,], h=h, f=f)
                    sum <- sum + gamma[n,k]*lnProbXn
                }
                return(sum)
            }
            qPropGrad <- function(h) { # h \in \Re^L
                sum <- 0
                for(n in 1:N) {
                    gradLogProbXn <- 
                     getGradLogProbXn(xn=x[,n], xnHist=xHist[,n,], h=h, f=f, 
                                                fDer=fDer)
                    sum <- sum + gamma[n,k]*gradLogProbXn
                }
                return(sum)
            }
            res <- optim(par=h[,k], fn=qProp, gr=qPropGrad, method=method,
                                    control=list(reltol=reltol, fnscale=-1))
            h[,k] <- res$par
        }
        phi <- list(h=h)
        return(phi)
    }
    getPoissonGLMProbabilities <- function(x, phi) {
        h <- phi$h
        L <- dim(xHist)[3]
        N <- ncol(x)
        D <- nrow(x)
        K <- ncol(phi$h)
        p <- matrix(0, nrow=N, ncol=K)
        for(n in 1:N) {
            for(k in 1:K) {
                p[n,k] <- exp(getLogProbXn(xn=x[,n], xnHist=xHist[,n,], 
                                                     h=h[,k], 
                                                     f=fEscola))
            }
        }
        return(p)
    }
    K <- 2
    convergenceTol <- 1e-5
    basisFilename <- "data/raisedCosBasisNBasis2Dt0.0100EndPoints_0.01_0.04B0.00ZFlag0.RData"
    dataFilename <- "data/simulationResPoissonGLMWithAutohistoryD10Delta0.01N3000Mu10.40C10.20_-0.10Mu20.20C20.00_0.00.RData"
    resultsFilenamePattern <- "results/simulationPoissonGLMWithAutohistoryD10Delta0.01N3000Mu10.40C10.20_-0.10Mu20.20C20.00_0.00_hmmEstimatedParams_%dDState.RData"

    res <- get(load(basisFilename))
    Phi <- res$nonOrthogonalBasis
    L <- nrow(Phi)
    resultsFilename <- sprintf(resultsFilenamePattern, K)

    loadRes <- get(load(dataFilename))
    x <- loadRes$x
    x <- matrix(x[1:3,], ncol=ncol(x))
    xHist <- buildXHist(x=x, L=L)
    xHist <- reduceDimensionalityXHist(xHist=xHist, Phi=Phi)
    estimationRes <- emEstimationHMM(x=x, K=K, 
                                          getInitEmissionModelParams=
                                           getInitPoissonGLMParams,
                                          getEmissionProbabilities=
                                           getPoissonGLMProbabilities,
                                          getUpdatedEmissionModelParams=
                                           getUpdatedPoissonGLMParams,
                                          convergenceTol=convergenceTol,
                                          iterDisplayFn=iterDisplayFn)
    save(estimationRes, file=resultsFilename)

    browser()
}

processAll()

rm(processAll)
