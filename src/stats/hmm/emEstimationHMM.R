initProbabilityVector <- function(K) {
    probabilityVector <- rep(NA, times=K)
    probabilityVectorBuilt <- FALSE
    while(TRUE) {
        probabilityVector[1:(K-1)] <- runif(n=K-1, min=0, max=2/K)
        sum <- sum(probabilityVector[1:(K-1)])
        if(sum<1) {
            probabilityVector[K] <- 1-sum
            return(probabilityVector)
        }
    }
}
initProbabilityMatrix <- function(K) {
    probabilityMatrix <- matrix(rep(NA, times=K^2), nrow=K)
    for(i in 1:K) {
        probabilityMatrix[i,] <- initProbabilityVector(K=K)
    }
    return(probabilityMatrix)
}
getAlphaHat <- function(Pi, A, p) { # (13.59)
    N <- nrow(p)
    K <- ncol(p)
    alphaHat <- matrix(NA, nrow=N, ncol=K)
    c <- rep(NA, times=N)
    c[1] <- t(p[1,])%*%Pi
    alphaHat[1,] <- p[1,]/c[1] 
    for(n in 2:N) {
        c[n] <- t(alphaHat[n-1,])%*%A%*%p[n,]
        alphaHat[n,] <- 1/c[n]*(t(alphaHat[n-1,])%*%A)%*%diag(p[n,])
    }
    answer <- list(alphaHat=alphaHat, c=c)
    return(answer)
}
getBetaHat <- function(A, p, c) {# (13.62)
    N <- nrow(p)
    K <- nrow(A)
    betaHat <- matrix(NA, nrow=N, ncol=K)
    betaHat[N,] <- rep(1, K)
    for(n in (N-1):1) {
        betaHat[n,] <- A%*%(betaHat[n+1,]*p[n+1,])/c[n+1] 
    }
    return(betaHat)
}
getXi <- function(alphaHat, betaHat, c, p, A) { # (13.65)
    K <- nrow(A)
    N <- nrow(p)
    xi <- array(NA, dim=c(K, K, N-1))

    for(n in 1:(N-1)) {
        xi[,,n] <- 1/c[n+1]*diag(alphaHat[n,])%*%A%*%diag(p[n+1,]*betaHat[n+1,])
    }
    return(xi)
}
getPi <- function(gamma) { # (13.18)
    Pi <- gamma[1,]/sum(gamma[1,])
    return(Pi)
}
getA <- function(xi) { # (13.19)
    K <- dim(xi)[1]
    N <- dim(xi)[3]
    A <- matrix(NA, nrow=K, ncol=K)
    for(j in 1:K) {
        denominator <- sum(xi[j,,])
        for(k in 1:K) {
            numerator <- sum(xi[j,k,])
            A[j,k] <- numerator/denominator
        }
    }
    return(A)
}
    
emEstimationHMM <- function(x, K, 
                               getInitEmissionModelParams,
                               getEmissionProbabilities,
                               getUpdatedEmissionModelParams,
                               convergenceTol=1e-4,
                               iterDisplayFn=NULL) {
    # Following notation from Bishop, 2006, Ch13
    # forward-backward algorithm or Baum-Welch algorithm
    N <- ncol(x)
    Pi <- initProbabilityVector(K=K)
    A <- initProbabilityMatrix(K=K)
    phi <- getInitEmissionModelParams(x=x, K=K)
    logLikelihood <- -Inf
    logLikelihoods <- c(logLikelihood)

    alpha <- matrix(NA, nrow=N, ncol=K)
    beta <- matrix(NA, nrow=N, ncol=K)
    while(TRUE) {
        if(!is.null(iterDisplayFn)) {
            iterDisplayFn(phi=phi)
        }
        p <- getEmissionProbabilities(x=x, phi=phi)
        res <- getAlphaHat(Pi=Pi, A=A, p=p)
        alphaHat <- res$alphaHat
        c <- res$c
        prevLogLikelihood <- logLikelihood
        logLikelihood <- sum(log(c)) # (13.63)
        if(is.nan(logLikelihood)) {
            return(NULL)
        }
        logLikelihoods <- c(logLikelihoods, logLikelihood) 
        if(logLikelihood<prevLogLikelihood) {
            warning(sprintf("Likelkihood decreased from %f to %f", prevLogLikelihood, logLikelihood))
        }
        logLikelihoodDiff <- logLikelihood-prevLogLikelihood
        show(sprintf("LogLikelihood difference %f (%f). LogLikelihood %f", 
                     logLikelihoodDiff, convergenceTol, logLikelihood))
        if(logLikelihoodDiff<convergenceTol) {
            break
        }
        betaHat <- getBetaHat(A=A, p=p, c=c)
        gamma <- alphaHat*betaHat # (13.64)
        xi <- getXi(alphaHat=alphaHat, betaHat=betaHat, c=c, p=p, A=A)
        Pi <- getPi(gamma=gamma)
        A <- getA(xi=xi)
        phi <- getUpdatedEmissionModelParams(x=x, gamma=gamma, phi=phi)
    }
    answer <- list(logLikelihoods=logLikelihoods, A=A, gamma=gamma, p=p,
                                                  Pi=Pi, phi=phi)
    show(sprintf("Exiting with logLikelihood difference %f and logLikelihood %f", logLikelihoodDiff, logLikelihood))
    return(answer)
}

