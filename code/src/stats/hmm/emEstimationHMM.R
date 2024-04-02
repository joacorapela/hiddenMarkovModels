initProbabilityMatrix <- function(K) {
    probabilityMatrix <- matrix(rep(NA, times=K^2), nrow=K)
    for(i in 1:K) {
        probabilityMatrix[i,] <- initProbabilityVector(K=K)
    }
    return(probabilityMatrix)
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
                               getInitialStateProbability=initProbabilityVector,
                               getInitialTransitionMatrix=initProbabilityMatrix,
                               getInitEmissionModelParams,
                               getEmissionProbabilities,
                               getUpdatedEmissionModelParams,
                               convergenceTol=1e-4,
                               iterDisplayFn=NULL,
                               saveEachNIter=.Machine$integer.max,
                               saveFilename=NULL) {
    # Following notation from Bishop, 2006, Ch13
    # forward-backward algorithm or Baum-Welch algorithm
    #
    # answers (N number of observations, K number of latents):
    #
    # logLikelihoods (list of length N): history of log likelihoods achieved during EM
    # optimization. Since EM guarantees to increase the log likelihood at every
    # step, this list should contain increasing values.
    #
    # A (KxK matrix): state transition probabilities.
    #
    # gamma (NxK matrix): marginal posterior distribution of latents given all the
    # observations (i.e., p(z_n | x_1, \ldots, x_N).
    #
    # p (NxK matrix): emission probabilities (i.e., p(x_n | z_n=k))
    #
    # Pi (K vector): initial state probabilities
    #
    # phi: emission probability parameters. For a Gaussian HMM
    # p(X|Z=k)=N(X|mu[:,k], sigma[:,:,k]) and phi is a list
    # with elements mu and sigma. mu is a matrix of size PxK and sigma a tensor
    # of size PxPxK

    N <- ncol(x)
    Pi <- getInitialStateProbability(K=K)
    A <- getInitialTransitionMatrix(K=K)
    phi <- getInitEmissionModelParams(x=x, K=K)
    logLikelihood <- -Inf
    logLikelihoods <- c(logLikelihood)

    alpha <- matrix(NA, nrow=N, ncol=K)
    beta <- matrix(NA, nrow=N, ncol=K)
    iter <- 1
    while(TRUE) {
        if(!is.null(iterDisplayFn)) {
            iterDisplayFn(phi=phi)
        }
        # start E step
        prevLogLikelihood <- logLikelihood
        res <- computeEStep(x=x, Pi=Pi, A=A, phi=phi, getEmissionProbabilities=getEmissionProbabilities)
        gamma <- res$gamma
        p <- res$p
        c <- res$c
        alphaHat <- res$alphaHat
        betaHat <- res$betaHat

        prevLogLikelihood <- logLikelihood
        logLikelihood <- sum(log(c)) # (13.63)
        if(is.nan(logLikelihood)) {
            print("LogLikelihood is NaN; returning NULL")
            browser()
            return(NULL)
        }
        logLikelihoods <- c(logLikelihoods, logLikelihood) 
        if(logLikelihood<prevLogLikelihood) {
            warning(sprintf("Likelkihood decreased from %f to %f", prevLogLikelihood, logLikelihood))
        }
        logLikelihoodDiff <- logLikelihood-prevLogLikelihood
        print(sprintf("LogLikelihood difference %f (%f). LogLikelihood %f", 
                     logLikelihoodDiff, convergenceTol, logLikelihood))
        if(logLikelihoodDiff<convergenceTol) {
            break
        }
        # end E step

        # start M step
        xi <- getXi(alphaHat=alphaHat, betaHat=betaHat, c=c, p=p, A=A)
        Pi <- getPi(gamma=gamma)
        A <- getA(xi=xi)
        phi <- getUpdatedEmissionModelParams(x=x, gamma=gamma, phi=phi)
        # end M step
        if(iter%%saveEachNIter==0) {
            answer <- list(logLikelihoods=logLikelihoods, A=A, gamma=gamma, p=p,
                                                          Pi=Pi, phi=phi)
            print(sprintf("Saving %s", saveFilename))
            save(answer, file=saveFilename)
        }
        iter <- iter + 1
    }
    answer <- list(logLikelihoods=logLikelihoods, A=A, gamma=gamma, p=p,
                                                  Pi=Pi, phi=phi)
    print(sprintf("Exiting with logLikelihood difference %f and logLikelihood %f", logLikelihoodDiff, logLikelihood))
    return(answer)
}

