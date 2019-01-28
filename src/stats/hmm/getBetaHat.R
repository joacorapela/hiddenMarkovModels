
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
