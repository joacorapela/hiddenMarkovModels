
getAlphaHat <- function(Pi, A, p, tol=1e-6) { # (13.59)
    N <- nrow(p)
    K <- ncol(p)
    alphaHat <- matrix(NA, nrow=N, ncol=K)
    c <- rep(NA, times=N)
    c[1] <- t(p[1,])%*%Pi
    alphaHat[1,] <- p[1,]/c[1] 
    for(n in 2:N) {
        c[n] <- t(alphaHat[n-1,])%*%A%*%p[n,]
        alphaHat[n,] <- 1/c[n]*(t(alphaHat[n-1,])%*%A)*p[n,]
    }
    answer <- list(alphaHat=alphaHat, c=c)
    return(answer)
}
