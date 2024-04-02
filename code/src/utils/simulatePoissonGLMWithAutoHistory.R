simulatePoissonGLMWithAutoHistory <- function(f, N, h) {
    g <- function(h, gamma) {
        return(f(sum(h*gamma)))
    }

    L <- length(h)-1
    xZeroPadded <- rep(x=0, times=L)
    lambda <- rep(x=NA, times=N)
    for(n in (1:N)) {
        rho <- xZeroPadded[n+L-(1:L)]
        gValue <- g(h=h, gamma=c(1, rho))
        lambda[n] <- gValue
        xZeroPadded <- c(xZeroPadded, rpois(n=1, lambda=lambda[n]))
    }
    x <- xZeroPadded[(L+1):length(xZeroPadded)]
    answer <- list(x=x, lambda=lambda)
    return(answer)
}
