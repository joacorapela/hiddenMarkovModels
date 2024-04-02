
weightedCov <- function(x, mu, weight) {
    N <- ncol(x)
    D <- nrow(x)
    wCov <- matrix(0, nrow=D, ncol=D)
    for(n in 1:N) {
        wCov <- wCov + weight[n]*outer(X=x[,n]-mu, Y=x[,n]-mu)
    }
    return(wCov)
}
