
buildXHist <- function(x, L) {
    # xHist \in Re^{N\times L}
    D <- nrow(x)
    N <- ncol(x)
    xHist <- array(NA, dim=c(D, N, L))
    zeroPad <- matrix(0, nrow=D, ncol=L)
    xZeroPadded <- cbind(zeroPad, x)
    for(n in 1:N) {
        xHist[,n,] <- xZeroPadded[,n+L-(1:L)]
    }
    return(xHist)
}
