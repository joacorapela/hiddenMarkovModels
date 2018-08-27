
fEscola <- function(u) {
    if(u[1]<=0) {
        return(exp(u))
    }
    return(1+u+.5*u^2)
}
