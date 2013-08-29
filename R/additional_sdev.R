#' @export
additional_sdev <- function(X, W) {
    nc <- ncol(W)
    sdev <- numeric(nc)
    Q <- qr.Q(qr(W))
    Xp <- X
    for (cc in seq_len(nc)) {
        sdev[cc] <- sd(Xp%*%W[ , cc])
        Xp <- Xp - Xp%*%Q[ , cc]%*%t(Q[ , cc])   
    }
    return(sdev)
}