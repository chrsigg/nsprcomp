#' @export
cardinality <- function(W) {
    return(colSums(abs(W) > 0))
}