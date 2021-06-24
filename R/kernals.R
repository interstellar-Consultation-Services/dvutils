#' @export
rbf <- function (X, sigma = 1, Y = NULL) {
  X <- as.matrix(X)
  if (!is.matrix(X)) {
    print("X must be a matrix containing samples in its rows")
    return()
  }
  if (length(sigma) != 1 || sigma <= 0) {
    print("sigma must be a number > 0 specifying the rbf-kernel width")
    return()
  }
  if (!is.null(Y)) {
    if (!is.matrix(Y)) {
      print("Y must be a matrix containing samples in its rows or NULL if it should not be used")
      return()
    }
    if (ncol(X) != ncol(Y)) {
      print("The samples in the rows of X and Y must be of same dimension")
      return()
    }
  }
  n <- nrow(X)
  if (is.null(Y)) {
    XtX <- tcrossprod(X)
    XX <- matrix(1, n) %*% diag(XtX)
    D <- XX - 2 * XtX + t(XX)
  }
  else {
    m <- nrow(Y)
    XX <- matrix(apply(X^2, 1, sum), n, m)
    YY <- matrix(apply(Y^2, 1, sum), n, m, byrow = TRUE)
    XY <- tcrossprod(X, Y)
    D <- XX - 2 * XY + YY
  }
  exp(-D/(2 * sigma))
}
