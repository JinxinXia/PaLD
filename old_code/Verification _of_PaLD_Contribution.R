#' Get Contribution Matrix
#'
#' @param distance_matrix A distance matrix. Can be obtained from the
#'    `get_distance_matrix()` function.
#' @param beta Numeric. A localizing parameter. By default, `beta = 1`, which reflects
#'   the natural choice of commonality foci. This should be a positive number.
#' @param digits Numeric. The number of allowed digits for your `distance_matrix`.
#'   Default is 15. If you do not want rounding, you can set this to `Inf`.
#'
#' @return A matrix of contributions.
#' @export
#'
#' @examples
#' d <- data.frame(
#'   X1 = sample(20, 10),
#'   X2 = sample(20, 10)
#'   )
#' distance_mat <- get_distance_matrix(d)
#' contribution_mat <- get_contribution_matrix(distance_mat)
get_contribution_matrix <- function(distance_matrix, beta = 1, digits = 15){
  if (beta < 0) {
    stop("The `beta` parameter must be positive.", .call = FALSE)
  }
  distance_matrix <- round(distance_matrix, digits = digits)
  n <- dim(distance_matrix)[1]
  a <- matrix(0, n, n)
  for(x in 1:(n - 1)){
    for(y in (x + 1):n){
      dx <- distance_matrix[x, ]
      dy <- distance_matrix[y, ]
      uxy <- which((dx <= beta * distance_matrix[x, y]) | (dy <= beta * distance_matrix[y, x]))
      wx <- 1 * (dx[uxy] < dy[uxy]) + 0.5 * ((dx[uxy] == dy[uxy]))
      
      wy <- 1 * (dy[uxy] < dx[uxy]) + 0.5 * ((dx[uxy] == dy[uxy]))
      a[x, uxy] <- a[x, uxy] + (1 / (length(uxy)) * wx)
      a[y, uxy] <- a[y, uxy] + (1 / (length(uxy)) * wy)
    }
  }
  diag(a) <- diag(a)
  if (is.null(rownames(distance_matrix))) {
    rownames(distance_matrix) <- 1:n
  }
  rownames(a) <- rownames(distance_matrix)
  colnames(a) <- 1:n
  return(a / (n - 1))
}


D <- matrix(1:9,nrow = 3,ncol=3)

D[1,1] <- 0
D[1,2] <- 0.5770
D[1,3] <- 0.4223
D[2,1] <- 0.5770
D[2,2] <- 0
D[2,3] <- 0.4130
D[3,1] <- 0.4223
D[3,2] <- 0.4130
D[3,3] <- 0


get_contribution_matrix(D,beta=1)


