func1 <- function() {
  a <- rep(NA, 10)
  func0 <- function(i, x) {
    a[i] <- x
    print(a)
    a
  }

  for (i in 1:10) {
    a <- func0(i, rpois(1, lambda=10))
  }
}
