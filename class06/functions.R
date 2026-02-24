add <- function(x, y=1) {
  x + y
}

add(x=4)

add(x=4, y=10)

add(c(1, 2, 5, 10), y=10)
add(c(1, 2, 5, 10), y=10, 20)
add(c(1, 2, 5, 10), y="string")

# Our 2nd function rescale()
rescale <- function(x) {
  rng <- range(x)
  (x-rng[1]) / (rng[2] - rng[1])
}

# Test on a simple example
rescale(1:10)

rescale(c(1, 4, 10, 20))

rescale(c(1, NA, 4, 10, 20))

# Fix the NA issue
rescale01 <- function(x) {
  rng <- range(x, na.rm=TRUE)
  (x-rng[1]) / (rng[2] - rng[1])
}

rescale01(c(1, NA, 4, 10, 20))

rescale01(c(1, 4, "word"))





