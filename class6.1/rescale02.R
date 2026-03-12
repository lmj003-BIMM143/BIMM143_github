rescale02 <- function(x, na.rm=TRUE) {
  rng <- range(x, na.rm=na.rm)
  (x-rng[1]) / (rng[2] - rng[1])
}


rescale <- function(x, na.rm=TRUE, plot=FALSE) {
  rng <- range(x, na.rm=na.rm)
  answer <- (x-rng[1]) / (rng[2] - rng[1])
  print ("is it me you are looking for?")
  if(plot) {
    plot(answer, type="b", lwd=4)
  }
  print("I can see it in...")
  return(answer)
}

rescale(1:10)

rescale(1:10, plot=TRUE)