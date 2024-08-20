### Find all timepoints when time series crosses zero
findZeros <- function(time, var, lim = 0.0005, timeSplit = 10, return_n1n2 = FALSE)
{
  # Data is interpolated to pinpoint the right moment 
  # where irradiance change direction.
  timeStep <- time[2] - time[1]
  x_out <- seq(min(time), max(time), timeStep / timeSplit)
  int <- approx(x = time, y = var, xout = x_out)
  sgf <- calDev(x = int$x, y = int$y)
  
  zeros <- vector()
  n1n2 <- vector()
  for(i in 1:length(sgf))
  {
    if((i+1) > length(sgf)){next()}
    
    n1 = sgf[i]
    n2 = sgf[i+1]
    if(n1 * n2 < 0)
    {
      # A limit is defined to ignore the extremely small changes of direction
      # i.e. when light is mostly stable
      if(abs(n1 * n2) > lim)
      {
        zeros <- append(zeros, values = (i + 1))
        n1n2 <- append(n1n2, values = n1*n2)
        if(n1 < 0)
        {
          names(zeros)[length(zeros)] <- "low"
        } else {
          names(zeros)[length(zeros)] <- "up"
        }
      }
    }
  }
  
  if(return_n1n2 == TRUE)
  {
    dfZ <- data.frame("zeros" = zeros, "n1n2" = n1n2)
    return(dfZ)
  } else {
    return(zeros)
  }
}

### Calculate numerical derivative of discrete time-series
calDev <- function(x, y)
{
  if(length(x) != length(y)){stop('x and y dont have the same length')}
  
  dev <- numeric(length(x)-1)
  for(i in 1:length(x))
  {
    if((i + 1) > length(x)){next()}
    
    iDev <- (y[i+1] - y[i]) / (x[i+1] - x[i]) 
    dev[i] <- iDev
  }
  return(dev)
}
