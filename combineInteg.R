### Combine time-series from higher resolution to lower resolution (e.g. 10ms to 30ms)
combineInteg <- function(time = df$Time, var = df$PAR_q, newInteg_ms = 20)
{
  timeStep <- time[2] - time[1]
  timeFrame <- tail(time, 1)
  
  # Return basic data frame is new int is same as old
  integ <- newInteg_ms / 1e3
  if(timeStep == integ){
    newdf <- data.frame("Time" = time, "Var" = var)
    newdf$Time <- newdf$Time - (newdf$Time[2] - newdf$Time[1])
    return(newdf)
  }
  
  # Check that original integ time is a multiple of new integration 
  timeRatio <- integ / timeStep
  newLength <- trunc(length(time) / timeRatio)
  if({timeRatio == round(timeRatio)} == FALSE){
    stop("Original integration time is not a multiple of new integration time...")
  }
  
  newTime <- (1:newLength)*timeRatio
  newdf <- data.frame("Time" = numeric(newLength), "Var" = numeric(newLength))
  for(i in newTime)
  {
    it <- timeRatio - 1
    iTime <- (i / 100) - integ
    iVar <- round(mean(var[(i-it):i], na.rm = TRUE),5)
    
    row <- c(iTime, iVar)
    newdf[(i / timeRatio),] <- row
  }
  return(newdf)
}
