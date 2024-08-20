### Plot time series of all found flecks
plotTSfleck <- function(time, var, zeroes, fleck_data = FALSE, timeSplit = 10)
{
  timeFrames <- length(time) / 100
  timeStep <- time[2] - time[1]
  timeTot <- max(time)
  
  x_out <- seq(0, max(time), timeStep / timeSplit)
  int <- approx(x = time, y = var, xout = x_out)
  dev <- calDev(x = int$x, y = int$y)
  
  for(iTime in 1:timeFrames)
  {
    minTime <- (iTime - 1) * (timeTot / timeFrames)
    maxTime <- (iTime) * (timeTot / timeFrames)
    
    par(mfrow = c(2,1), mar = c(0,4.5,0.2,0.2), oma = c(4,0,0,0), bty = "L")
    plot(-500, ylim = c(min(var,na.rm = T),max(var,na.rm = T)), xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = c(minTime,maxTime))
    abline(v = z * timeStep / timeSplit, col = "gray80", lty = 2)
    points(var~time, type = "o", pch = 20, lwd = 1, col = "gray20")
    if(fleck_data != FALSE){
      points(baseline1~baselineTime1, data = fleck_data, type = "p", pch = 20, col = "chartreuse3")
      points(baseline2~baselineTime2, data = fleck_data, type = "p", pch = 20, col = "chartreuse3")
      points(peak~peakTime, data = fleck_data, type = "p", pch = 20, col = "coral2")
      text(x = fleck_data$peakTime, y = fleck_data$peak, labels = fleck_data$no, pos = 1, font = 2, col = "coral2", cex = 0.8)
    }
    axis(side = 1, labels = FALSE)
    axis(side = 2, las = 2, font = 2, cex.axis = 0.8)
    mtext(outer = F, side = 2, text = "Variable", cex = 1.2, line = 3, font = 2)
    
    plot(-500, ylab = "", ylim = c(min(dev),max(dev)), xlim = c(minTime, maxTime), xaxt = "n", yaxt = "n")
    abline(v = z * timeStep / timeSplit, col = "gray80", lty = 2)
    abline(h = 0, col = "antiquewhite4")
    points(dev~head(seq(0,timeTot,timeStep/timeSplit),-1), type = "l", lwd = 1)
    axis(side = 2, las = 2, font = 2, cex.axis = 0.8)
    axis(side = 1, font = 2)
    mtext(outer = T, side = 1, text = "Time (s)", cex = 1.2, line = 2.2, font = 2)
    mtext(outer = F, side = 2, text = "1st derivative", cex = 1.2, line = 3, font = 2)
  }
}

### Plot time series of all found flecks, simpler version
plotTSfleckEz <- function(time, var, zeroes, fleck_data = FALSE, timeSplit = 10)
{
  timeFrames <- length(time) / 100
  timeStep <- time[2] - time[1]
  timeTot <- max(time)
  
  x_out <- seq(0, max(time), timeStep / timeSplit)
  int <- approx(x = time, y = var, xout = x_out)
  dev <- calDev(x = int$x, y = int$y)
  
  if(length(fleck_data) > 1){ifleck_data = TRUE} else {ifleck_data = FALSE}
  
  for(iTime in 1:timeFrames)
  {
    minTime <- (iTime - 1) * (timeTot / timeFrames)
    maxTime <- (iTime) * (timeTot / timeFrames)
    
    par(mfrow = c(1,1), mar = c(0,4.5,0.2,0.2), oma = c(4,0,0,0), bty = "L")
    plot(-500, ylim = c(min(var,na.rm = T),max(var,na.rm = T)), xaxt = "n", yaxt = "n", xlab = "", ylab = "", xlim = c(minTime,maxTime))
    #abline(v = z * timeStep / timeSplit, col = "gray80", lty = 2)
    points(var~time, type = "o", pch = 20, lwd = 1, col = "gray20")
    if(ifleck_data == TRUE){
      points(baseline1~baselineTime1, data = fleck_data, type = "p", pch = 20, col = "chartreuse3")
      points(baseline2~baselineTime2, data = fleck_data, type = "p", pch = 20, col = "chartreuse3")
      points(peak~peakTime, data = fleck_data, type = "p", pch = 20, col = "coral2")
      text(x = fleck_data$peakTime, y = fleck_data$peak, labels = fleck_data$no, pos = 1, font = 2, col = "coral2", cex = 0.8)
    }
    axis(side = 1, labels = TRUE)
    axis(side = 2, las = 2, font = 2, cex.axis = 0.8)
    mtext(outer = F, side = 2, text = "Variable", cex = 1.2, line = 3, font = 2)
  }
}
