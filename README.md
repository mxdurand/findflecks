# R Algorithm to find sunflecks in irradiance time-series  

## Description of functions

### findzeroes
First function to use on discrete time-series. Returns a vector giving the position of every zero crossing within the time-series. For every timepoints, the numerical derivative is calculated using `calDev` and each sequential timepoints of the numerical derivative are multiplied. When the results is negative, the time-series crosses zero. 

_Arguments:_  
- **time** [no defaults]: Vector of the time of the time-series (x-axis).   
- **var** [no defaults]:  Variable of interested linked with time of the time-series (y-axis).   
- **lim** [default = 0.0005]: Limit for multiplication, only values higher than lim are kept (removes noise and prevents recording zeroes when time-series is flat.).   
- **timeSplit** [default = 10]: Increase time-series frequency (value at t = 0 are copied from t = 1 to t = 9, etc.). 10 usually garantees accuracy.  
- **return_n1n2** [default = FALSE]:  If true, also return the value of the multiplication.  

### findflecks
Detects flecks from time series and zero vector given by `findzeroes`. A sunfleck is caracterized by an increased followed by a decrease. Function search zero vector for such events. Once found, checks for asymetry between baselines. If found, tries to extend baseline a bit. If still asymeric, behaviour is as defined by `asmMethod`. Then fleck is checks for criteria as defined by `minTime` `minAmp` `minPdiff`. If passed, the baselines are trimmed. Conditions are checked once more, and trimming is reversed if conditions are not passed anymore. Finally, log the fleck in a table returned at the end. At the end, remove flecks overlapping, and calcualte interval in time between two flecks. 

_Arguments:_  
- **time** [no defaults]: Vector of the time of the time-series (x-axis).  
- **var** [no defaults]:  Variable of interested linked with time of the time-series (y-axis).  
- **zeroes** [no defaults]: Vector of zeroes from `findzeroes` function.  
- **minTime** [default = 0]: Flecks found with duration below `minTime` will be discarded. Keep at 0 to keep all flecks.  
- **minAmp** [default = 0]: Flecks found with amplitude below `minAmp` will be discarded. Keep at 0 to keep all flecks.  
- **minPdiff** [default = 0]: Flecks found with a percent difference between peak and baseline that is below `minPdiff` will be discarded. Keep at 0 to keep all flecks.  
- **asymmetry** [default = 1/4]: Threshold value to qualify fleck as asymeric. Asymetry happens when the baseline have large difference in their value.  
- **trimCV** [default = 0.05]: Control trimming threshold. Fleck baselines are trimmed iteratively based on the coefficient of variation between two points  at each baseline side. 
- **asmMethod** [no default]: One of ["mean", "max", "rm"]. Decides what to do with asymetric flecks `mean` averages the two baselines, `max` keeps the largest baseline, `rm` discard the asymetric fleck.   
- **bounds** [default = c(0,1)]: For relative amplitude calculations,  normalize between 0 and 1 by default.  
- **timeSplit** [default = 10]: Increase time-series frequency (value at t = 0 are copied from t = 1 to t = 9, etc.). 10 usually garantees accuracy.  
- **shadeflecks** [default = FALSE]: If true, run the function in shadefleck mode, ie. the function will find troughs instead of peaks in the data.  
- **verbose** [default = TRUE]: If true, provides more information while running.    

### plotflecks
Plot time series of all found flecks.  

_Arguments:_  
- **time** [no defaults]: Vector of the time of the time-series (x-axis).  
- **var** [no defaults]:  Variable of interested linked with time of the time-series (y-axis).  
- **zeroes** [no defaults]: Vector of zeroes from `findzeroes` function.  
- **fleck_data** [no defaults]:  Returned data frame from  `findflecks` function.  
- **timeSplit** [default = 10]: Increase time-series frequency (value at t = 0 are copied from t = 1 to t = 9, etc.). 10 usually garantees accuracy.  

### combineInteg  
Combine time-series from higher resolution to lower resolution (e.g. 10ms to 30ms). Higher resolution needs to be a multiple of lower resolution.  

_Arguments:_  
- **time** [no defaults]: Vector of the time of the time-series (x-axis).  
- **var** [no defaults]:  Variable of interested linked with time of the time-series (y-axis).  
- **newInteg_ms** [default = 20]:  New integration time for the returned time-series. Use same units as main time-series.    



