# R Algorithm to find sunflecks in irradiance time-series  

## Description of functions

### findzeroes
First function to use on discrete time-series. Returns a vector giving the position of every zero crossing within the time-series. For every timepoints, the numerical derivative is calculated using `calDev` and each sequential timepoints of the numerical derivative are multiplied. When the results is negative, the time-series crosses zero. 

Arguments:
**time** [no defaults]: Vector of the time of the time-series (x-axis)
**var** [no defaults]:  Variable of interested linked with time of the time-series (y-axis)
**lim** [default = 0.0005]: Limit for multiplication, only values higher than lim are kept (removes noise and prevents recording zeroes when time-series is flat.)
**timeSplit** [default = 10]: Increase time-series frequency (value at t = 0 are copied from t = 1 to t = 9, etc.). 10 usually garantees accuracy.
**return_n1n2** [default = FALSE]:  If true, also return the value of the multiplication

### findflecks
### plotflecks
### combineInteg

