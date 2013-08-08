Bullet Points
=========

This document is a short list of the basic findings of my simulations with the local level model. Throughoout, *V* denotes the observational variance, *W* denotes the system variance, and *R = W/V* denotes the signal-to-noise ratio.


Basic Samplers
----------

### State Sampler
* When *R << 1* the state sampler tends to have high autocorrelation for *W* and low autocorrelation for *V*. This results in poor mixing for *W* though convergence still happens pretty quickly.

* When *R >> 1* the state sampler tends to have high autocorrelation for *V* and low autocorrelation for *W*. This results in both poor mixing and poor convergence for *V*.

* When *R* is near 1, the state sampler tends to have acceptable levels of autocorrelation for both *V* and *W* and no problems with mixing or convergence with either.

* Another way to think about the state sampler is that when *V* and *W* are significantly different from each other, the state sampler has trouble with the smaller one.

### Scaled Disturbance Sampler
* When *R << 1* the scaled disturbance sampler has low autocorrelation, good mixing and good convergence for both *V* and *W*.

* When *R* is near 1 or larger the scaled disturbance sampler has high autocorrelation for both *V* and *W*. *W* exhibits awful mixing and convergence in this situation. *V* exhibits good mixing but slow convergence for *R* near 1 and slow mixing and convergence for *R* larger than 1.

* The quick and dirty way to remember this is that the scaled disturbance sampler has trouble with both *V* and *W* when *V* is small (note *V* is the variance parameter that is NOT involved in the transformation of the states).

### Scaled Error Scampler
* When *R* is near 1 or larger the scaled error sampler has low autocorrelation, good mixing and good convergence for both *V* and *W*.

* When *R* << 1 the scaled error sampler has high autocorrelation for both *V* and *W*. Both *V* and *W* mix slowly but still converge quickly, at least in the range of the parameter space tested.

* The quick and dirty way to remember this is that the scaled error sampler has trouble with both *V* and *W* when *W* is small (note *W* is the variance parameter that is NOT involved in the transformation of the states).


### Across all samplers
* Increasing the length of the time series, *T*, exacerbates all autocorrelation, mixing, and convergence problems.

* Changing the values of *V* and *W* doesn't seem to affect convergence unless it changes the signal to noise ratio, *R = W/V*.

Posterior Correlation
------------

### Between Variance Parameters
* The posterior correlation between *V* and *W* tends to be highest when *R* is near 1, and is typically less than 0.5.

### Between Variance Parameters and States
* The posterior correlation between *V* and the states tends to be highest when *R* is near 1, and is typically less than 0.5.

* The posterior correlation between *W* and the states tends to be highest when *R* is near 1, and is typically less than 0.5.

### Between Variance Parameters and Scaled Disturbances
* The posterior correlation between *V* and the scaled disturbances tends to be highest when *R* is near 1, and is typically less than 0.5.

* The posterior correlation between *W* and the scaled disturbances tends to be highest when *R >> 1*, and is as high as 0.99 or so.

### Between Variance Parameters and Scaled Errors
* The posterior correlation between *V* and the scaled errors tends to be highest when *R >> 1*, and is as high as 0.5 or so.

* The posterior correlation between *W* and the scaled errors tends to be highest when *R* is near 1, and is as high as 0.5 or so.

Interweaving and Alternating Samplers
----------------------

* Whether you alternate or interweave doesn't appear to affect convergence or mixing when looking at trace plots or first order autocorrelation.

### State-Disturbance Samplers
* The state-disturbance samplers have poor mixing for *V* when *R* > 1, poor convergence when *R* >> 1, and good mixing and covergence everywhere else.

* The state disturbance samplers appear to always have good mixing and convergence for *W*.

### State-Error Samplers
* The state-error samplers have poor mixing for *W* when *R << 1* and good mixing everywhere else. Convergence appears to be good everywhere, at least in the region of the parmaeter space I explored.

* The state-error samplers always have good mixing and convergence for *V*.

### Generally for State-X Samplers
* A state-X sampler seems to only have problems when both the state sampler and the X sampler has problems, as expected. E.g. when *V* is small relative to *W*, the state sampler will have trouble with *V* and the scaled-disturbance sampler will have trouble with both *V* and *W*, but the state-disturbance samplers (which alternative or interweave between the state and scaled-disturbance samplers) will only have trouble with *V*. When *V* is large relative to *W*, the state sampler will have trouble with *W* while the scaled-disturbance sampler won't have trouble with *V* or *W*. In this case, the state-disturbance sampler won't have trouble with *V* or *W*.

* So generally, the state-error samplers only ever have trouble for *W* and the state-disturbance samplers only ever have trouble for *V*.

* The state-error samplers appear to work well for *W* more often than the state-disturbance samplers work well for *V*, at least relative to the values of *R* we tend to expect.

* Notably, both the state-error and state-disturbance samplers appear to work well for larger regions of the parameter space than any of the more basic samplers they are composed of.

### Dist-Error and Triple Samplers
* The dist-error samplers appear to always have good convergence and mixing for both *V* and *W*, except with *R* is significantly different from 1 (indicating that *V* and *W* are significantly different) the smaller of *V* and *W* will have worse mixing. I.e. when *V << W* *V* will have relatively poor mixing, and when *V >> W* *W* will have relatively poor mixing. Convergence, however, appears to always be good - at least in the region of the parameter space I looked at.

* The triple samplers appear identical to the dist-error samplers (except, of course, computation time!)
