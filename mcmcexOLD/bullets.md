Bullet Points
=========

This document is a short list of the basic findings of my simulations with the local level model. Throughoout, *V* denotes the observational variance, *W* denotes the system variance, and *R = W/V* denotes the signal-to-noise ratio.


Basic Samplers
----------

### State Sampler
* When *R << 1* the state sampler tends to have high autocorrelation for *W* and low autocorrelation for *V*. This results in poor mixing for *W* though convergence still happens pretty quickly.

* When *R >> 1* the state sampler tends to have high autocorrelation for *V* and low autocorrelation for *W*. This results in both poor mixing and poor convergence for *V*.

* When *R* is near 1, the state sampler tends to have acceptable levels of autocorrelation for both *V* and *W* and no problems with mixing or convergence with either.

### Scaled Disturbance Sampler
* When *R << 1* the scaled disturbance sampler has low autocorrelation, good mixing and good convergence for both *V* and *W*.

* When *R* is near 1 or larger the scaled disturbance sampler has high autocorrelation for both *V* and *W*. *W* exhibits awful mixing and convergence in this situation. *V* exhibits good mixing but slow convergence for *R* near 1 and slow mixing and convergence for *R* larger than 1.

### Scaled Error Scampler
* When *R* is near 1 or larger the scaled error sampler has low autocorrelation, good mixing and good convergence for both *V* and *W*.

* When *R* << 1 the scaled error sampler has high autocorrelation for both *V* and *W*. Both *V* and *W* mix slowly but still converge quickly, at least in the range of the parameter space tested.

### Across all samplers
* Increasing the length of the time series, *T*, exacerates all autocorrelation, mixing, and convergence problems.

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

* Whether you alternate or interweave doesn't appear to affect convergence or mixing when looking at trace plots.

### State-Disturbance Samplers
* The state-disturbance samplers have poor mixing for *V* when *R* > 1, poor convergence when *R* >> 1, and good mixing and covergence everywhere else.

* The state disturbance samplers appear to always have good mixing and convergence for *W*.

### State-Error Samplers
* The state-error samplers have poor mixing for *W* when *R << 1* and good mixing everywhere else. Convergence appears to be good everywhere, at least in the region of the parmaeter space I explored.

* The state-error samplers always have good mixing and convergence for *V*.

### Generally for State-X Samplers
* The state-error samplers appear to work well for *W* more often than the state-disturbance samplers work well for *V*, at least relative to the values of *R* we tend to expect.

* Notably, both the state-error and state-disturbance samplers appear to work well for larger regions of the parameter space than any of the more basic samplers.

### Dist-Error and Triple Samplers
* The dist-error samplers appear to always have good convergence and mixing for both *V* and *W*.

* The triple samplers appear identical to the dist-error samplers.
