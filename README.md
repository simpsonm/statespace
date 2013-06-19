statespace
==========

This repository is a collection of a couple of state space modeling projects originating in an RA with Jarad Niemi in Summer, 2013. Currently there are two main components:

1. An exploration of different sampling methods for linear, gaussian state space models. The ultimate goal is to come up with an interweaving type strategy that (hopefully) converges faster than other approaches.

2. An applied project analyzing data from an ecological experiment exploring the impacts of global warming. The data from this experiment is 65 univariate time series each of around 3000 observations.

Note on .Rnw files
-------------
I use knitr instead of Seave to compile .Rnw files into latex. Some of the syntax is different, so Seave won't work with my .Rnw files. To use knitr, run the following code in R:

    library(knitr)
    knit("filename.Rnw")

This will output filename.tex to the current working directory. Simply run latex on this file to get the pdf output.


Directories
-------------

1. **distsmooth**: Some code for sampling from the system disturbances directly rather than using FFBS then transforming. Note: this samples from the *unscaled* disturbances. Might be useful for decreasing computation time in some algorithms, but not for improving convergence.

2. **exper**: Data from the ecological experiment and some code for manipulating and fitting models to that data. (Note: only what is aloud to be publicly available is here)

   a. **data**: Contains the data along with some documents explaining and analyzing the data.

   b. **replicate**: Contains code attempting to replicate some previous analysis performed on the data.

3. **mcmcex**: A bunch of MCMC examples for the local level model. Constains: 1) code for implementing a few different samplers for the local level model, 2) simulations fitting the local level model to a bunch of fake data, and 3) explanations of when each sampler works well in terms of convergence.

4. **papers**: A variety of useful resources, in pdf form.

5. **random**: Random bits of code that aren't used anywhere but may be useful in the future - from before this was turned into a git repository.