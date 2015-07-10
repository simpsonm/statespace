This zip file contains the following files:

mcfa.R contains R code for the MCFA smoother in the local level model. 

dlmasisfun.R contains R functions for obtaining simulations from a variety of MCMC algorithms in the local level model, summarizing those simulations, and for simulating data from the model. It depends on wscalerej.R and mcfa.R.

wscalerej.R contains several R functions for obtaining simulations from the wrongly-scaled DA algorithms in the local level model. 

dlmasismixrun.R contains R code that uses the functions in dlmasisfun.R and mcfa.R to obtain the simulations discussed in Section 6. Note: running this script will take a very long time -- on the order of weeks.

dlmasislongrun.R contains R code that uses the functions in dlmasisfun.R to obtain simulations for long time series in the local level model used in Appendix N. 

plots.R contains R code that uses output from dlmasismixrun.R to create the plots in the main document and the appendices.
