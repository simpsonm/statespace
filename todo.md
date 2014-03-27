Things so do to get interweaving paper ready
==========

This is just for the interweaving paper. I'm putting off the "which covariance prior?" paper for now.

1. Finish simulations with the alternate prior.

2. See if it's possible to unite the two priors in any way. (Looks hard / impossible. Can save for later paper.)

3. Get some intuition about WHY certain DAs work better in some situations instead of others. Will require reading a bunch of papers about DAs, parameterizations, and EM algorithm stuff.

4. Implement AWOL smoothing? Not needed for relative comparisons, though AWOL does work directly with the error sampler, so maybe I should. In which case, rerun all sims. And improve speed of everything else too.

5. Wrongly scaled DAs? Check their properties at least. So get a better sampler.
   Better sampler gotten. ARS for log(VW), and pick a proposal when that fails, again for log(VW)
   Just a worse version of state sampler, it turns out

6. Vivek's idea: scaled by BOTH standard deviations, no? (Probably in a later paper)

7. Real data example? Perhaps with a different model? E.g. dynamic regression so we have to augment F in order to do the scaled errors. Or a vector autoregression (where dim(F) = dim(G), F=I, G is time constant and unknown - VAR(1) anyway). (Maybe don't need this yet - could just wait until the journals ask for more).
   Idea for augmenting F: gramm - schmidt algorithm (assuming the rows of F are already linearly independent, just need more rows!!) 
   (don't need gramm schmidt for dynamic regression - easy, just make the diagonal 1s and 0s everwhere else except the first row is F = [1, x1, ..., xk])

8. Tighten up the paper into a cohesive whole.

9. Scale theta instead of disturbances/errors ? Do a quick check.
   This is the same in the local level model, but not the same whenever G_t is a matrix