Things so do to get interweaving paper ready
==========

This is just for the interweaving paper. I'm putting off the "which covariance prior?" paper for now.

1. Finish simulations with the alternate prior.

2. See if it's possible to unite the two priors in any way. (Looks hard / impossible. Can save for later paper.)

3. Get some intuition about WHY certain DAs work better in some situations instead of others. Will require reading a bunch of papers about DAs, parameterizations, and EM algorithm stuff.

4. Implement AWOL smoothing? Not needed for relative comparisons, though AWOL does work directly with the error sampler, so maybe I should. In which case, rerun all sims. 

5. Wrongly scaled DAs? Check their properties at least. So get a better sampler.

6. Vivek's idea: scaled by BOTH standard deviations, no? (Probably in a later paper)

7. Real data example? Perhaps with a different model? E.g. dynamic regression so we have to augment F in order to do the scaled errors. Or a vector autoregression (where dim(F) = dim(G), F=I, G is time constant and unknown - VAR(1) anyway). (Maybe don't need this yet - could just wait until the journals ask for more).

8. Tighten up the paper into a cohesive whole.