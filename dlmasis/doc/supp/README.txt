This zip file contains the supplementary materials to the paper Interweaving Markov Chain Monte Carlo Strategies for Efficient Estimation of Dynamic Linear Models:

appendices.pdf: This file provides 7 appendices cited in the main article:
 Appendix A: Proof of lemma 1.
 Appendix B: Derivations of joint distributions and full conditional distributions for each DA in the DLM.
 Appendix C: Mixed Cholesky factorization algorithm (MCFA) for simulation smoothing.
 Appendix D: Further augmentation for non-invertible $F_t$.
 Appendix E: Efficiently drawing from $p(W|V,\gamma,y)$ and $p(V|W,\psi,y)$ in the LLM.
 Appendix F: Efficiently drawing from $p(W|V,\tilde{\gamma},y)$ and $p(V|W,\tilde{\psi},y)$ in the LLM.
 Appendix G: Using posterior correlations to understand patterns of ESP.
 Appendix H: Further plots of ESP and log time per 1000 effective draws for different values of $T$.

mcfa.R: This file contains R code to perform the MCFA for obtaining a draw from $p(\theta|V,W,y)$ and from $p(\psi|V,W,y)$ in the local level model.

dlmasisfun.R: This file contains several R functions for obtaining the simulations from a variety of MCMC algorithms in the local level model, summarizing those simulations appropriately, and for simulating data from the model.

dlmasismixrun.R: This file contains R code that uses the functions in \verb0dlmasisfun.R0 and \verb0mcfa.R0 to obtain the simulations used in the local level model example from Section \ref{sec:LLM}.

plots.R: This file contains R code that uses the output from \verb0dlmasismixrun.R0 to create the various plots in the main document and the appendices.
