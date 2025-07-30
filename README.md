FW over Bergman Fans is an R Project that contains all scripts, data, and figures found in "Tropical Fermat-Weber Points over Bergman Fans" by (Cox et al.) https://arxiv.org/abs/2505.09584.

The project itself contains four folders (Scripts, Data, Results, Figs).

"Scripts" - Contains (4) R script files:
1) Functions.R - Provides a list of various functions used in computing matroidal information, maximal cones of various fan structures, and Fermat-Weber computations. These functions are loaded at the start of every other script.
2) Bergman_Fan.R - Where one defines the matroid M and computes all relevant information (circuits, flats, maximal cones, etc.). The output is saved to "Data" folder.
3) Generate_Samples.R - Generates and saves the samples generated randomly on the Bergman fan. Saves samples into "Data" folder.
4) Stochastic_SR.R - Implements Algorithm 1 and saves the results into the "Results" folder.

"Data" - R data files (.rda) for the matroid and the samples generated on its Bergman Fan.
"Results" - .rda files for the computations involving checks of Lemma 17 and Theorem 18, and the stochastic results described in Section 4.
"Figs" - .png graphics found in the paper.
