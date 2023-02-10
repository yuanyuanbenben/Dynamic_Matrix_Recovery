## Code Part

All R codes used in "*Dynamic Matrix Recovery*" are in this directory.
 
- The folder "**simulation**" contains codes used for "*Section 4: Simulation Studies*"
    - "***DFISTA.R***" and "***baseline_FISTA.R***": functions to implement the optimal algorithm used for methods presented in article under the dynamic matrix completion setting.
    - "***simulations.R***": functions to present simulations including the performance compared to benchmarks, the effect of sample size and number of time points and the dependence across time for corvariance and noise.
    - "***help_functions.R***": some auxiliary functions for above implement.
- The folder "**real_data**" contains codes used for "*Section 5: Real Data Examples*"
    - "***realdata1_preprocess.R***" and "***realdata2_preprocess.R***": functions for the curation and preparation process of the two real dataset.
    - "***realdata1_test.R***" and "***realdata2_test.R***": functions to implement DLR methods and benchmarks in the two real dataset.
    - "***cs_DFISTA.R***" and "***cs_baseline_FISTA.R***": functions to implement the optimal algorithm used for methods presented in article under the dynamic compressed sensing setting.
    - "***robust_pca.R***": functions used for the preprocess of the vedio example which seperates a matrix into a low rank and a sparse part.
    - "***help_functions.R***": some auxiliary functions for above implement.
- "***plot.R***": functions to collect results in simulations and real data examples and draw the pictures used in the article.
