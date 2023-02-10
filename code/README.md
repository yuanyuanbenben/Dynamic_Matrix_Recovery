## Code Part

All R codes used in "*Dynamic Matrix Recovery*" are in this dictionary.

Files in this dictionary implement the following functions. 
- The folder "**simulation**" contains codes used for "*Section 4: Simulation Studies*"
    - The files "**DFISTA.R**" and "**baseline_FISTA.R**" implement the optimal algorithm used for methods presented in article under the dynamic matrix completion setting.
    - The file "**simulations.R**" presents simulations including the performance compared to benchmarks, the effect of sample size and number of time points and the dependence across time for corvariance and noise.
    - The file "**help_functions.R**" contains some auxiliary functions for above implement.
- The folder "**real_data**" contains codes used for "*Section 5: Real Data Examples*"
    - The files "**realdata1_preprocess.R**" and "**realdata2_preprocess.R**" are the curation and preparation process for the two real dataset.
    - The files "**realdata1_test.R**" and "** realdata2_test.R**" implement DLR methods and benchmarks in the two dataset.
    - The files "**cs_DFISTA.R**" and "**cs_baseline_FISTA.R**" implement the optimal algorithm used for methods presented in article under the dynamic compressed sensing setting.
    - The flie "**robust_pca.R**" is used for the preprocess of the vedio example which seperates a matrix into a low rank and a sparse part.
    - The file "**help_functions.R**" contains some auxiliary functions for above implement.
- The file "**plot.R**" collects results in simulations and real data examples and draws the pictures used in the article.
