## Materials for Dynamic Matrix Recovery

### Overview

<details>
    <summary>Directory <strong><em>code</em></strong> contains all codes used in the simulation experiments and real data examples </summary>
    <ul>
        <li>The folder <strong><em>simulation</em></strong> contains codes used for "<em>Section 4: Simulation Studies</em>"
        <ul>
            <li>"<strong><em>DFISTA.R</em></strong>" and "<strong><em>baseline_FISTA.R</em></strong>": functions to implement the optimal algorithm used for methods presented in article under the dynamic matrix completion setting.</li>
            <li> "<strong><em>simulations.R</em></strong>": functions to present simulations including the performance compared to benchmarks, the effect of sample size and number of time points and the dependence across time for corvariance and noise.</li>
            <li> "<strong><em>help_functions.R</em></strong>": some auxiliary functions for above implement.</li>
        </ul>
        </li>
        <li> The folder <strong><em>real_data</em></strong> contains codes used for "<em>Section 5: Real Data Examples</em>"
        <ul>
            <li>"<strong><em>realdata1_preprocess.R</em></strong>" and "<strong><em>realdata2_preprocess.R</em></strong>": functions for the curation and preparation process of the two real dataset.</li>
            <li>"<strong><em>realdata1_test.R</em></strong>" and "<strong><em>realdata2_test.R</em></strong>": functions to implement DLR methods and benchmarks in the two real dataset.</li>
            <li>"<strong><em>cs_DFISTA.R</em></strong>" and "<strong><em>cs_baseline_FISTA.R</em></strong>": functions to implement the optimal algorithm used for methods presented in article under the dynamic compressed sensing setting.</li>
            <li>"<strong><em>robust_pca.R</em></strong>": functions used for the preprocess of the vedio example which seperates a matrix into a low rank and a sparse part.</li>
            <li> "<strong><em>help_functions.R</em></strong>": some auxiliary functions for above implement.</li>
        </ul>
        </li>
        <li>"<strong><em>plot.R</em></strong>": functions to collect results in simulations and real data examples and draw the pictures used in the article.</li>
    </ul>
</details>

<details>
    <summary>Directory <strong><em>data</em></strong> contains original real data and environment after preprocessing used in main context.</summary>
    <ul>
        <li>"<strong><em>netflix_data.zip</em></strong>": the raw netflix dataset.</li>
        <li>"<strong><em>lions/0000(0095).jpg</em></strong>": the raw 1-96 frames of the lions video.</li>
        <li>"<strong><em>netflix_data.RData</em></strong>": the working environment after preprocess the netflix dataset.</li>
        <li>"<strong><em>video_data.RData</em></strong>": the working environment after process the lions video dataset.</li>
    </ul>
</details>


<details>
    <summary>Directory <strong><em>manuscript</em></strong> contains figures generated by numerical simulation and real data experiments used in the main context.</summary>
    <ul>
        <li>"<strong><em>baseline.png</em></strong>": The comparisons of MSE for DLR and three benchmarks, which corresponds to <em>Figure 1</em>
        in the main context.</li>
        <li>"<strong><em>phase_transition.png</em></strong>": The MSE of DLR estimates under different settings of sample size and number of time point, which corresponds to <em>Figure 2 left</em> in the main context.</li>
        <li>"<strong><em>precise_nt.png</em></strong>": The relationship between the logarithm of AMSE and the logarithm of $\rho/\tau$ show a linear trend with slope -4/5, which corresponds to <em>Figure 2 right</em> in the main context.</li>
        <li>"<strong><em>dependent_xi.png</em></strong>": Under the dependent noise case, the influence of different noise variances and dependent coefficients $\beta$ to MSE, which corresponds to <em>Figure 3 left</em> in the main context.</li>
        <li>"<strong><em>dependent_whit_phiy.png</em></strong>": The relationship between AMSE and the mixing coefficient of noise $\xi$, which corresponds to <em>Figure 3 right</em> in the main context.</li>
        <li> "<strong><em>dependent_X.png</em></strong>": Under the dependent noise case, the influence of different noise variances and dependent coefficients $\alpha$ to MSE, which corresponds to <em>Figure 4 left</em> in the main context.</li>
        <li>"<strong><em>dependent_whit_phix.png</em></strong>": The relationship between AMSE and the mixing coefficient of covariance matrix X, which corresponds to <em>Figure 4 right</em> in the main context.</li>
        <li>"<strong><em>netflix_data.png</em></strong>": The MSE at differnet time points for DLR and three benchmarks, which corresponds to <em>Figure 5</em> in the main context.</li>
        <li>"<strong><em>real_lions_5(25,45,65,85).png</em></strong>, <strong><em>lions_5(25,45,65,85).png</em></strong>, <strong><em>baseline_lions_5(25,45,65,85).png</em></strong>, <strong><em>twostep_lions_5(25,45,65,85).png</em></strong>": The original, DLR estimated, Static estimated and TwoStep estimated frames in the lions video, which correspond to <em>Figure 6</em> in the main context.</li>
    </ul>
</details>

<details>
    <summary>Directory <strong><em>output</em></strong> contains the results and outputs of those simulation experiments and real data examples</summary>
    <ul>
        <li><strong>Simulation</strong>
        <ul>
            <li><strong>independent case</strong>
            <ul>
                <li>"<strong><em>dmc_5000_30000.csv</em></strong>": the output MSE of simulation 1 using our proposed DLR method with T=100 and sample size from 5000 to 30000.</li>
                <li> "<strong><em>baseline_120000.csv</em></strong>": the output MSE of simulation 1 using Static benchmark with T=100 and sample size 120000.</li>
                <li> "<strong><em>local_smoth_120000.csv</em></strong>": the output MSE of simulation 1 using TwoStep benchmark with T=100 and sample size 120000.</li>
                <li> "<strong><em>tensor_30000.csv</em></strong>": the output MSE of simulation 1 using Tensor benchmark with T=100 and sample size 30000.</li>
                <li> "<strong><em>phase_transition.csv, phase_transition_precise.csv</em></strong>": the output MSE of simulation 1 using DLR method with different settings of the number of time points and sample size.</li>
            </ul>
            </li>
            <li><strong>dependent case</strong>
            <ul>
                <li>"<strong><em>dependent_mc.csv</em></strong>": the output MSE of simulation2 under the setting that noise is dependent while covariance X is independent.</li>
                <li>"<strong><em>dependent_X_mc.csv</em></strong>": the output MSE of simulation2 under the setting that covariance X is dependent while noise is independent.</li>
            </ul>
            </li>
        </ul>
        </li>
        <li><strong>Real data example</strong>
        <ul>
            <li><strong>Netflix dataset</strong>
            <ul>
                <li>"<strong><em>netflix/mse_1(100).csv</em></strong>": the output MSE of netflix data example using our DLR method.</li>
                <li>"<strong><em>netflix/baseline_mse_1(100).csv</em></strong>": the output MSE of netflix data example using Static method. </li>
                <li>"<strong><em>netflix/baseline_mse_twostep.csv</em></strong>": the output MSE of netflix data example using TwoStep method. </li>
                <li>"<strong><em>netflix/baseline_mse_tensor.csv</em></strong>": the output MSE of netflix data example using Tensor method. </li>
            </ul>
            </li>
            <li><strong>Davis 2016 lions video</strong>
            <ul>
                <li>"<strong><em>lions_ren(blue,green)_5(25,45,65,85).csv</em></strong>": the output rgb values of corresponding fames using DLR method. </li>
                <li>"<strong><em>baseline_lions_ren(blue,green)_5(25,45,65,85).csv</em></strong>": the output rgb values of corresponding famres using Static method. </li>
            </ul>
            </li>
        </ul>
        </li>
    </ul>
</details>

### Workflows

<details>
    <summary><strong>Simulation 1 and 2</strong></summary>
    <ol>
        <li>Open the floder <em>Dynamic_Matrix_Recovery/code/simulation/</em>.</li>
        <li>Make sure that the R scripts <strong><em>DFISTA.R</em></strong>, <strong><em>baseline_FISTA.R</em></strong> and <strong><em>help_functions.R</em></strong> are available in the path <em>~/Dynamic_Matrix_Recovery/code/simulation/...</em>. Or revise the path in <strong><em>simulations.R</em></strong> line 15-17 to match the path of those three R scripts.</li>
        <li>In the R scripts <strong><em>simulations.R</em></strong>simulations.R, independent case and dependent case simulations are performed. Running corresponding codes of different methods including our DLR method and three benchmarks Static, TwoStep and Tensor can obtain corresponding results. Parameters including the dimension of matrix, number of time points and sample size can easily changed if needed. All codes that stores output is commented out to avoid overwriting existing output.</li>
        <li> For <em>Figure 1-4</em>, return to the previous directory and run the corresponding part in <strong><em>plot.R</em></strong>.</li>
    </ol>
</details>
<!-- #### Simulation 1 and 2

1. Open the floder *Dynamic_Matrix_Recovery/code/simulation/*.
2. Make sure that the R scripts ***DFISTA.R***, ***baseline_FISTA.R*** and ***help_functions.R*** are available in the path *~/Dynamic_Matrix_Recovery/code/simulation/...*. Or revise the path in ***simulations.R*** line 15-17 to match the path of those three R scripts.
3. In the R scripts ***simulations.R***, independent case and dependent case simulations are performed. Running corresponding codes of different methods including our DLR method and three benchmarks Static, TwoStep and Tensor can obtain corresponding results. Parameters including the dimension of matrix, number of time points and sample size can easily changed if needed. All codes that stores output is commented out to avoid overwriting existing output.
4. For *Figure 1-4*, return to the previous directory and run the corresponding part in ***plot.R***. -->

<!-- #### Real data example 1: the netflix dataset

1. Open the floder *Dynamic_Matrix_Recovery/code/real_data/*.
2. Make sure the R scripts ***DFISTA.R***, ***baseline_FISTA.R*** and ***help_functions.R*** are available in the path *~/Dynamic_Matrix_Recovery/code/simulation/...* and the R script ***help_functions.R*** is available in the path *~/Dynamic_Matrix_Recovery/code/real_data/...*. Or revise the path in ***realdata1_test.R*** line 15-18 to match the path of those three R scripts.
3. Download the ***netflix_data.zip*** in the directory *~/Dynamic_Matrix_recovery/data/*, unzip it and save it as *~/Dynamic_Matrix_recovery/data/netflix_data.csv*(Or download the dataset [here](https://www.kaggle.com/datasets/netflix-inc/netflix-prize-data)). Run the R script ***realdata1_preprocess.R*** to preprocess the dataset. Another alternative is directly using the environment ***netflix_data.RData*** in the directory *Dynamic_Matrix_Recovery/data/*.
4. Run the R script ***realdata1_test.R***, in which the DLR method and benchmarks can be applied to the dataset and codes that stores output is commented out to avoid overwriting existing output. 
5. For *Figure 5*, return to the previous directory and run the corresponding part in ***plot.R***.

 -->
 

<!-- #### Real data example 2: the lion video in Davis 2017 dataset

1. Open the floder *Dynamic_Matrix_Recovery/code/real_data/*.
2. Make sure the R scripts ***cs_DFISTA.R***, ***cs_baseline_FISTA.R***, ***help_functions.R*** and ***robust_pca.R*** are available in the path *~/Dynamic_Matrix_Recovery/code/real_data/...*. Or revise the path in ***simulations.R*** line 13-16 to match the path of those four R scripts.
3. Download the ***0000.jpg*** to ***0095.jpg*** in the directory *~/Dynamic_Matrix_recovery/data/lions*(Or download the dataset [here](https://davischallenge.org/davis2017/code.html#unsupervised)). Run the R script ***realdata2_preprocess.R*** to preprocess the dataset. Another alternative is directly using the environment ***vedio_data.RData*** in the directory *Dynamic_Matrix_Recovery/data/*.
4. Run the R script ***realdata2_test.R***, in which the DLR method and benchmarks can be applied to the dataset and codes that stores output is commented out to avoid overwriting existing output.
5. For *Figure 6*, return to the previous directory and run the corresponding part in ***plot.R***.
-->



 <details>
    <summary><strong>Real data example 1: the netflix dataset</strong></summary>
    <ol>
        <li>Open the floder <em>Dynamic_Matrix_Recovery/code/real_data/</em>.</li>
        <li>Make sure the R scripts <strong><em>DFISTA.R</em></strong>, <strong><em>baseline_FISTA.R</em></strong> and <strong><em>help_functions.R</em></strong> are available in the path <em>~/Dynamic_Matrix_Recovery/code/simulation/...</em> and the R script <strong><em>help_functions.R</em></strong> is available in the path <em>~/Dynamic_Matrix_Recovery/code/real_data/...</em>. Or revise the path in <strong><em>realdata1_test.R</em></strong> line 15-18 to match the path of those four R scripts.</li>
        <li>Download the <strong><em>netflix_data.zip</em></strong> in the directory <em>~/Dynamic_Matrix_recovery/data/</em>, unzip it and save it as <em>~/Dynamic_Matrix_recovery/data/netflix_data.csv</em>(or download the dataset <a href="https://www.kaggle.com/datasets/netflix-inc/netflix-prize-data">here</a>). Run the R script <strong><em>realdata1_preprocess.R</em></strong> to preprocess the dataset. Another alternative is directly using the environment <strong><em>netflix_data.RData</em></strong> in the directory <em>Dynamic_Matrix_Recovery/data/</em>.</li>
        <li> Run the R script <strong><em>realdata1_test.R</em></strong>, in which the DLR method and benchmarks can be applied to the dataset and codes that stores output is commented out to avoid overwriting existing output.</li>
        <li>For <em>Figure 5</em>, return to the previous directory and run the corresponding part in <strong><em>plot.R</em></strong>.</li>
    </ol>
</details>
<details>
    <summary><strong>Real data example 2: the lion video in Davis 2017 dataset</strong></summary>
    <ol>
        <li>Open the floder <em>Dynamic_Matrix_Recovery/code/real_data/</em>.</li>
        <li>Make sure the R scripts <strong><em>cs_DFISTA.R</em></strong>, <strong><em>cs_baseline_FISTA.R</em></strong>, <strong><em>help_functions.R</em></strong> and <strong><em>robust_pca.R</em></strong> are available in the path <em>~/Dynamic_Matrix_Recovery/code/real_data/...</em>. Or revise the path in <strong><em>simulations.R</em></strong> line 13-16 to match the path of those four R scripts.</li>
        <li>Download the <strong><em>0000.jpg</em></strong> to <strong><em>0095.jpg</em></strong> in the directory <em>~/Dynamic_Matrix_recovery/data/lions</em>(or download the dataset <a href="https://davischallenge.org/davis2017/code.html#unsupervised">here</a>). Run the R script <strong><em>realdata2_preprocess.R</em></strong> to preprocess the dataset. Another alternative is directly using the environment <strong><em>vedio_data.RData</em></strong> in the directory <em>Dynamic_Matrix_Recovery/data/</em>.
        <li> Run the R script <strong><em>realdata2_test.R</em></strong>, in which the DLR method and benchmarks can be applied to the dataset and codes that stores output is commented out to avoid overwriting existing output.</li>
        <li>For <em>Figure 6</em>, return to the previous directory and run the corresponding part in <strong><em>plot.R</em></strong>.</li>
    </ol>
</details>
