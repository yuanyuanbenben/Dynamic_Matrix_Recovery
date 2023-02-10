## Materials for Dynamic Matrix Recovery

### Overview

<details>
    <summary>Directory <strong><em>Code</em></strong> contains all codes used in the simulation experiments and real data examples </summary>
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
    <summary>Directory <strong><em>Output</em></strong> contains the results and outputs of those simulation experiments and real data examples</summary>
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

