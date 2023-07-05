## Code Part


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
            <li>The folder <strong><em>netflix_data</em></strong> contains flies about the Netflix data experiment.
	        <ul>
		        <li>"<strong><em>realdata1_preprocess.R</em></strong>" and "<strong><em>realdata1_test.R</em></strong>": functions for the preparation and test processes of the real dataset.</li>
		        <li>"<strong><em>netflix_DFISTA.R</em></strong>", "<strong><em>netflix_DFISTA_link.R</em></strong>" and "<strong><em>netflix_baseline_FISTA.R</em></strong>": functions to implement the optimal algorithm used for methods presented in article for Netflix dataset.</li>
	        </ul>
            </li>
            <li>The folder <strong><em>video_data</em></strong> contains flies about the lions video in Davis 2017 data experiment.
	        <ul>
		       <li>"<strong><em>realdata2_preprocess.R</em></strong>" and "<strong><em>realdata2_test.R</em></strong>": functions for the preparation and test processes of the real dataset.</li>
		        <li>"<strong><em>cs_DFISTA.R</em></strong>" and "<strong><em>cs_baseline_FISTA.R</em></strong>": functions to implement the optimal algorithm used for methods presented in article under the dynamic compressed sensing setting.</li>
		        <li>"<strong><em>robust_pca.R</em></strong>": functions used for the preprocess of the vedio example which seperates a matrix into a low rank and a sparse part.</li>
	        </ul>
            </li>
            <li>The folder <strong><em>cifar10</em></strong> contains flies about the Cifar10 data experiment.
	        <ul>
		        <li>"<strong><em>compress.py</em></strong>": the preparetion and compressed prcesses for real data.</li>
		        <li>"<strong><em>lenet.py</em></strong>" and "<strong><em>resnet.py</em></strong>": network model used.</li>
		        <li>"<strong><em>utils.py</em></strong>":some auxiliary functions.</li>
		        <li>"<strong><em>main.py</em></strong>": main function.</li>
	        </ul>
            </li>
            <li> "<strong><em>help_functions.R</em></strong>": some auxiliary functions for above implement.</li>   
        </ul>
        </li>
        <li>"<strong><em>plot.R</em></strong>": functions to collect results in simulations and real data examples and draw the pictures used in the article.</li>
    </ul>
</details>
            
