# Two-sample testing with local community depth

This repository contains the code and data to reproduce the results in "Two-sample testing with local community depth" (Evans and Berenhaut). Hypothesis testing is done using the `depthtestr` package, available at [https://github.com/ciaran-evans/depthtestr](https://github.com/ciaran-evans/depthtestr).

## File list

* `depth_plots.R` contains the code for creating Figure 1 in the paper (a comparison between local community depth and Mahalanobis depth on a univariate mixture of Gaussians)
* The `simulations` folder contains everything needed to reproduce the simulations in Section 4 of the paper, and the supplementary materials. Each simulation is included as a `.R` file, with its output included as a `.csv` file.
  * Section 4.1 (Table 1): `normal_location_shift.R` and `normal_location_shift_results.csv`
  * Section 4.2 (Table 2): `normal_scale_shift.R` and `normal_scale_shift_results.csv`
  * Section 4.3 (Table 3): `lognormal_shift.R` and `lognormal_shift_results.csv`
  * Section 4.4 (Table 4): `mixture_scale_shift_2_components.R` and `mixture_results_2_components.csv`
  * Section 4.4 (Figure 2): `mixture_power_vs_sd.R`, `mixture_power_vs_distance.R`, `mixture_results_power_vs_sd.csv`, `mixture_results_power_vs_distance.csv`
  * Section 4.4.1 (Table 5): `mixture_scale_shift_3_components.R` and `mixture_results_3_components.csv`
  * Supplementary materials (comparison with PaLD): the supplementary materials compares local community depth with partitioned local depth [(PaLD)](https://www.pnas.org/doi/10.1073/pnas.2003634119). The power and type I error for PaLD can be found in the same files as the other depth functions
  * Supplementary materials (small sample correction): Section 3.5 of the paper proposes a correction to the Liu-Singh test when the samples are relatively small. The type I error and power for the corrected tests are presented in the main text. The results for the *uncorrected* tests can be found in the supplementary materials, and are included in the `simulations` folders in the files marked `uncorrected`
* The `dengue-analysis` folder contains everything needed to reproduce the real data example in Section 5 of the paper. 
  * `dengue_original.csv`: the original dengue data from [Tuan *et al.* (2015)](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0003638)
  * `dengue_example.R`: the code for performing the full dengue analysis and simulations (Figure 3 and Figure 4)
  * `dengue_power_vs_positives.csv` and `dengue_power_vs_positives_subset.csv`: results used to create Figure 4 (left and right, respectively)
