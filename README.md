# Two-sample testing with local community depth

This repository contains the code and data to reproduce the results in "Two-sample testing with local community depth" (Evans and Berenhaut). Hypothesis testing is done using the `depthtestr` package, available at [https://github.com/ciaran-evans/depthtestr](https://github.com/ciaran-evans/depthtestr).

## File list

* The `examples` folder contains everything needed to reproduce the two examples from Section 3 (Section 3.3 and Section 3.4, Figures 1 and 2 respectively) 
  * `example_1.R` reproduces Figure 1 (a comparison between local community depth, density as depth, Mahalanobis depth, and the local depth from Paindevaine and van Bever)
  * `example_2.R` reproduces Figure 2 (a comparison between local community depth and an unweighted pairwise depth)
* The `simulations` folder contains everything needed to reproduce the simulations in Section 4 of the paper, and the supplementary materials. Each simulation is included as a `.R` file, with its output included as a `.csv` file.
  * Section 4.1 (Table 1): `lognormal_shift.R` and `lognormal_shift_results.csv`
  * Section 4.2 (Table 2): `mixture_scale_shift_2_components.R` and `mixture_results_2_components.csv`
  * Section 4.2 (Figure 3): `mixture_power_vs_sd.R`, `mixture_power_vs_distance.R`, `mixture_results_power_vs_sd.csv`, `mixture_results_power_vs_distance.csv`
  * Section 4.3 (Table 3): `mixture_scale_shift_3_components.R` and `mixture_results_3_components.csv`
  * Section 4.4 (Table 4): `uniform_ball_scale_shift_3_components.R` and `uniform_ball_scale_shift_3_components_results.csv`
  * Appendix (Table 5): `normal_location_shift.R` and `normal_location_shift_results.csv`
  * Appendix (Table 6): `normal_scale_shift.R` and `normal_scale_shift_results.csv`
  * Appendix (Table 7): `uniform_ball_scale_shift.R` and `uniform_ball_scale_shift_results.csv`
  * Appendix A.2 (comparison with LCD variant): the supplementary materials compares local community depth with the local community depth variant. The power and type I error for the variant can be found in the same files as the other depth functions. Figure 6 (comparison of type I error and power between LCD and the variant) is created by `lcd_vs_pald.R`
  * Appendix A.3 (small sample correction): Section 3.7 of the paper proposes a correction to the Liu-Singh test when the samples are relatively small. The type I error and power for the corrected tests are presented in the main text. The results for the *uncorrected* tests can be found in the supplementary materials, and are included in the `simulations` folders in the files marked `uncorrected`
* The `dengue-analysis` folder contains everything needed to reproduce the real data example in Section 5 of the paper. 
  * `dengue_original.csv`: the original dengue data from [Tuan *et al.* (2015)](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0003638)
  * `dengue_example.R`: the code for performing the full dengue analysis and simulations (Figure 4 and Figure 5)
  * `dengue_power_vs_positives.csv` and `dengue_power_vs_positives_subset.csv`: results used to create Figure 5 (left and right, respectively)
