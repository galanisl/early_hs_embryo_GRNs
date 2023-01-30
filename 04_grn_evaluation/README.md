# Scripts for statistical and biological evaluation of the inferred GRNs

The scripts provided in this folder are computationally demanding and time consuming. We recommend running in a high-performance computing environment. The scripts in the `grid_engine` folder were written specifically for the [UCL Advanced Research Computing platform](https://www.rc.ucl.ac.uk/docs/Experienced_Users/#batch-system).

- `eval_cv.R`: Evaluates the reproducibility of the GRNs predicted by different methods using a k-fold cross-validation approach. Note that this script requires access to data provided in the `../02_data_processing/final_data` folder.
- `kfold_cv_evaluation.R`: Functions to facilitate GRN reproducibility assessments via k-fold cross-validation.
- `TFnetSimGen.R`: Generation of synthetic networks for the study TF and sample size impact on GRN inferences.
- `grid_engine`:
    - `TFnetRviolsAll.R`: R and V stat calculations for the different GRNs inferred.
    - `TFnetSimGENIE3HPC.R`: Application of GENIE3 to simulated GRNs.
    - `TFnetSimL0L2HPC.R`: Application of L0L2 to simulated GRNs
    - `TFnetSimMIHPC.R`: Application of MI to simulated GRNs
    - `TFnetSimSpearman.R`: Application of Spearman correlation to simulated GRNs
