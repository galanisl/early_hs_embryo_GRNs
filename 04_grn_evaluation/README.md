# Scripts for statistical and biological evaluation of the inferred GRNs

The scripts provided in this folder are computationally and time consuming. We recommend running in a high-performance computing environment.

- `eval_cv.R`: Evaluates the reproducibility of the GRNs predicted by different methods using a k-fold cross-validation approach. Note that this script requires access to data provided in the `../02_data_processing/final_data` folder.
- `kfold_cv_evaluation.R`: Functions to facilitate GRN reproducibility assessments via k-fold cross-validation.
