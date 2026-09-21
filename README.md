# PLSIM-GOF

Code and simulation results for goodness-of-fit testing in partially linear single-index models.

## Files

- `code/`: core R and Python functions:
  - `SIM_functions.R`: functions for single-index model estimation and simulation.
  - `Orthogonalization.R`: functions for the orthogonalization procedure.
  - `python_code.py`: Python functions for estimating conditional expectations using neural networks.
- `Low_dimension_example123.R`: low-dimensional simulation example.
- `High_dimension_example_HD123.R`: high-dimensional simulation example.
- `Correlated_settings.R`: simulation under correlated covariates.
- `real_data.R`: real-data analysis.
- `results/`: numerical results from the low-dimensional, high-dimensional, and correlated-covariate simulation studies.
- `plot/`: plotting code and figures for the low-dimensional simulation examples.
- `requirements.txt`: Python package requirements.
- `R_session_info.txt`: information on the R environment used for the numerical experiments.

## Usage

The R scripts call Python through the `reticulate` package. Please change the Python interpreter path in the scripts to the Python environment on your local machine if needed.

Numerical experiments were conducted on multiple computing devices. The provided `R_session_info.txt` and `requirements.txt` document one representative R/Python environment used in the experiments, in which Python 3.13.5 (64-bit) was used.

The implementation of the comparison method `Ln` is not included in this repository. The code used to obtain the reported `Ln` results was provided privately by the authors of the original paper upon request and is not redistributed here. The authors' original implementation and parameter settings were used without modification.

## Reproducing Examples 1--3

The baseline single-index link for Examples 1--3 is selected manually in `Low_dimension_example123.R` and `High_dimension_example_HD123.R`. Before running a given example, uncomment the corresponding `Y_base` line and comment out the other two. The `scenario` argument specifies the null or alternative setting. The corresponding configurations are as follows:

| Example | Baseline single-index link | Null scenario | Alternative scenario |
| --- | --- | --- | --- |
| Example 1 / HD-1 | `t_val^2` | `H0` | `case1` |
| Example 2 / HD-2 | `cos(2*t_val)` | `H0` | `case2` |
| Example 3 / HD-3 | `exp(-t_val^2)` | `H0` | `case3` |

The reported low-dimensional simulations use N=800, with p=q=s=10 or 20, whereas the high-dimensional simulations use N=1000, s=20, and p=q=600 or 2000. The master random seed is 123. The corresponding results are stored in `results/different_c/` and `results/high/`, respectively.

## Simulation replications

Each simulation configuration was attempted 1000 times. One replication in the alternative setting of Example HD-2 with \(p=q=600\) and two replications in the high-dimensional correlated-covariate alternative setting did not yield valid test p-values and were recorded as missing values. The corresponding rejection rates were calculated using the successful replications only; all other reported configurations used the full 1000 replications.


