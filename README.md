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
