## Contents

There are 4 subdirectories, corresponding to the 4 variants of the model (see paper for details).

For each of the variants:

- `functions_global.py` is a module containing the global parameters of the model, as well as the PDEs and the cost functions.

- `parameter_estimation.py` is the main program for parameter estimation. It can be run as:

    ```
        $ python parameter_estimation.py
    ```

    It outputs the best fit parameters to the `params` folder in the parent directory. Two files are generated: `best_fit_params_base_(model_variant).csv` for the non-adjuvanted vaccine, and `best_fit_params_adj_(model_variant).csv` for the adjuvanted formulations. Note that model `A*` uses the same base parameters as model `A`. 

- `results.ipynb` is a Jupyter notebook that can be used to interactively explore the output of the main program for the best parameter values. In particular, reproduce all model plots from the paper.

---

## Requirements

Packages needed for the main program:

- `numpy` (version `1.15.0` used)

- `scipy` (version `1.1.0` used)

- `pandas` (version `0.20.1` used)

Additional requirements for interactive exploration and figure generation:

- Jupyter Notebook (version `5.0.0` used)

- `matplotlib` (version `2.0.2`)

- LaTeX (for plot labels).

**Note:** All these were written in Python 2.7, but should work in Python 3.x with either minor or no modifications.
