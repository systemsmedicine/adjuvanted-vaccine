# Mathematical model of adjuvanted vaccine dynamics

This repository contains different variants of a mathematical model of antibody production dynamics in response to influenza vaccination with adjuvanted and non-adjuvanted vaccine formulations. This code accompanies the paper [Adjuvanted influenza vaccine dynamics](https://www.nature.com/articles/s41598-018-36426-9), by C. Parra-Rojas, V. von Messling, and E. A. Hernández-Vargas.

Vaccination data correspond to [Schmidt *et al.*, Vaccine **34(44)**, 5329 (2016)](http://www.sciencedirect.com/science/article/pii/S0264410X16307897).

---

## Brief model description

The system considers two variables: `B(x,t)` and `Ab(x,t)`, representing the distributions of B cells and antibodies in a one-dimensional shape space, where the NA and HA proteins of the virus are located at positions `x = 0` and `x = 1`, respectively.

Euclidean distance to the antigen determines affinity which, in turn, influences the rates of proliferation and decay of B cells, as well as their ability to produce antibodies. The effects of the adjuvants in the system are introduced as boosts in glycoprotein immunogenicities and, in one variant of the model, the overall rate of antibody production independent of protein specificity.

This simple model is able to capture the essential feature observed in the data, correponding to a shift in protein-specificity and a large enhancement in antibody response for the vaccine formulated with squalene-containing adjuvants (MF59 and AS03), with respect to the non-adjuvanted vaccine and a vaccine formulated with an adjuvant based on vitamin E (Diluvac).

---


## Repository structure

- `code` contains all the files needed to reproduce the results from the paper for the different variants of the model. See the `README.md` file inside this folder for further details.
- `params` contains the best fit parameters obtained by the estimation procedure.
- `figures` contains the outputs from the model with the best fit parameters.
- `data` contains the raw data employed for the fitting procedure.
