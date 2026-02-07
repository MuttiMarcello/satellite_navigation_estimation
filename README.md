# Satellite Navigation and Estimation

Individual project implementing and assessing uncertainty propagation and navigation filters for LEO satellite and formation-flying scenarios.

## Overview

This repository contains three simulation-based studies developed as part of a university satellite navigation project.

The work is organized in three parts:
1. Comparison of linearized, unscented, and Monte Carlo uncertainty propagation methods for a formation-flying satellite pair
2. Batch filter-based orbit estimation using simulated spacecraft motion and measurements
3. Relative CW motion navigation using an Unscented Kalman Filter

## Results and Validation

Key results:
- Linear Transform effective in capturing uncertainty propagation for circular orbits, fails with non-zero eccentricity. Unscented Transform effective over long integration periods, approximating Monte Carlo validation results.
- Weighted least squares batch filter performance proportional to number of measurements and accuracy of dynamical model. Estimate uncertainty decreases of two order of magnitude adding J2 perturbation to 2BP dynamical model.
- Absolute motion 2BP+J2 UKF 3-sigma uncertainty converges to ~ 0.3 km, ~ 0.5 m/s in 10 minutes, with measurement period of 5s. Relative motion CW UKF 3-sigma uncertainty converges to ~ 0.001 km, ~ 0.001 m/s in 20 minutes, with measurement period of 5s.

Representative outputs:
- Covariance ellipses plots, with Monte Carlo samples
- Batch filter least squares residuals
- UKF uncertainty history

Representative figures are available in 'results/' (PNG format).
See 'results/results.txt' for figure-by-figure explanations.
The full methodology and results are documented in 'docs/report.pdf'.

## Repository structure

- 'src/' - MATLAB implementations of each study
- 'docs/' - Project report
- 'results/' - Key result figures (PNG)
- 'figures/' - Source figures (EPS)

## Development notes

This repository was uploaded after project completion.
Commit history does not reflect the original development timeline.

## Reproducibility and external dependencies

This project relies on external SPICE kernels for ephemerides and constants, not included in this repository.
To run the simulations, the following NAIF kernels are required:
- 'naif0012.tls'
- 'pck00010.tpc'
- 'earth_720101_070426.bpc'
- 'earth_070425_370426_predict.bpc'
- 'earth_fixed.tf'
- 'de425s.bsp'
- 'gm_de431.tpc'
- 'estrack_v04.bsp'
- 'estrack_v04.tf'

The code calls the kernels through meta-kernel named 'assignment02.tm', not included in this repository.

Some scripts call proprietary helper functions provided by the course staff, not included in this repository.

Missing dependencies:
- 'twoline2rv.m' (used in 'batch_filters.m', 'relative_motion_ukf.m' to convert the two line element set character string data to variables and initialize the sgp4 variables)
- 'invjday.m' (used in 'batch_filters.m', 'relative_motion_ukf.m' to finds the year, month, day, hour, minute and second given the julian date)
- 'sgp4.m' (used in 'batch_filters.m', 'relative_motion_ukf.m' to find position, velocity given sgp4 prediction model from space command)
- 'teme2eci.m' (used in 'batch_filters.m', 'relative_motion_ukf.m' to transform a vector from the true equator mean equinox system (teme) to the mean equator mean equinox (j2000) system)
