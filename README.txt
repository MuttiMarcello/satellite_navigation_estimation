# Satellite Navigation and Estimation

Academic project implementing and assessing uncertainty propagation and navigation filters.

## Overview

This repository contains three simulation-based studies developed as part of a university satellite navigation project.

The work is organized in three parts:
1. Comparison of linearized, unscented, and Monte Carlo uncertainty propagation methods
2. Batch filter-based orbit estimation using simulated spacecraft motion and measurements
3. Relative motion navigation using an Unscented Kalman Filter

The full methodology and results (including plots) are available in 'docs/report.pdf'.

## Repository structure

- 'src/' - MATLAB implementations of each study
- 'docs/' - Project report

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
