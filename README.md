# OzoneExposure: An R package for modeling the impact of human exposure to ozone.

This package implements two models for predicting the change in FEV1
over time after exposure to ozone. Both models are written in
[Stan](mc-stan.org) for speed, and are fit to the same dataset.

This package started with both models implemented in R
[here](https://github.com/dsidavis/Ozone), and this code has been
preserved here for reference. However, the Stan versions have a couple
of notable advantages:

1. Both models are significantly faster
1. The models have been written to use the exact same input data,
   allowing users to subset the original dataset and compare each
   model.
1. The models calculate the AIC, using the same calculation and
   functions for the log-likelihood.

## Installation

The easiest method is to install directly from github using the
devtools package by running:

```
devtools::install_github("dsidavis/Ozone2")
```

