# Diagnosing Positivity Violations in Data

A comparison of the Positivity Regression Tree (PoRT) algorithm developed by Chatton et al., and a kernel-based sparsity diagnostic developed by Ring et al., that aim to diagnose positivity violations in observational data.

The performance is compared on simulated data for different scenarios, i.e. with varying distributional assumptions and number of confounders in the adjustment set, and different forms of positivity violations (see files `port_kbsd_<#confounders>.R`).

It is further compared on data on antiretroviral therapy from the study by Kerschberger et al. (see file `art_data_comparison.R`).

