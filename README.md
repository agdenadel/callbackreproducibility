# callback (Calibrated Clustering via Knockoffs) Reproducibility <img src="man/figures/callback_logo.png" align="right" alt="" width="120"/>


## Introduction

We hope to make it as simple as possible to reproduce the results found in `A knockoff calibration method to avoid over-clustering in single-cell RNA-sequencing`.
To that end, we have organized our analysis scripts as a series of `R` vignettes in this repository.

## Reproducing Figures from the `callback` manuscript

Clone the repository:

```bash
git clone https://github.com/lcrawlab/callbackreproducibility
```

Then, navigate to the repo directory and launch R
```bash
cd callbackreproducibility
R
```

You can build the entire website using the following command.
```r
pkgdown::build_site()
```

Note that each `Rmarkdown` file is not fully run by `R` by default. To properly run an have an `Rmarkdown` file run during the website building process, you need to remove 

```R
knitr::opts_chunk$set(eval = FALSE)
```

from the file header. We do not recommend doing this for the vignettes that actually use `callback`, `sc-SHC`, and `CHOIR` for clustering. Rather, the R portions of these files should be put in a script and run using `Rscript`.


## Relevant Citations
A. DenAdel, M. Ramseier, A. Navia, A. Shalek, S. Raghavan, P. Winter, A. Amini, and L. Crawford. A knockoff calibration method to avoid over-clustering in single-cell RNA-sequencing. _bioRxiv_.

## Questions and Feedback
For questions or concerns with `callbackreproducibility` or the `callback` `R` package, please contact
[Alan DenAdel](mailto:alan_denadel@brown.edu) or [Lorin Crawford](lcrawford@microsoft.com). Any feedback on the manuscript or figure reproducibility is appreciated.
