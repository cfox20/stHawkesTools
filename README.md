# stHawkesTools

stHawkesTools provides utilities for simulating, fitting, and diagnosing general spatio-temporal Hawkes processes. The package wraps estimation routines based on the EM algorithm together with simulation helpers, visualization utilities, and bootstrap-based uncertainty quantification to support end-to-end analysis workflows.

## Getting started

The development version of `stHawkesTools` is hosted on GitHub.  Install it
with `remotes::install_github()` (or the equivalent `devtools::install_github()`).

After installing, use the [Getting Started vignette](vignettes/getting-started.Rmd) to walk through an example analysis. You can also open the vignette from an R session once the package is installed by running:
```r
vignette("getting-started", package = "stHawkesTools")
```
