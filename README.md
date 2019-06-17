SSDM: Stacked species distribution modelling
================

[![Travis-CI Build Status](https://travis-ci.org/sylvainschmitt/SSDM.svg?branch=master)](https://travis-ci.org/sylvainschmitt/SSDM)[![CRAN](https://www.r-pkg.org/badges/version/SSDM)](https://CRAN.R-project.org/package=SSDM) [![Downloads](http://cranlogs.r-pkg.org/badges/SSDM?color=brightgreen)](http://www.r-pkg.org/pkg/SSDM) [![Coverage Status](https://img.shields.io/codecov/c/github/sylvainschmitt/SSDM/master.svg)](https://codecov.io/github/sylvainschmitt/SSDM?branch=master) [![Research software impact](http://depsy.org/api/package/cran/SSDM/badge.svg)](http://depsy.org/package/r/SSDM)

SSDM is a package to map species richness and endemism based on stacked species distribution models (SSDM). Individual SDMs can be created using a single or multiple algorithms (ensemble SDMs). For each species, an SDM can yield a habitat suitability map, a binary map, a between-algorithm variance map, and can assess variable importance, algorithm accuracy, and between-algorithm correlation. Methods to stack individual SDMs include summing individual probabilities and thresholding then summing. Thresholding can be based on a specific evaluation metric or by drawing repeatedly from a Bernouilli distribution. The SSDM package also provides a user-friendly interface `gui`.

For a full list of changes see [`NEWS`](./NEWS.md).

Installation
============

Please be aware that SSDM package use a lot of dependencies (see [`DESCRIPTION`](./DESCRIPTION))

### Install from Github

You can install the latest version of **SSDM** from Github using the [`devtools`](https://github.com/hadley/devtools) package:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github("sylvainschmitt/SSDM")
```

### Install from CRAN

The stable version of **SSDM**, is available on CRAN:

``` r
install.packages("SSDM")
```

*We advise users to install from github. Due to CRAN policies and the development of SSDM, many new features and bugfixes may be available on CRAN later.*

Usage
=====

After installing, **SSDM** package, you can launch the graphical user interface by typing gui() in the console.

<center>
[**Click to enlarge**](https://raw.githubusercontent.com/sylvainschmitt/SSDM/master/examples/SSDM.gif)
<center>
![Screenshot](https://raw.githubusercontent.com/sylvainschmitt/SSDM/master/examples/SSDM.gif)

Functionnalities
================

SSDM provides five categories of functions (that you can find in details below): Data preparation, Modelling main functions, Model main methods, Model classes, and Miscellaneous.

### Data preparation

-   `load_occ`: Load occurrence data
-   `load_var`: Load environmental variables

### Modelling main functions

-   `modelling`: Build an SDM using a single algorithm
-   `ensemble_modelling`: Build an SDM that assembles multiple algorithms
-   `stack_modelling`: Build an SSDMs that assembles multiple algorithms and species

### Model main methods

-   `ensemble,Algorithm.SDM-method`: Build an ensemble SDM
-   `stacking,Ensemble.SDM-method`: Build an SSDM
-   `update,Stacked.SDM-method`: Update a previous SSDM with new occurrence data

### Model classes

-   `Algorithm.SDM`: S4 class to represent SDMs
-   `Ensemble.SDM`: S4 class to represent ensemble SDMs
-   `Stacked.SDM`: S4 class to represent SSDMs

### Miscellanous

-   `gui`: user-friendly interface for SSDM package
-   `plot.model`: Plot SDMs
-   `save.model`: Save SDMs
-   `load.model`: Load SDMs
