NEWS
================

SSDM 0.2.6
===============

- CRAN v0.2.6 submission

SSDM 0.2.5.9001
===============

- #36 and #49 fixed

SSDM 0.2.5.9000
===============

- #46 fixed
- Fixed division by zero bug (NA propagation) in `project`

SSDM 0.2.5
===============

- CRAN v0.2.5 submission

SSDM 0.2.3.9011
===============

- Fixed dependency compatibility with raster 2.9-5
- Changed maintainter mail adress (permanent adress)
- Removed rgdal from vignette builder as asked by mail

SSDM 0.2.3.9010
===============

- solved compatibility issue with new raster version

SSDM 0.2.3.9009
===============

- evaluation unactivation with `eval=F` in stacking function for Robin issue

SSDM 0.2.3.9008
===============

- issue #32 fixed related to `shinyFiles` version 0.7.0

SSDM 0.2.3.9007
===============

-   fixed GBM and removed MAXENT for MEM issue \#24 from @BoiMau

SSDM 0.2.3.9006
===============

-   fixed second issue \#24 from @BoiMau

SSDM 0.2.3.9005
===============

-   fixed issue \#24 from @BoiMau

SSDM 0.2.3.9004
===============

-   fixed issue \#20 from @Rekyt

SSDM 0.2.3.9003
===============

-   fixed issue \#22 from @Rekyt

SSDM 0.2.3.9002
===============

-   `rmarkdown` added to VignetteBuilder field in DESCRIPTION
-   CTATION file added following article publication

SSDM 0.2.3.9001
===============

-   axes contribution evaluation when only one variable

SSDM 0.2.3
==========

-   CRAN submission following article submission in MEE

SSDM 0.2.2.9002
===============

-   Occurrences.csv and TRAVIS gdal
-   Vignettes 2, misspelling, and TRAVIS
-   Flo & Dim check

SSDM 0.2.2.9001
===============

-   SSDM and GUI vignettes

SSDM 0.2.9000
=============

-   AppVeyor integration

SSDM 0.1.9040
=============

-   `mapDiversity` S4 methods for SSDM with pSSDM, bSSDM, Bernoulli, MaximumLikelyhood, PRR.MEM, PRR.pSSDM

SSDM 0.1.9039
=============

-   stylistic rules correction with formatR and goodpractice

SSDM 0.1.9038
=============

-   Travis-CI 0.1.9037 fixed
-   Pre formatR test

SSDM 0.1.9037
=============

-   community evaluations for SSDM (*see Pottier et al*) in `evaluate.Stack.SDM`
-   SSDM evaluation in doc, plot and gui

SSDM 0.1.9036
=============

-   doc about new stacking methods (including literature)
-   include new stacking method in GUI

SSDM 0.1.9036
=============

-   further testing of probability ranking stacking method with real data
-   `project.R` with MEM bug fixed
-   `stacking.R` with MEM bug fixed
-   Travis-CI note removed with .Rbuildignore

SSDM 0.1.9035
=============

-   Adjusted binaries with probability ranking method
-   Add binary raster slot in SDM, ESDM and SSDM methods
-   Add binary computation in modelling, ensemble, stacking
-   Add binary tmp save in ensemble and stack modelling
-   Save binary in save function
-   Load binary in load function
-   Create binary if not in file with load function (for backward compatibility)
-   Adapt plot methods to new binary slot

SSDM 0.1.9034
=============

-   `spThin` package responsible for the doc issue, not fixed but currently just commenting the pacakge to unactivate it in waiting. thread open on `Roxygen2` development repository: <https://github.com/klutometis/roxygen/issues/597>.

SSDM 0.1.9033
=============

-   Algortihm.SDM.R file splitted in multiple files (searching for doc issue)

SSDM 0.1.9032
=============

-   richness input in function stacking
-   MEM computing in function stacking
-   backbone to implement all new stacking functions (warning implemented everywhere to alert not implemented parts)
-   update in doc
-   update in check args

SSDM 0.1.9031
=============

-   `rgdal`issue on travis due to `test_load_occ`

SSDM 0.1.9030
=============

-   `rgdal`issue on travis due to `load_var`

SSDM 0.1.9029
=============

-   `rgdal` in DEPENDENCIES for testthat in Travis
-   cran-comment.md
-   NEWS link in README
-   .Rbuildignore

SSDM 0.1.9028
=============

-   Travis test

SSDM 0.1.9027
=============

-   formal testing with `testthat` (39%)
-   `quit` button in gui
-   `stack_modelling`example fix
-   warning messages corrected or removed
-   T/F replace by `TRUE` and `FALSE`
-   `goodpractice` package check
-   `length(x)`replaced by `seq_len(x)`
-   `<-` instead of `=` in examples
-   Package startup message
-   `shinyFiles` in DEPENDENCIES
-   Spelling checked in `plot(SSDM)`
-   Raw `Env`and `Occurrences` data
-   Example data in `gui`
-   Fixed working directory in `gui`
-   Repository and bug report in DESCRIPTION
-   Travis support

SSDM 0.1.9026
=============

-   checkargs
-   stack\_modelling cores
-   NEWS & README

SSDM 0.1.9025
=============

Dimitri Justeau: endemism parameter bug in the gui fixed

SSDM 0.1.9024
=============

-   Enable host and port option when launching the gui [Dimitri Justeau](https://github.com/dimitri-justeau)
-   README

SSDM 0.1.9023
=============

-   tmmpath in stack\_modelling()
-   split with '.Ensemble.SDM' instead of '.' in save, plot, and gui

SSDM 0.1.9022
=============

-   ESDM list in plot.model bug fixed

SSDM 0.1.9021
=============

-   CRAN second submit changes

SSDM 0.1.9020
=============

-   CRAN first submit changes

SSDM 0.1.9019
=============

-   Occurrences man url corrected
-   Check -as--cran
-   Check on CRAN win-build

SSDM 0.1.9018
=============

-   save tab in gui debug for ESDM (Jay)

SSDM 0.1.9017
=============

-   tmppath in ensemble\_modelling dir creation (Jay)
-   tmp doc in ensemble\_ and stack\_modelling [Florian de Boissieu](https://github.com/floriandeboissieu)

SSDM 0.1.9016
=============

-   Weighted endemism index major change : range definition
-   load\_occ major beug fixed, wrong rows were removed
-   tmppath change in ensemble and stack\_modelling (overwrite of path var issue) [Florian de Boissieu](https://github.com/floriandeboissieu)
-   raster::readAll to force loading Env in memory (load\_var) [Florian de Boissieu](https://github.com/floriandeboissieu)
-   Raster\[Raster\[\]&lt;=-900\]=NA instead of reclassify [Florian de Boissieu](https://github.com/floriandeboissieu)
-   Resample with more condition to avoid extra computing [Florian de Boissieu](https://github.com/floriandeboissieu)
-   Add timestampt to ensemble\_modelling [Florian de Boissieu](https://github.com/floriandeboissieu)

SSDM 0.1.9015
=============

-   n.cores parameter in GBM auto-adapted to parallel computing in stack\_modelling
-   gc() added at the end of each modelling function to avoid memory loss
-   ensemble.AlgoSDM() weight corrected
-   .cran-comment in .Rbuildignore
-   shinyFiles in Suggests and checked when gui() is used
-   Check with 0 Notes 0 Warnings 0 Errors

SSDM 0.1.9014
=============

-   Parallel computing with parallel::parLapply
-   .csv in save\_model [Dimitri Justeau](https://github.com/dimitri-justeau)
-   Species in ./Species in save\_stack [Dimitri Justeau](https://github.com/dimitri-justeau)
-   Stack results in ./Stack in save\_stack [Dimitri Justeau](https://github.com/dimitri-justeau)
-   Pearson cor = NA or NULL in evaluate.axis

SSDM 0.1.9013
=============

-   Parallel computing: for loop correctly removed
-   load\_var tested

SSDM 0.1.9012
=============

-   Parallel computing for stack\_modelling with parallel::mclapply

SSDM 0.1.9011
=============

-   load\_var format and reso + extent bugged fixed [Dimitri Justeau](https://github.com/dimitri-justeau)

SSDM 0.1.9010
=============

-   uncertainties stack bug catch in try
-   load\_var with only folder path debug [Dimitri Justeau](https://github.com/dimitri-justeau)
-   load\_var parameter: factors -&gt; categorical [Dimitri Justeau](https://github.com/dimitri-justeau)
-   Presence / absence occurrences modelling (Robin Pouteau)

SSDM 0.1.9009
=============

-Stacking Algo. Corr. row names duplicate (strsplit) - Duplicated '.tif' in load\_occ doc[Dimitri Justeau](https://github.com/dimitri-justeau) - Null supported by Spcol (default) in .checkargs [Dimitri Justeau](https://github.com/dimitri-justeau)

SSDM 0.1.9008
=============

-   Windows gui data input fix (default switch)

SSDM 0.1.9007
=============

-   load\_occ porr data load improvement (occ not in spatial extent + less than 3 occ)

SSDM 0.1.9006
=============

-   ShinyApp structure in inst and shinyFiles working

SSDM 0.1.9005
=============

-   All inputs changed with shinyFiles package + Pcol text + Description (gui and pkg doc)

SSDM 0.1.9004
=============

-   Previous model input change (windows comptibilty)

SSDM 0.1.0
==========

-   Initial private beta release!
