## Test environments
* local Ubuntu 14.04 LTS, R 3.2.2
* local Windows 7, R 3.2.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Resubmission
This is a resubmission. In this version I have:

* Remove the term package in the first sentence of Description field in DESCRIPTION

* Actualized Date field in DESCRIPTION to today

* Import following functions in NAMESPACE : aggregate.data.frame arrows barplot cor glm glm.control heat.colors is.raster legend lsf.str rainbow rbinom read.csv read.csv2 runif sd tail terrain.colors var write.csv

* Unactivated examples with an elapsed time > 5s
