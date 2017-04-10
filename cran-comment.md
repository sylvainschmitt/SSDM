## Test environments
* local Ubuntu 14.04 LTS, R 3.2.2
* local Windows 7, R 3.2.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There was one NOTE:
Possibly mis-spelled words in DESCRIPTION:
  SDM (11:73)
  SDMs (10:53, 11:45, 14:56)
  SSDM (10:34, 17:23)
  Thresholding (15:61)
  endemism (9:49)
  thresholding (15:34)
  
Those are scientific correctly spelled words of species distribution modelling field (SDM). It follows litterature included in pacakge documentation.

## Resubmission
This is a resubmission. In this version I have:

* Still an issue with graphical user interface, shiny was not loaded by the UI. Bug is now fixed.
