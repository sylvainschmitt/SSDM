## Test environments
* local Ubuntu 16.04 LTS, R 3.4.3
* Travis CI Ubuntu 14.04.5 LTS, R 3.4.2
* AppVeyor Windows Server 2012, R 3.4.3
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

* Added `rmarkdown` to VignetteBuilder field in DESCRIPTION following Prof Brian Ripley mail
* Added CTATION file following article publication
