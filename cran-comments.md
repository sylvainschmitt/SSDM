## Test environments
* local Ubuntu 16.04 LTS, R 3.4.3
* Travis CI Ubuntu 14.04.5 LTS, R 3.4.2
* AppVeyor Windows Server 2012, R 3.4.3
* R-hub Fedora Linux, R-devel, clang
* R-hub Ubuntu Linux 16.04 LTS, R-release
* R-hub Windows Server 2008 R2 SP1, R-devel
* Win-builder Windows, R-3.7.0

## R CMD check results
There were no NOTEs, ERRORs or WARNINGs.

## Resubmission
This is a resubmission. In this version I have:

* Fixed dependency compatibility with raster 2.9-5
* Changed maintainter mail adress (permanent adress)
* Removed rgdal from vignette builder as asked by mail
