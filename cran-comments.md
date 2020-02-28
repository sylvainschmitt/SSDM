## Test environments 

* local Ubuntu 16.04 LTS, R 3.6.0
* Travis CI Ubuntu 14.04.5 LTS, R-release
* R-hub Fedora Linux, R-devel, clang
* R-hub Ubuntu Linux 16.04 LTS, R-release
* R-hub Windows Server 2008 R2 SP1, R-devel
* Win-builder Windows, R-4.0.0

## R CMD check results

There were no NOTEs, ERRORs or WARNINGs.

## Resubmission

This is a resubmission. In this version I have:

# * Added stringsAsFactors = TRUE in read.csv2 from load_occ.R folowwing the new defaulot in R 4.0.0
