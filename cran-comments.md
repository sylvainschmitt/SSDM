## Test environments 

* local Ubuntu 22.04 LTS, R 4.3.1
* GitHub Actions Ubuntu 20.04 LTS, R-release

## R CMD check results

There were no NOTEs, ERRORs or WARNINGs.

## Resubmission

This is a resubmission. In this version I have:

- removed superseded package snow
- removed rgdal and sp from the package dependencies following r-spatial package archiving
- updated invalid URLs
- removed lost braces from documentation
