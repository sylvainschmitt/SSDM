## Test environments 

* local Ubuntu 22.04 LTS, R 4.3.1
* GitHub Actions Ubuntu 20.04 LTS, R-release

## R CMD check results

There were no ERRORs or WARNINGs.

There was one NOTE:

  Imports includes 24 non-default packages.
   Importing from so many packages makes the package vulnerable to any of
   them becoming unavailable.  Move as many as possible to Suggests and
   use conditionally.
   
The fact that we quickly updated the package following several dependency
changes shows that we are able to assume all the needed dependencies.

## Resubmission

This is a resubmission. In this version I have:

- removed rgdal and sp from the package dependencies following r-spatial package archiving
- updated invalid URLs
- removed lost braces from documentation
