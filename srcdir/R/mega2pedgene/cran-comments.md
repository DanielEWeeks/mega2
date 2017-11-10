## Test environments
* Ubuntu 16.04.1, R devel (2017-08-17 r73102)
* OS X 10.12.5, R 3.4.1 (2017-06-30)

## R CMD check results
On OS X, there were no ERRORs or WARNINGs.

On Ubuntu, the status is reported as "1 NOTE" about a
possibly invalid URL.

The URL warning is a function of the University of Pittsburgh not
updating their security certificates in a timely manner.  This is
beyond our control, and no such URL warning is generated when tested
from the OS X system.

The required  'mega2r' package was submitted to CRAN today, before this
package was submitted.


## Excerpts from the Ubuntu 00check.log:

* using R Under development (unstable) (2017-08-17 r73102)
* using platform: x86_64-pc-linux-gnu (64-bit)
* this is package ‘mega2pedgene’ version ‘1.0.0’
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Daniel E. Weeks <weeks@pitt.edu>’

New submission

Strong dependencies not in mainstream repositories:
  mega2r

Found the following (possibly) invalid URLs:
  URL: https://watson.hgen.pitt.edu/mega2/mega2r/
    From: DESCRIPTION
    Status: Error
    Message: libcurl error code 60:
    	SSL certificate problem: unable to get local issuer certificate
    	(Status without verification: OK)

Status: 1 NOTE
