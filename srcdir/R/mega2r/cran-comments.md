## Thank you

Thank you to Uwe Ligges for helpful feedback on the
previous submission.  As recommended, we have combined
everything into a single R package.

We have also added single quote software names as recommended.

And we have added a web reference to Mega2 in the Description
field in the form <http...>.

Thank you.


## Test environments
* Ubuntu 16.04.1, R devel (2017-08-22 r73111)
* OS X 10.12.5, R 3.4.1 (2017-06-30)

## R CMD check results
On OS X, there were no ERRORs or WARNINGs.

On Ubuntu, the status is reported as "1 NOTE" about a
possibly invalid URL.

The URL warning is a function of the University of Pittsburgh not
updating their security certificates in a timely manner.  This is
beyond our control, and no such URL warning is generated when tested
from the OS X system.


## Excerpts from the 00check.log generated on the Ubuntu system:

* using R Under development (unstable) (2017-08-22 r73111)
* using platform: x86_64-pc-linux-gnu (64-bit)
* this is package ‘mega2r’ version ‘1.0.0’
* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Daniel E. Weeks <weeks@pitt.edu>’

New submission

Found the following (possibly) invalid URLs:
  URL: https://watson.hgen.pitt.edu/mega2/mega2r/
    From: DESCRIPTION
    Status: Error
    Message: libcurl error code 60:
    	SSL certificate problem: unable to get local issuer certificate
    	(Status without verification: OK)

Status: 1 NOTE
