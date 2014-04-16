## Documentation ##

Please see the Mega2 documentation, which is available in this
distribution as a PDF in the mega2_html folder:

>mega2_html/Mega2_Documentation.pdf

It is also available online as html <http://watson.hgen.pitt.edu/docs/mega2_html/mega2.html> and as a pdf <http://watson.hgen.pitt.edu/docs/mega2_html/Mega2_Documentation.pdf>.

## Mega2 Releases ##

 1. 4.6.3 -- fixes a bug in PLINK format introduced by 4.6.2

 2. 4.6.2 -- adds VCF (Variant Call Format) to Mega2

 3. 4.6.1 -- generates PLINK/SEQ output; adds support for native and MINGW Windows 8

 4. 4.6.0 -- improves memory efficency: 2 bits for markers with two allele

## A note about Mega2 versions available via git: ##
  
  + The master branch reflects our latest working code; **it may contain bugs.**

  + The current release is named v4.6.2\_Release; you should "git checkout v4.6.2\_Release" after cloning Mega2.

  + *There is a bug in v4.6.2\_Release relating to specifying PLINK input files.*  This bug is fixed in the v4.6.2\_Branch; you should "git checkout v4.6.2\_Branch" after cloning to get it.  To avoid confusion, the internal Mega2 version number of this new code is 4.6.3.
