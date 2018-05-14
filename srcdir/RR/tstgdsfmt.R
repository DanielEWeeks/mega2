#   Mega2R: Mega2 for R.
#
#   Copyright 2018, University of Pittsburgh. All Rights Reserved.
#
#   Contributors to Mega2R: Robert V. Baron and Daniel E. Weeks.
#
#   This file is part of the Mega2R program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later
#   version.
#
#   Mega2R is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
#
# ===========================================================================

## For installation do:
## source("http://bioconductor.org/biocLite.R")
## biocLite("gdsfmt")


library(Mega2R)
library(gdsfmt)


Mega2gdsfmtseqtst = function(aa, bb, cc="aa bb") {

    print(cc)
    gdsn = function(aa, bb, name, fn=NULL) {
        tryCatch({
        a1 = read.gdsn(index.gdsn(aa, name))
        a2 = read.gdsn(index.gdsn(bb, name))
        cat("\n")
        print(index.gdsn(aa, name))
        print(index.gdsn(bb, name))
        if (! is.null(fn)) {
            a3 = fn(a1, a2)
        } else {
            a3 = a1==a2
        }
        cat(name, all(a3), "\n")
        print(head(a1))
        print(head(a2))
        }, error = function(e) { cat(name, " ", conditionMessage(e), "\n")} )
    }
    gdsn(aa, bb, "sample.id")
    gdsn(aa, bb, "variant.id")
    gdsn(aa, bb, "position")
    gdsn(aa, bb, "chromosome", function(a,b) {a[a=="X"] = "23"; a==b })

    fn = function(a,b) {
        a=sub(",", "/", a)
        b[b=="G/<CN0>"] = "<CN0>/G";   b[b=="A/<CN0>"] = "<CN0>/A"
        (a==b)
    }
    gdsn(aa, bb, "allele", fn)

    gdsn(aa, bb, "genotype/data")
    gdsn(aa, bb, "genotype/extra.index")
    gdsn(aa, bb, "genotype/extra")

    gdsn(aa, bb, "phase/data", function(a,b) {b=rep(1,length(b));a==b})
    gdsn(aa, bb, "phase/extra.index")
    gdsn(aa, bb, "phase/extra")
    
    gdsn(aa, bb, "annotation/id", function(a,b) {b[grep("chr", b)]="";a==b})
    gdsn(aa, bb, "annotation/qual")
    gdsn(aa, bb, "annotation/filter", function(a,b) {as.character(a)==as.character(b)})
    gdsn(aa, bb, "annotation/info")
    gdsn(aa, bb, "annotation/format")

    gdsn(aa, bb, "sample.annotation/family")
}

Mega2gdsfmtsnptst = function(aa, bb, ss="aa bb") {

    print(ss)
    gdsn = function(aa, bb, name, fn=NULL) {
        tryCatch({
        a1 = read.gdsn(index.gdsn(aa, name))
        a2 = read.gdsn(index.gdsn(bb, name))
        cat("\n")
        print(index.gdsn(aa, name))
        print(index.gdsn(bb, name))
        if (! is.null(fn)) {
            a3 = fn(a1, a2)
        } else {
            a3 = a1==a2
        }
        cat(name, all(a3), "\n")
        print(head(a1))
        print(head(a2))
        }, error = function(e) { cat(name, " ", conditionMessage(e), "\n")} )
    }
    gdsn(aa, bb, "sample.id")
    gdsn(aa, bb, "snp.id")
    gdsn(aa, bb, "snp.rs.id", function(a,b) {b[grep("chr", b)]="";a==b})
    gdsn(aa, bb, "snp.position")
    gdsn(aa, bb, "snp.chromosome", function(a,b) {a[a=="X"] = "23"; a==b })

    fn = function(a,b) {
        a[a=="T/A"] = "A/T"; a[a=="T/C"] = "C/T"; a[a=="T/G"] = "G/T"
        a[a=="G/C"] = "C/G"; a[a=="G/A"] = "A/G"; a[a=="C/A"] = "A/C"
        b[b=="T/A"] = "A/T"; b[b=="T/C"] = "C/T"; b[b=="T/G"] = "G/T"
        b[b=="G/C"] = "C/G"; b[b=="G/A"] = "A/G"; b[b=="C/A"] = "A/C"
        (a==b)
    }
##  gdsn(aa, bb, "snp.allele", fn)
    gdsn(aa, bb, "snp.allele")

    gdsn(aa, bb, "genotype")
}

t2 = function() {
    ENV=read.Mega2DB("~/mega2/data/momenly.db")
    setfam(uniqueFamMember())

    closefn.gds(gn)
    unlink("test.gds")
    gn = Mega2gdsfmt(filename = "test.gds")
      
    if (0) {  ## use AGCT for alleles
        closefn.gds(gn)
        gn=openfn.gds("test.gds", readonly=F)
        a=read.gdsn(index.gdsn(gn, "snp.allele"))
        a[a[] == "A/B"]="A/T"
        a[a[] == "B/A"]="A/T"
        delete.gdsn(index.gdsn(gn, "snp.allele"))
        add.gdsn(gn, "snp.allele", a)
    }
  
    snpset = snpgdsLDpruning(gn, ld.threshold=0.2)
    snpset.id = unlist(snpset)
    pca = snpgdsPCA(gn, snp.id=snpset.id, num.thread=2)

    sample.id = read.gdsn(index.gdsn(gn, "sample.id"))
    cc_code   = read.gdsn(index.gdsn(gn, "sample.annot/cc.id"))
    pop_code = c("CEU", "HCB", "JPT", "YRI")[cc_code + 1]
    YRI.id    = sample.id[pop_code == "CEU"]

    ibd       = snpgdsIBDMoM(gn, sample.id=YRI.id, snp.id=snpset.id,
                             maf=0.05, missing.rate=0.40, num.thread=1)


    set.seed(100)
    snp.id = sample(snpset.id, 1500)
    ibd    = snpgdsIBDMLE(gn, sample.id=YRI.id, snp.id=snpset.id,
                          maf=0.05, missing.rate=0.40, num.thread=1)
    ibd.coeff = snpgdsIBDSelection(ibd)

    family.id  = read.gdsn(index.gdsn(gn, "sample.annot/family.id"))
    family.id  = family.id[match(YRI.id, sample.id)]
    table(family.id)
    ibd.robust = snpgdsIBDKING(gn, sample.id=YRI.id, family.id=family.id, num.thread=1)
    
    dat = snpgdsIBDSelection(ibd.robust)
    head(dat)

}

tx = function() {
  seq=openfn.gds("Gdsfmt/aa.seq.gds")
  snp=openfn.gds("Gdsfmt/aa.snp.gds")
  on.exit(closefn.gds(seq))
  on.exit(closefn.gds(snp))
  ENV=read.Mega2DB("Gdsfmt/dbmega2.db", bp=0)

  sqma=Mega2gdsfmt("sqma.gds", SeqArray=T)
  sqna=Mega2gdsfmt("sqna.gds", SeqArray=T, snp.order=T)
  ssma=Mega2gdsfmt("ssma.gds", SeqArray=F)
  ssna=Mega2gdsfmt("ssna.gds", SeqArray=F, snp.order=T)
  on.exit(closefn.gds(sqma))
  on.exit(closefn.gds(sqna))
  on.exit(closefn.gds(ssma))
  on.exit(closefn.gds(ssna))
  Mega2gdsfmtseqtst(seq, sqma, "native Seq vs mega2 Seq Sam")
  Mega2gdsfmtseqtst(seq, sqna, "native Seq vs mega2 Seq Snp")
  Mega2gdsfmtsnptst(snp, ssma, "native Snp vs mega2 Snp Sam")
  Mega2gdsfmtsnptst(snp, ssna, "native Snp vs mega2 Snp Snp")

  append_genotype_a = append_genotype_b = append_genotype_c = F
  append_genotype_b = T

  sqmb=Mega2gdsfmt("sqmb.gds", SeqArray=T)
  sqnb=Mega2gdsfmt("sqnb.gds", SeqArray=T, snp.order=T)
  ssmb=Mega2gdsfmt("ssmb.gds", SeqArray=F)
  ssnb=Mega2gdsfmt("ssnb.gds", SeqArray=F, snp.order=T)
  on.exit(closefn.gds(sqmb))
  on.exit(closefn.gds(sqnb))
  on.exit(closefn.gds(ssmb))
  on.exit(closefn.gds(ssnb))
  Mega2gdsfmtseqtst(seq, sqmb, "native Seq vs mega2 Seq Sam")
  Mega2gdsfmtseqtst(seq, sqnb, "native Seq vs mega2 Seq Snp")
  Mega2gdsfmtsnptst(snp, ssmb, "native Snp vs mega2 Snp Sam")
  Mega2gdsfmtsnptst(snp, ssnb, "native Snp vs mega2 Snp Snp")

  append_genotype_a = append_genotype_b = append_genotype_c = F
  append_genotype_c = T

  sqmc=Mega2gdsfmt("sqmc.gds", SeqArray=T)
  sqnc=Mega2gdsfmt("sqnc.gds", SeqArray=T, snp.order=T)
  ssmc=Mega2gdsfmt("ssmc.gds", SeqArray=)
  ssnc=Mega2gdsfmt("ssnc.gds", SeqArray=F, snp.order=T)
  on.exit(closefn.gds(sqmc))
  on.exit(closefn.gds(sqnc))
  on.exit(closefn.gds(ssmc))
  on.exit(closefn.gds(ssnc))
  Mega2gdsfmtseqtst(seq, sqmc, "native Seq vs mega2 Seq Sam")
  Mega2gdsfmtseqtst(seq, sqnc, "native Seq vs mega2 Seq Snp")
  Mega2gdsfmtsnptst(snp, ssmc, "native Snp vs mega2 Snp Sam")
  Mega2gdsfmtsnptst(snp, ssnc, "native Snp vs mega2 Snp Snp")

  append_genotype_a = append_genotype_b = append_genotype_c = F
  append_genotype_a = T

  showfile.gds(closeall=T, verbose=T)

}

tx()

function() {
    snp=openfn.gds("Gdsfmt/aa.snp.gds")
    ssma=Mega2gdsfmt("ssma.gds", SeqArray=F)
    ENV=read.Mega2DB("Gdsfmt/dbmega2.db", bp=0)
    Mega2gdsfmtsnptst(snp, ssma, "native Snp vs mega2 Snp Sam")

    o = read.gdsn(index.gdsn(snp, "genotype"))
    n = read.gdsn(index.gdsn(ssma, "genotype"))

}
