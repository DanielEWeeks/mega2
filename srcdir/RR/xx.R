go = function(db = 0) {
#    if (db) browser()

    library(mega2)

    if (TRUE) {
##  ENV = dbmega2_import("/Users/rbaron/mega2/test/samoan_GWAS/dbmega2.db", verbose = 1)
        source("../R/mega2/tests/mega2_test.R")
        tst1()

        source("../R/mega2/tests/mega2gene_test.R")
        tst10()
        tst11()
        tst12()
    }

    library(mega2pedgene)

    aa=setwd("~/mega2/test/R/yj/realDataAnalysis/mega2")
    init_pedgene(verbose = 1)
    run(gs=1:10)
    aa=setwd(aa)
}

go(1)

