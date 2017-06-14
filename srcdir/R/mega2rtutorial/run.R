run = function(export = "") {

    library(mega2vcf)

    DIR = getSrcDirectory(run)

    read.Mega2DB(paste0(DIR, "/seqsimr.db"), verbose = TRUE)

    NDIR = paste0(DIR, "/svcf")
    if (! dir.exists(NDIR))
        dir.create(NDIR)

    Mega2VCF(paste0(NDIR, "/vcf.01"))

}

run()                
