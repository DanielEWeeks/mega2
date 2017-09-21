# This file contains pieces of code from GenABEL
# https://CRAN.R-project.org/package=GenABEL 
# GenABEL is GPL-licensed.
#
# ~/rvb/Work/R/pkg/GenABEL/R/alleleID.R lines 1:39
alleleID.alleles <- function() {
	a <- list();
	a[[1]] <- c("1","2")
	a[[2]] <- c("A","B")
	alleles <- c("A","T","G","C","-")
	idx <- 3
	for (i in alleles) {
		for (j in alleles) {
			if (i==j) next;
			a[[idx]] <- c(i,j)
			idx <- idx + 1
		}
	}
	a[[idx]] <- c("2","1")
	idx <- idx + 1
	a[[idx]] <- c("B","A")
	idx <- idx + 1
	a[[idx]] <- c("I","D")
	idx <- idx + 1
	a[[idx]] <- c("D","I")
	idx <- idx + 1
	allalleles <- c("1","2","B","I","D","A","T","G","C","-")
	for (jj in allalleles) {
		a[[idx]] <- c(jj,jj)
		idx <- idx + 1
	}
	a
}

alleleID.codes <- function() {
	a <- alleleID.alleles()
	out <- c("OPPA")
	idx <- 1
	for (i in a) {
		out[idx] <- paste(i[1],i[2],sep="")
		idx <- idx + 1
	}
	out
}

# ~/rvb/Work/R/pkg/GenABEL/R/Xcheck.R
#' @importFrom GenABEL crnames perid.summary
"Xcheck" <-
function(data,Pgte=0.01,Pssw=0.01,Pmsw=0.01,odds=1000,tabonly=FALSE,Fmale=0.8,Ffemale=0.2) {
	if (!is(data,"snp.data")) stop("data argument should be of snp.data-class")
	if (any(data@chromosome != "X")) stop("All markers should be X-linked")
	male <- (data@male==1)
	out <- list()
	out$xerr <- 0
	q <- summary(data)[,"Q.2"]
	if (sum(male)) {
		xdat <- as.numeric(data[male,])
		if (any(xdat==1,na.rm=T)) {
			out$xerr <- 1
			out$Xerrtab <- crnames(dimnames(xdat),which(xdat==1))
			colnames(out$Xerrtab) <- c("ID","SNP")
		}
	}
	if (tabonly) return(out)
	xdat <- as.numeric(data)
	xdat <- t(xdat)
	ll.female <- log(q*q*(1-Pgte)*(xdat==2)+(2*q*(1-q)*(1-Pgte)+Pgte)*(xdat==1)+(1-q)*(1-q)*(1-Pgte)*(xdat==0))
	ll.male <- log(q*(1-Pgte)*(xdat==2)+Pgte*(xdat==1)+(1-q)*(1-Pgte)*(xdat==0))
	ll.sex <- ll.female
	ll.sex[,male] <- ll.male[,male]
# find SNPs with are likely to be (pseudo)autosomal
	snpprob.0 <- apply(ll.sex,MARGIN=1,FUN=sum,na.rm=T)+log(1-Pmsw)
	snpprob.1 <- apply(ll.female,MARGIN=1,FUN=sum,na.rm=T)+log(Pmsw)
	snpODDs <- snpprob.1-snpprob.0
	out$Xmrkfail <- names(snpprob.1[snpODDs>log(odds)])
# find male which are likely to be female
	idprob.0 <- apply(ll.sex,MARGIN=2,FUN=sum,na.rm=T)+log(1-Pssw)
	idprob.1 <- apply(ll.female,MARGIN=2,FUN=sum,na.rm=T)+log(Pssw)
	idODDs <- idprob.1-idprob.0
	out$isfemale <- names(idprob.1[idODDs>log(odds)])
# find female which are likely to be male
	idprob.1 <- apply(ll.male,MARGIN=2,FUN=sum,na.rm=T)+log(Pssw)
	idODDs <- idprob.1-idprob.0
	out$ismale <- names(idprob.1[idODDs>log(odds)])
# find people with strange F
	pis <- perid.summary(data)
	out$otherSexErr <- rownames(pis)[pis$F > Ffemale & pis$F < Fmale]
# return object
	out$Xidfail <- unique(c(out$ismale,out$isfemale,out$othersexErr))
	out
}


# 
#load.gwaa.data =
gwaaCheckPhe = function(dta, id) {
    coln <- names(dta)
    idColumn <- match(id, coln)
    names(dta)[idColumn] <- "id"
    if (!any(names(dta) == "id", na.rm = TRUE)) 
        stop("the filed named \"id\", containing the identifier presented in both pheno- and geno- files was not found in the phenofile")
    class(dta$id) <- "character"
    if (!any(names(dta) == "sex", na.rm = TRUE)) 
        stop("the column named \"sex\", containing the male identifier was not found in the phenofile")
    v <- version
    if (as.numeric(v$major) > 2 || ((as.numeric(v$major) == 2) && 
        (as.numeric(v$minor) >= 8))) {
        a <- table(dta$sex, useNA = "ifany")
    }
    else {
        a <- table(dta$sex, exclude = NULL)
    }
    if (length(a) > 2) 
        stop("column named \"sex\" contains more than 2 codes")
    if (length(a) == 1 && !(names(a)[1] == 0 || names(a)[1] == 
        1)) 
        stop("the column named \"sex\" contains 1 code which is neither 0 (=female) or 1 (=male)")
    if (length(a) == 2 && names(a)[1] != 0 && names(a)[2] != 
        1) 
        stop("the column named \"sex\" is not coded as 0=female and 1=male")
    rm(a)
    gc(verbose = FALSE)
    if (any(table(dta$id) > 1)) 
        stop("there are duplicated IDs in the phenotypic data file")
    if (any(is.na(dta$id))) 
        stop("there are missing IDs in the phenotypic data file")
    if (any(is.na(dta$sex))) 
        stop("there are missing sex values in the phenotypic data file")
    rownames(dta) <- dta$id
    dta
}

gwaaCheckPersons = function(dta, ids) {
    mlst0 <- match(as.character(dta$id), ids)
    for (i in 1:length(mlst0)) {
        cid <- mlst0[i]
        if (is.na(cid)) 
            cat("person with id =", as.character(dta$id)[i], 
                "was not found in genotypic file; excluded\n")
    }
    mlst <- match(ids, as.character(dta$id))
    oerr <- 0
    for (i in 1:length(mlst)) {
        cid <- mlst[i]
        if (is.na(cid)) {
            cat("person with id =", ids[i], "was not found in phenotypic file!!! - FATAL\n")
            oerr <- oerr + 1
        }
    }
    if (oerr) 
        stop("fatal error. update pheno-file")
    newdta <- data.frame(dta[mlst, ])
    newdta
}

#' @importFrom GenABEL sortmap.internal
gwaaEpilog = function(a, dta, force, makemap, sort) {
    if (any(a@chromosome == "X") && any(a@male == 1) && !force) {
        xmrk <- (a@chromosome == "X")
        mlst <- (a@male == 1)
        rxm <- a[mlst, xmrk]
        Xch <- Xcheck(rxm)
        rm(rxm, xmrk, mlst)
        gc(verbose = FALSE)
        if (Xch$xerr) {
            cat("Wrong male X genotypes (heterozygous) found in", 
                dim(Xch$tab)[1], "occasions\n")
            cat("Error table is saved as the output object\n")
            return(Xch$tab)
        }
    }
    if (force) 
        cat("assignment of gwaa.data object FORCED; X-errors were not checked!\n")
    if (makemap) {
        cat("increase in map order FORCED\n")
        chun <- levels(a@chromosome)
        if (any(chun != "X")) {
            numchun <- sort(as.numeric(chun[chun != "X"]))
            gsize <- max(a@map[a@chromosome == as.character(numchun[1])])/5
            if (length(numchun) > 1) {
                for (i in c(2:(length(numchun)))) {
                  inc <- max(a@map[a@chromosome == as.character(numchun[i - 
                    1])]) + gsize
                  a@map[a@chromosome == as.character(numchun[i])] <- a@map[a@chromosome == 
                    as.character(numchun[i])] + inc
                }
            }
            if (any(chun == "X")) {
                inc <- max(a@map[a@chromosome == as.character(numchun[length(numchun)])]) + 
                  gsize
                a@map[a@chromosome == "X"] <- a@map[a@chromosome == 
                  "X"] + inc
            }
        }
    }
    out <- new("gwaa.data", phdata = dta, gtdata = a)
    rm(a, dta)
    gc(verbose = FALSE)
    if (sort) {
        ord <- sortmap.internal(out@gtdata@chromosome, out@gtdata@map)
        out <- out[, ord$ix]
    }
    out
}

#' @importFrom GenABEL snp.data
#' @importFrom utils read.table
gwaaO = function (phenofile = "pheno.dat", genofile = "geno.raw", force = TRUE, 
                  makemap = FALSE, sort = TRUE, id = "id", envir = ENV) 
{

    dta <- read.table(phenofile, header = TRUE, as.is = TRUE)
    dta = gwaaCheckPhe(dta, id)

    ifile <- file(genofile, "r")
    header <- scan(file = ifile, what = character(), nlines = 1, 
        quiet = TRUE)
    vver <- grep(x = header, pattern = "version")
    if (length(vver) > 0) {
        ver <- as.numeric(header[vver + 1])
    }
    else {
        ver <- 0
    }
    if (is.na(ver)) 
        warning("Incorrect data format version number")
    if (ver > 0) {
        ids <- scan(file = ifile, what = character(), nlines = 1, 
            quiet = TRUE)
    }
    else {
        ids <- header
    }
    nids <- length(ids)
    cat("ids loaded...\n")
    mnams <- scan(file = ifile, what = character(), nlines = 1, 
        quiet = TRUE)
    cat("marker names loaded...\n")
    chrom <- scan(file = ifile, what = character(), nlines = 1, 
        quiet = TRUE)
    chrom <- as.factor(chrom)
    gc(verbose = FALSE)
    cat("chromosome data loaded...\n")
    pos <- scan(file = ifile, what = double(), nlines = 1, quiet = TRUE)
    cat("map data loaded...\n")
    if (ver == 0) {
        coding <- new("snp.coding", as.raw(rep(1, length(pos))))
        strand <- new("snp.strand", as.raw(rep(0, length(pos))))
    }
    else {
        coding <- scan(file = ifile, what = raw(), nlines = 1, 
            quiet = TRUE)
        class(coding) <- "snp.coding"
        cat("allele coding data loaded...\n")
        strand <- scan(file = ifile, what = raw(), nlines = 1, 
            quiet = TRUE)
        class(strand) <- "snp.strand"
        cat("strand data loaded...\n")
    }
    nsnps <- length(mnams)
    nbytes <- ceiling(nids/4)
    rdta <- scan(file = ifile, what = raw(), quiet = TRUE)
    cat("genotype data loaded...\n")
    close(ifile)
    dim(rdta) <- c(nbytes, nsnps)
    rdta <- new("snp.mx", rdta)

    gc(verbose = FALSE)
    dta = gwaaCheckPersons(dta, ids)

    gc(verbose = FALSE)
    a <- snp.data(nids = nids, rawdata = rdta, idnames = ids, 
        snpnames = mnams, chromosome = chrom, map = pos, coding = coding, 
        strand = strand, male = dta$sex)
    cat("snp.data object created...\n")
    rm(rdta, ids, mnams, chrom, pos, coding, strand)
    gc(verbose = FALSE)

    gwaaEpilog(a, dta, force, makemap, sort)
}
