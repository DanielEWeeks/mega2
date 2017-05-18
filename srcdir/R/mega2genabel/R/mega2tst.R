
#   Mega2: Manipulation Environment for Genetic Analysis
#   Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   Daniel E. Weeks, and University of Pittsburgh
#  
#   This file is part of the Mega2 program, which is free software you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation
#   either version 3 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================

#library(GenABEL)


#' generate required .ped, .fam and .map  for PLINK PED files
#'
#' @description
#'  Use provided gwaa.class-object and create a PLINK PED file, PLINK MAP file
#'  and a PLINK phe file to hold the phenotypes.
#'
#' @param gwaa_ name of gwaa.data-class object
#'
#' @param pfx prefix for PLINK ped file names
#'
#' @return None
#'
#' @export
#' @importFrom utils data
#'
#' @examples
#'\dontrun{
#' dmpPed()
#'}
dmpPed = function(gwaa_ = srdta, pfx = "srdta") {

    dfphe = data.frame(gwaa_@phdata)
    dfphe$sex = dfphe$id
    names(dfphe)[1:2] = c("FID", "IID")
    write.table(dfphe, file=paste0(pfx, ".phe"), sep="\t", quote=FALSE,
                row.names=FALSE, col.names=TRUE)


    dfmap = data.frame(Chromosome=gwaa_@gtdata@chromosome, Name=gwaa_@gtdata@snpnames,
                       Map.k.a=0, BP.p=gwaa_@gtdata@map)
    write.table(dfmap, file=paste0(pfx, ".map"), sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE)

    dfped = data.frame(pid=gwaa_@gtdata@idnames, person=gwaa_@gtdata@idnames )
    dfped$father = 0
    dfped$mother = 0
    dfped$sex    = gwaa_@phdata$sex
    dfped$sex[dfped$sex == 0] = 2
    dfped$default    = gwaa_@phdata$bt


    predfped = sub("/", " ", as.character(gwaa_@gtdata))
    predfped[is.na(predfped)] = "0 0"
    dfped = data.frame(dfped, predfped, stringsAsFactors=FALSE)
    write.table(dfped, file=paste0(pfx, ".ped"), sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
}

#' compare two gwaa.data-class objects
#'
#' @description
#'  verify that the fields in two gwaa.data-class objects
#'
#' @param mega_ name of first gwaa.data-class object
#'
#' @param gwaa_ name of first gwaa.data-class object
#'
#' @return None
#'
#' @export
#' @importFrom utils data
#'
#' @examples
#'\dontrun{
#' tst()
#'}
tst = function (mega_ = mega, gwaa_ = srdta) {
    tstph()
    tstgt()
}

#' test phenotype fields
#'
#' @description
#'  verify that the fields in two gwaa.data-class objects have the same phenotype fields
#'
#' @param mega_ name of first gwaa.data-class object
#'
#' @param gwaa_ name of first gwaa.data-class object
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' tstph()
#'}
tstph = function (mega_ = mega, gwaa_ = srdta) {
    print(all(! is.na(mega_@phdata$sex) && ! is.na(gwaa_@phdata$sex) &&
              mega_@phdata$sex == gwaa_@phdata$sex))

    print(all(! is.na(mega_@phdata$age) && ! is.na(gwaa_@phdata$age) &&
              mega_@phdata$age == gwaa_@phdata$age))

    print(all(! is.na(mega_@phdata$qt1) && ! is.na(gwaa_@phdata$qt1) &&
              mega_@phdata$qt1 == gwaa_@phdata$qt1))

    print(all(! is.na(mega_@phdata$qt2) && ! is.na(gwaa_@phdata$qt2) &&
              mega_@phdata$qt2 == gwaa_@phdata$qt2))

    print(all(! is.na(mega_@phdata$qt3) && ! is.na(gwaa_@phdata$qt3) &&
              mega_@phdata$qt3 == gwaa_@phdata$qt3))

    print(all(! is.na(mega_@phdata$bt) && ! is.na(gwaa_@phdata$bt) &&
              mega_@phdata$bt == gwaa_@phdata$bt))

    print(all(! is.na(mega_@phdata$default) && ! is.na(gwaa_@phdata$bt) &&
              mega_@phdata$default == gwaa_@phdata$bt))
}

#' test genotype fields
#'
#' @description
#'  verify that the fields in two gwaa.data-class objects have the same genotype fields
#'
#' @param mega_ name of first gwaa.data-class object
#'
#' @param gwaa_ name of first gwaa.data-class object
#'
#' @return None
#'
#' @export
#'
#' @examples
#'\dontrun{
#' tstgt()
#'}
tstgt = function (mega_ = mega, gwaa_ = srdta) {
#    print(all(! is.na(mega_@gtdata$bt) && ! is.na(gwaa_@gtdata$bt) &&
#              mega_@gtdata$bt == gwaa_@gtdata$bt))

    print(all(mega_@gtdata@nids == gwaa_@gtdata@nids))
    print(all(mega_@gtdata@nsnps == gwaa_@gtdata@nsnps))
    print(all(mega_@gtdata@nbytes == gwaa_@gtdata@nbytes))

#    print(all(mega_@gtdata@idnames == gwaa_@gtdata@idnames))

    print(all(mega_@gtdata@snpnames == gwaa_@gtdata@snpnames))
    print(all(mega_@gtdata@chromosome == gwaa_@gtdata@chromosome))
#   print(all(mega_@gtdata@coding == gwaa_@gtdata@coding))
#   print(all(mega_@gtdata@strand == gwaa_@gtdata@strand))

    print(all(mega_@gtdata@map == gwaa_@gtdata@map))
    print(all(mega_@gtdata@male == gwaa_@gtdata@male))

    ms = as.character(srdta@gtdata)
    mm = as.character(mega@gtdata)

    ms[ms == "T/G"] = "G/T"
    ms[ms == "T/C"] = "C/T"
    ms[ms == "T/A"] = "A/T"
    ms[ms == "G/C"] = "C/G"
    ms[ms == "G/A"] = "A/G"
    ms[ms == "C/A"] = "A/C"
    ms[is.na(ms)]   = "0/0"

    mm[mm == "T/G"] = "G/T"
    mm[mm == "T/C"] = "C/T"
    mm[mm == "T/A"] = "A/T"
    mm[mm == "G/C"] = "C/G"
    mm[mm == "G/A"] = "A/G"
    mm[mm == "C/A"] = "A/C"
    mm[is.na(mm)]   = "0/0"

    print(all(mm == ms))
}
