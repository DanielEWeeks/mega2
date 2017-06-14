library(mega2genabel)

go = function() {

    print("## generate genabel for srdta")

    ENV = read.Mega2DB("srdta.db", verbose = 1)

    mega = Mega2GenABEL("xyz", envir = ENV)

    Mega2GenABELtst()

    mega
  
}

cat("First, call data(srdta) to expose genabel's data object.", "\n")
cat("Second, call dmpPed() to write the srdta to srdta.ped/.map", "\n")
cat("Then at the command line type: 'mega2 MEGA2.BATCH.srdta' to get a .db file.", "\n")
cat("The 'mega2 MEGA2.BATCH.srdta' file is in the mega2tutorial directory.", "\n")
cat("Finally, run the go() function.", "\n")
