library(mega2genabel)

go = function() {

    print("## generate genabel for srdta")

    ENV = read.Mega2DB("srdta.db", verbose = 1)

    mega = Mega2GenABEL("xyz", envir = ENV)

    Mega2GenABELtst()

}

cat("First, type:", "\n\tdata(srdta)", "\nto expose genabel's data object.", "\n")
cat("Second, type:", "\n\tdmpPed()", "\nto write the srdta to srdta.ped/.map", "\n")
cat("Then at the command line type:", '\n\tmega2 MEGA2.BATCH.srdta', "\nto get a .db file.", "\n")
cat("Note the 'mega2 MEGA2.BATCH.srdta' file is in the mega2tutorial directory.", "\n")
cat("Finally, run the:", "\n\tgo()", "\nfunction.", "\n")
