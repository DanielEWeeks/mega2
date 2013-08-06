   Some additional nawk scripts for creating map files.

mapmarsh.awk:  Script to convert a frameworkMap file created by
	       the 'Build your own map' option of the Marshfield
	       Center for Medical Genetics web page
	       (http://www.marshmed.org/genetics/) to a Mega2-format
	       map file.  However, note that this will create a
	       map file in Kosambi cM, whereas Mega2 expects a
	       map file in Haldane cM.

mapout.awk:    Script to convert from the Mega2-format map file
	       into a format that is easier to edit.  This format
	       gives the distances between each pair of adjacent
	       markers, so we will call it the 'intermarker-format'.

kos2hal.awk:   Script to convert from an intermarker-format map
	       file in Kosambi cM to a Mega2-format file in
	       Haldane cM.

Example:
  Here I illustrate the commands I used to convert some framework
  map files downloaded from the Marshfield site (which are in
  Kosambi cM) into Mega2-format map files in Haldane cM.

  Comments are in '[]' brackets.

Script started on Fri Dec 03 16:20:33 1999
[ NOTES:
    1) Make sure you save the files from Marshfield
    as text-only.  The 'mapmarsh.awk' script will
    not work properly otherwise.
    2) You may have to edit the files from Marshfield,
    as sometime some markers are missing or alternative
    marker names are used. ]
watson%	mapmarsh.awk frameworkMap.20 20 > mapk.20
watson%	mapmarsh.awk frameworkMap.21 21 > mapk.21
watson%	mapmarsh.awk frameworkMap.22 22 > mapk.22
[So now we have three chromosome specific map files
in Kosambi cM. We will place these in one file
using the 'cat' command: ] 
watson%	cat mapk.* > map.kos
[Now we must edit the 'map.kos' file to remove
the extra internal title lines.]
watson%	vi map.kos
[ Convert this file into intermarker-format:]
watson%	mapout.awk map.kos > mapin.kos
[ Convert the intermarker-format file into
a Mega2-format file in Haldane cM: ]
watson%	kos2hal.awk mapin.kos > map.hal
[ And here it is: ]
watson%	more map.hal
Chr      Haldane cM  Name        Haldane  Theta  Kosambi
20           0.000   D20S103      10.981 0.09859   9.990
20          10.981   GATA149E11   14.146 0.12321  12.580
20          25.128   D20S851       8.916 0.08166   8.240
20          34.043   D20S604       6.707 0.06277   6.310
20          40.751   D20S470      17.011 0.14419  14.840
20          57.761   D20S478       8.904 0.08156   8.230
20          66.666   D20S481      20.622 0.16899  17.590
20          87.288   D20S480      18.243 0.15285  15.790
20         105.531   D20S171     
21           0.000   D21S1432     11.065 0.09926  10.060
21          11.065   D21S1437     13.032 0.11472  11.680
21          24.097   D21S1442     13.476 0.11813  12.040
21          37.573   D21S1440      3.858 0.03713   3.720
21          41.431   D21S2055     20.208 0.16623  17.280
21          61.640   D21S1446    
22           0.000   D22S420      17.553 0.14803  15.260
22          17.553   D22S1174     10.101 0.09146   9.250
22          27.654   D22S689       3.966 0.03813   3.820
22          31.620   D22S685       3.977 0.03823   3.830
22          35.597   D22S683      10.516 0.09484   9.600
22          46.113   D22S445     
watson%	
script done on Fri Dec 03 16:23:02 1999

Other scripts:

mapin.awk   : Script to convert an intermarker-format map file
	      (without the chromosome numbers indicated) into
	      a Mega2-format map file.  Handles only files with
	      markers from a single chromosome.
mapin2.awk  : Script to convert an intermarker-format map file
	      into a Mega2-format map file.  Handles files
	      with markers from multiple chromosomes.
mapmen.awk  : Script to convert an intermarker-format map file
	      into a Mendel 4 format map file.
mapsage.awk : Script to convert a Mega2-format map file into
	      a map file in the new SAGE format.
mapsum.awk  : Script to rearrange a Mega2-format map file
	      so that the locus name is in the first column.


NOTE:
 1) These script were developed quickly for my personal use, and
 so may not be very robust to unexpected input file formats.  Please
 check the files that are made by these scripts carefully.
 2) Under linux, you will need to change the word 'nawk' to the word
 'gawk' at the top of each file in order for these scripts to run.
 3) The scripts must be executable.  Use 'chmod +x script.awk' to make
 the 'script.awk' file executable.
