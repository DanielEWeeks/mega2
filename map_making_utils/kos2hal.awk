#! /bin/gawk -f
# Script to convert a Kosambi map file
# into a Haldane map file.
#
#  Expects an input file in the format:
# MarkerName1 Chr#
# 10.2
# MarkerName2 Chr#
# 5.1
# MarkerName3 Chr#
# 20.2 
# MarkerName4  Chr#
#
# where the numbers are the inter-marker distances
# in Kosambi cM.  Will convert this to a map file
# in Haldane cM.
#
# Assumes that all markers on a given chromosome are
# grouped together.  
#
BEGIN { n = 0; i = 0; j = 0; pos = 0.0;
	if (ARGC != 2)
	{
	print "USAGE: kos2hal.awk file.name"
	exit(1);
	}
	printf("Chr      Haldane cM  Name        Haldane  Theta  Kosambi\n"); 
}
{ n = n + 1
  if (n%2 == 0) 
   {
   j = j + 1;
   dist_k[j] = $1
   x = exp(4.0*(dist_k[j]/100.0));
   theta[j] = 0.50*(x - 1.0)/(x + 1.0);
   dist[j] = 100.0*(-0.50)*log(1.0 - 2.0*theta[j]);
   }
  else
  {
   i = i + 1
   name[i] = $1
   chr[i] = $2
  }
}
END { 
 for (k=1;k<=i; k++)
 {
  if (k <= j)
  {
  if (chr[k+1] == chr[k])
  printf("%-2d         %7.3f   %-10s  %7.3f %7.5f %7.3f\n",chr[k], pos, name[k],dist[k], theta[k], dist_k[k]);
  else
  printf("%-2d         %7.3f   %-10s  \n",chr[k],pos, name[k]);
  }
  else
  printf("%-2d         %7.3f   %-10s  \n",chr[k],pos, name[k]);
  if (k<=j)
   if (chr[k+1] != chr[k])
    pos = 0;
   else
    pos = pos + dist[k];
 }
}
