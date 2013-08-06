#! /bin/gawk -f
#  Expects an input file in the format:
# MarkerName1 Chr#
# 10.2
# MarkerName2 Chr#
# 5.1
# MarkerName3 Chr#
# 20.2 
# MarkerName4 Chr#
#
# where the numbers are the inter-marker distances
# in Haldane cM.  Will convert this to a map file
# in Kosambi cM.
#
BEGIN { n = 0; i = 0; j = 0; pos = 0.0;
	if (ARGC != 3)
	{
	print "ERROR: Please type the file name followed by the chromosome number"
	print "USAGE: mapin.awk file.name chrnum"
	exit(1);
	}
	chrnum = ARGV[2];
	delete ARGV[2];
	printf("CHROMOSOME  Kosambi cM  NAME          Kos cM  Hal cM  Theta\n");
}
{ n = n + 1
  if (n%2 == 0) 
   {
   j = j + 1;
   dist_h[j] = $1;
   x = dist_h[j]/100.0;
   theta[j] = 0.50*(1.0 - exp(-2.0*x));
   dist[j] = 100.0*0.25*log( (1 + 2.0*theta[j])/(1.0 - 2.0*theta[j]) );
   }
  else
  {
   i = i + 1
   name[i] = $1
  }
}
END { 
 for (k=1;k<=i; k++)
 {
  if (k <= j)
  printf("%-2d         %7.3f      %-10s  %7.3f %7.3f   %7.5f\n",chrnum, pos, name[k],dist[k], dist_h[k],theta[k]);
  else
  printf("%-2d         %7.3f      %-10s  \n",chrnum,pos, name[k]);
  pos = pos + dist[k];
 }
}
