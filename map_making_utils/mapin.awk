#! /bin/nawk -f
#  Expects an input file in the format:
# MarkerName1
# 10.2
# MarkerName2
# 5.1
# MarkerName3
# 20.2 
# MarkerName4 
#
# where the numbers are the inter-marker distances
# in cM.  Will convert this to a map file.
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
	printf("CHROMOSOME    cM   NAME\n");
}
{ n = n + 1
  if (n%2 == 0) 
   {
   j = j + 1;
   dist[j] = $1
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
  printf("%-2d         %7.3f   %-10s  %7.3f\n",chrnum, pos, name[k],dist[k]);
  else
  printf("%-2d         %7.3f   %-10s  \n",chrnum,pos, name[k]);
  pos = pos + dist[k];
 }
}
