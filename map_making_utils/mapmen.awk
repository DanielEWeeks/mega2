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
	if (ARGC != 2)
	{
	print "USAGE: mapmen.awk file.name"
	exit(1);
	}
}
{ n = n + 1
  if (n%2 == 0) 
   {
   j = j + 1;
   dist[j] = $1
   theta[j] = 0.5*(1.0 - exp(-2.0*(dist[j]/100.0)));
   }
  else
  {
   i = i + 1
   name[i] = $1
  }
}
END { 
 printf("%-8s\n",name[1]);
 for (k=1;k<=i; k++)
 {
  if (k <= j)
  printf("       %8.4f%8.4f \n%-8s\n",theta[k],theta[k],name[k+1]);
 }
}
