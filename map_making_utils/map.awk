#! /bin/nawk -f
BEGIN { ndist = 0; print "Intermarker distances from the file ",FILENAME }
{ if (NR > 1)
  {
  pos[NR-1] = $2
  name[NR-1] = $3
  ndist += 1
  }
}
END { 
 for (i=2;i<=ndist;i++)
  printf("%10s --- %5.1f --- %10s (%5.1f)\n", name[i-1],pos[i]-pos[i-1],name[i],pos[i])
  }
