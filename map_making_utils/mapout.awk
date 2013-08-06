#! /bin/nawk -f
BEGIN { ndist = 0; }
{ if (NR > 1)
  {
  pos[NR-1] = $2
  name[NR-1] = $3
  chr[NR-1] = $1
  ndist += 1
  }
}
END { 
 for (i=2;i<=ndist;i++)
  if (chr[i] == chr[i-1])
  printf("%-10s   %d\n%-5.2f\n", name[i-1],chr[i-1],pos[i]-pos[i-1])
  else
  printf("%-10s   %d\n%-5.2f\n", name[i-1],chr[i-1],-99)
 printf("%-10s   %d\n",name[ndist],chr[ndist]);
  }
