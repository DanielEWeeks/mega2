#! /bin/nawk -f
BEGIN { ndist = 0; }
{ if (NR > 1)
  {
  pos[NR-1] = $2
  name[NR-1] = $3
  ndist += 1
  }
}
END { 
 printf("genome\n{\n region\n {\n");
 for (i=2;i<=ndist;i++)
  printf("marker = %-10s\ndistance = %-5.2f\n", name[i-1],pos[i]-pos[i-1])
 printf("marker = %-10s\n",name[ndist]);
 printf(" }\n}");
  }
