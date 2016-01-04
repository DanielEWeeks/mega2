/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2016 Robert Baron, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  Daniel E. Weeks, and University of Pittsburgh

  This file is part of the Mega2 program, which is free software; you
  can redistribute it and/or modify it under the terms of the GNU
  General Public License as published by the Free Software Foundation;
  either version 3 of the License, or (at your option) any later
  version.

  Mega2 is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

  For further information contact:
      Daniel E. Weeks
      e-mail: weeks@pitt.edu

*/

#ifndef WRITE_SHAPEIT_ARGS_H
#define WRITE_SHAPEIT_ARGS_H

static char SHAPEIT_ARGS[] = "\n\
\n\
# To set any of these shapeit parameters, remove the # from the start of the line and\n\
# adjust the parameters to your requirements.\n\
\n\
set MoreArgs=\"\"\n\
\n\
# run shapeit in check mode; no other argument from this listing should be enabled\n\
#set MoreArgs=\"$MoreArgs -check\"\n\
\n\
# for faster results, replace 2 with the number of cores on your computer \n\
#set MoreArgs=\"$MoreArgs --thread 2\"\n\
\n\
\n\
#    ALGORITHM  PARAMETERS\n\
\n\
\n\
# burn is the number of burn-in iterations\n\
#set MoreArgs=\"$MoreArgs --burn 7\"\n\
\n\
# prune is the iterations of pruning stage\n\
#set MoreArgs=\"$MoreArgs --prune 8\"\n\
\n\
# main is the number of main iterations\n\
#set MoreArgs=\"$MoreArgs --main 20\"\n\
\n\
# set the start value for the pseudo random number generator\n\
#set MoreArgs=\"$MoreArgs --seed 123456789\"\n\
\n\
\n\
#    MODEL PARAMETERS\n\
\n\
\n\
# more conditioning states improve accuracy but increase run time\n\
#set MoreArgs=\"$MoreArgs --states 100\"\n\
\n\
# window size in Mb.  For sequence data 0.5 is better\n\
#set MoreArgs=\"$MoreArgs --window 2\"\n\
\n\
# revert to shapeit version 1 haplotype conditioning\n\
#set MoreArgs=\"$MoreArgs --model-version1\"\n\
\n\
# if there is no recombination data, use a uniform value\n\
#set MoreArgs=\"$MoreArgs --rho 0.0004\"\n\
\n\
# effective population size\n\
#set MoreArgs=\"$MoreArgs --effective-size 15000\"\n\
\n\
\n\
#    FILTER PARAMETERS\n\
\n\
\n\
# restrict start position bp (default shown)\n\
#set MoreArgs=\"$MoreArgs --input-from 0\"\n\
\n\
# restrict end position bp (default shown)\n\
#set MoreArgs=\"$MoreArgs --input-to 1e9\"\n\
\n\
# exclude markers by bp from analysis\n\
#set MoreArgs=\"$MoreArgs --exclude-snp exclusion_file\"\n\
\n\
# include markers by bp and only these markers\n\
#set MoreArgs=\"$MoreArgs --include-snp inclusion_file\"\n\
\n\
# exclude samples from analysis\n\
#set MoreArgs=\"$MoreArgs --exclude-ind exclusion_file\"\n\
\n\
# include only these samples\n\
#set MoreArgs=\"$MoreArgs --include-ind inclusion_file\"\n\
\n\
\n\
#    MISC PARAMETERS\n\
\n\
# do not analyze pedigree for duos/trios\n\
# This is necessary for complex pedigrees.  You might want to turn this off by\n\
# prepending the #\n\
set MoreArgs=\"$MoreArgs --duohmm\"\n\
\n\
\n\
#\n\
# This covers the first third of \n\
# https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gettingstarted\n\
# If you need more specialized features of shapeit, look there for the parameters to use.\n\
#\n\
";

#endif
