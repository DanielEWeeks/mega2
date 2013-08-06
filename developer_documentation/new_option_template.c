/* This is a template C file for adding in a new option */
/* We are calling the option NEWOPT */

#include "typedefs.h"
#include "fdefs.h"


/* 1) function declarations:
   These are optional but avoid compiler warnings, if functions are
   called before they are declared, and are good for clarity. */

/* 1a) Single entry function that is called from mega2.c */
void create_NEWOPT_files(linkage_ped_top **LPedTop, char *file_names[],
		       int untyped_ped_opt);

/* 1b) Local function definitions 
static f1()
static f2() etc.

*/
/* This is a function that implements the file-names menu for NEWOPT,
   arguments: */
static void NEWOPT_file_names(file_names, Top->OrigIds, Top->UniqueIds, ...);

/* This has to come before the entry level function, it allows the
   pointer to the LPedTop to be renamed for convinience */

#ifdef Top
#undef Top
#endif
#define Top (*LPedTop)

/* 2) Template for the entry function */

void create_NEWOPT_files(linkage_ped_top **LPedTop, char *file_names[],
		       int untyped_ped_opt)

{

  /* get the output file names, also use this function to set output
     related parameters, like whether an overall shell script should
     be created for all chromosome etc. */

  NEWOPT_file_names(file_names, Top->OrigIds, Top->UniqueIds, ...);
  /* ... refers to other output file related parameters */

  /* This displays a line on the screen stating that the following
     files will be created for NEWOPT */

  create_mssg(analysis);

  /* This function sets the formatting for pedigree and individual
     ids, as well as locus Ids */

  field_widths(Top, LTop, &fwid, &pwid, NULL, &mwid);
    
  /* Add in code for setting other run-realted parameters that are
     specific to this option, and if these are going to remain
     constant across all chromosomes */

  /* Now write the main loop for creating chromosome-specific files:
   */

  for (i=0; i < main_chromocnt; i++) {
    numchr=global_chromo_entries[i]; 
    
    /* This creates a list of marker loci on the selected chromosome */
    if (global_chromo_entries[i] == UNKNOWN_CHROMO) {
      get_unmapped_loci(0);
    }
    else {
      get_loci_on_chromosome(numchr);
    }
    
    replace_chr_number(&analysis, file_names, &numchr);
    omit_peds(untyped_ped_opt, Top);
    
       ...
  }

}

/* other functions that are usually relevant */
write_NEWOPT_pedigrees()
{
}

write_NEWOPT_locus()
{
}

write_NEWOPT_map()
{
}

write_NEWOPT_shell_script()
{
}
