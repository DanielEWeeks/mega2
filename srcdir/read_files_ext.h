/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2017, University of Pittsburgh. All Rights Reserved.

  Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  and Daniel E. Weeks.

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

===========================================================================
*/

#ifndef READ_FILES_EXT_H
#define READ_FILES_EXT_H

extern int check_dos_file(FILE *fp);

extern int check_empty(char *file_name);

extern void check_size_dos(void);

extern void free_lpedtop_pedinfo(linkage_ped_top *PTop);

extern int plink_annot_string_aff_phen(int line, pheno_rec *locus,
                                       pheno_pedrec_data *pedrec, const char *cstatus);

extern int plink_annot_string_quant_phen(int line, pheno_rec *locus,
                                         pheno_pedrec_data *pedrec, const char *quantstr);

extern void premakeped_omit_file(linkage_ped_top *Top, const char *omitfl_name,
				 const int raw_allele);


extern int read_aff_phen(FILE *filep, int locusnm,
                         pheno_rec *locus, void *ventry,
			 record_type rec);

extern linkage_ped_top *read_linkage2(char *pedfl_name,
				      char *locusfl_name,
				      char *omitfl_name,
				      char *mapfl_name,
				      int untyped_ped_opt, analysis_type analysis);

extern linkage_ped_top *read_linkage_ped_file(FILE *filep,
					      linkage_locus_top *LTop,
                                              int *col2locus);

extern void count_Missing_Quant_consistency();

extern void annot_ignore_numbered_data(int line, linkage_locus_rec *locus,
                                       void *pedrec, int loc);

extern char *canonical_allele(const char *ra);

extern int read_numbered_data(FILE *filep, int locusnm,
                              linkage_locus_rec *locus, void *entry,
			      record_type rec, int last);


extern int read_premakeped_affec(FILE *filep, int locusnm, pheno_rec *locus,
				 pre_makeped_record *entry);

extern int read_premakeped_bin(FILE *filep, int locusnm, linkage_locus_rec *locus,
			       pre_makeped_record *entry);

extern int read_premakeped_quant(FILE *filep, int locusnm,
				 pheno_rec *locus, pre_makeped_record *entry);


extern int read_quant_phen(FILE *filep, int locusnm,
                           pheno_rec *locus, void *ventry,
			   record_type rec);

extern void omit_file_data_processing(linkage_ped_top *Top,
                                      const char *omitfl_name, FILE *omitfp,
                                      const int ped_i, const int per_i, const int loci_i,
                                      const int raw_allele);

extern void premakeped_omit_file(linkage_ped_top *Top, const char *omitfl_name, const int raw_allele);

/*
extern int read_map_file(char *mapfl_name, linkage_locus_top *LTop);

extern int read_premakeped_num(FILE *filep, int locusnm, linkage_locus_rec *locus1,
			       pre_makeped_record *entry);

*/

#endif
