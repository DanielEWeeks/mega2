/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

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

#ifndef OUTPUT_ROUTINES_EXT_H
#define OUTPUT_ROUTINES_EXT_H

extern void create_formats_no_space(char *fformat, char *pformat);
extern void create_formats(const int famwid, const int perwid,
			   char *fformat, char *pformat);

extern void field_widths(linkage_ped_top *TTop, linkage_locus_top *LTop,
                         int *FidWidth, int *PidWidth, int *AlleleWidth,
			 int *LocWidth);


extern void locus_name_width(linkage_locus_top *TTop, int *LocWidth);

extern void ped_id_width(linkage_ped_top *TTop, int *FidWidth);

extern void ped_name_width(linkage_ped_top *TTop, int *FidWidth);

extern void person_id_width(linkage_ped_top *TTop, int *PidWidth);

extern void prID_fam(FILE *fp, const char *format, linkage_ped_rec *tpme, linkage_ped_rec *tpe);

extern const char *prID_ped_type(int Id);
extern void prID_ped(FILE *fp, int ped, const char *format, linkage_ped_tree *pedp);

extern const char *prID_per_type(int Id);
extern void prID_per(FILE *fp, const char *format, linkage_ped_rec *tpme);

extern void prID_rel(FILE *fp, const char *format, int key, linkage_ped_rec *tpe);

extern void write_key_file(char *ID_file, linkage_ped_top *Top);

extern void write_nuc_ped_key_file(char *ID_file, linkage_ped_top *Top,
				   linkage_ped_top *NukeTop);

extern void prID_fam(FILE *fp, int w, linkage_ped_rec *tpme, linkage_ped_rec *tpe, const char *beg, const char *end);

extern void prID_ped(FILE *fp, int ped, int w, linkage_ped_tree *pedp, const char *beg, const char *end);

extern void prID_per(FILE *fp, int w, linkage_ped_rec *tpme, const char *beg, const char *end);

extern void prID_rel(FILE *fp, int w, int key, linkage_ped_rec *tpe, const char *beg, const char *end);


/*
extern void allele_number_width(linkage_locus_top *TTop, int *AlleleWidth);

extern void person_name_width(linkage_ped_top *TTop, int *PidWidth);

extern void unique_id_width(linkage_ped_top *TTop, int *UidWidth);

*/

#endif
