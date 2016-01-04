/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2016 Robert Baron, Charles P. Kollar,
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

===========================================================================
*/

/* Logs from output_file_names.c

   Revision 1.5.2.13.4.10.2.3  2002/05/23 19:40:42  nandita
   1) Rewrote the field_widths subroutine with each field separately
   calculated in separate functions.

   Revision 1.5.2.13.4.5  2001/08/28 16:59:40  nandita
   1) Added unique id field to the field_widths routine.
   2) Imported bug fix for SAGE fcor from bugfixes2_2 branch.

   Logs from write_files.c

   Revision 1.10.2.12.2.2.2.8.2.5  2002/06/10 20:57:15  nandita
   New function to write out ID-Unique ID map
*/


#include <stdio.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "grow_string_ext.h"
#include "output_routines_ext.h"
#include "utils_ext.h"
/*
        grow_string_ext.h:  grow
    output_routines_ext.h:  create_formats
              utils_ext.h:  imax summary_time_stamp
*/



/* The following functions are to find out how many columns
   are necessary to output data such as ids, alleles etc.
*/

void locus_name_width(linkage_locus_top *TTop, int *LocWidth)


{
    int i, new_loc_wid;

    *LocWidth = (int) strlen(TTop->Locus[0].LocusName);
    for (i = 1; i < TTop->LocusCnt; i++) {
        new_loc_wid= (int) strlen(TTop->Locus[i].LocusName);
        *LocWidth = ((new_loc_wid > *LocWidth)? new_loc_wid : *LocWidth);
    }
    return;
}



void allele_number_width(linkage_locus_top *TTop, int *AlleleWidth)

{
    int i, new_width;
    int all_cnt;

    all_cnt = TTop->Locus[TTop->PhenoCnt].AlleleCnt;
    if (all_cnt) {
        *AlleleWidth= (int) floor(log10((double)all_cnt)) + 1;
    } else {
        *AlleleWidth=2;
    }
    for (i = TTop->PhenoCnt + 1; i < TTop->LocusCnt; i++) {
        if (TTop->Locus[i].AlleleCnt > 0) {
            new_width =
                (int) floor(log10((double)TTop->Locus[i].AlleleCnt)) + 1;
            *AlleleWidth = ((new_width > *AlleleWidth)? new_width : *AlleleWidth);
        }
    }
    return;
}

void unique_id_width(linkage_ped_top *TTop, int *UidWidth)

{
    int i, j;
    int new_width;

    *UidWidth=1;
    for (i = 0; i < TTop->PedCnt; i++) {
        for (j=0; j < TTop->Ped[i].EntryCnt; j++){
            new_width = (int) strlen(TTop->Ped[i].Entry[j].UniqueID);
            *UidWidth = ((new_width > *UidWidth)? new_width : *UidWidth);
        }
    }
    (*UidWidth)++;
    return;
}

void person_name_width(linkage_ped_top *TTop, int *PidWidth)

{

    int i, j;
    int new_width = 0; //silly compiler

    *PidWidth=1;
    for (i = 0; i < TTop->PedCnt; i++) {
        for (j=0; j < TTop->Ped[i].EntryCnt; j++){
            new_width = (int) strlen(TTop->Ped[i].Entry[j].OrigID);
            *PidWidth = ((new_width > *PidWidth)? new_width : *PidWidth);
        }
    }

    (*PidWidth)++;
    return;
}

void person_pre_name_width(linkage_ped_top *TTop, int *PidWidth)

{

    int i, j;
    int new_width = 0; //silly compiler

    *PidWidth=1;
    for (i = 0; i < TTop->PedCnt; i++) {
        for (j=0; j < TTop->Ped[i].EntryCnt; j++){
            new_width = (int) strlen(TTop->Ped[i].Entry[j].PerPre);
            *PidWidth = ((new_width > *PidWidth)? new_width : *PidWidth);
        }
    }

    (*PidWidth)++;
    return;
}

void person_id_width(linkage_ped_top *TTop, int *PidWidth)

{
    int i, j;
    int new_width = 0; // silly compiler

    *PidWidth=1;
    for (i = 0; i < TTop->PedCnt; i++) {
        for (j=0; j < TTop->Ped[i].EntryCnt; j++){
            new_width =
                (int) floor(log10((double)TTop->Ped[i].Entry[j].ID))+ 1;
        }
        *PidWidth = ((new_width > *PidWidth)? new_width : *PidWidth);
    }

    (*PidWidth)++;
    return;
}

void ped_name_width(linkage_ped_top *TTop, int *FidWidth)

{
    int i, new_width;
    *FidWidth=1;
    for (i = 0; i < TTop->PedCnt; i++) {
        new_width = (int) strlen(TTop->Ped[i].Name);
        *FidWidth = ((new_width > *FidWidth)? new_width : *FidWidth);
    }

    (*FidWidth)++;
    return;
}

void ped_pre_name_width(linkage_ped_top *TTop, int *FidWidth)

{
    int i, new_width;
    *FidWidth=1;
    for (i = 0; i < TTop->PedCnt; i++) {
        new_width = (int) strlen(TTop->Ped[i].PedPre);
        *FidWidth = ((new_width > *FidWidth)? new_width : *FidWidth);
    }

    (*FidWidth)++;
    return;
}

void ped_id_width(linkage_ped_top *TTop, int *FidWidth)

{
    int i, new_width;
    *FidWidth=1;
    for (i = 0; i < TTop->PedCnt; i++) {
        new_width = (int) floor(log10((double)TTop->Ped[i].Num)) + 1;
        *FidWidth = ((new_width > *FidWidth)? new_width : *FidWidth);
    }

    (*FidWidth)++;
    return;
}

/*----------------------------------------------------------------------------
  Function to compute necessary field widths for the TCL parameter files
  Modified: 4/7/00 to include Orig IDs
  Modified: 7/30/01, shortened the widths returned by strlen() by 1.0,
  this was overcounting by 1, whereas the log10 method of finding field_widths
  gives the exact length:
  e.g. floor(log10(1000)) +1 = 4,
  floor(log10(999)) +1 = 3
  returns 1 more than the width in order to leave one column gap.
  ---------------------------------------------------------------------------- */

void field_widths(linkage_ped_top *TTop, linkage_locus_top *LTop,
		  int *FidWidth, int *PidWidth, int *AlleleWidth,
		  int *LocWidth)
{
    if (FidWidth != NULL) {
        if (OrigIds[1] == 2 || OrigIds[1] == 4) {
            ped_name_width(TTop, FidWidth);
        } else if (OrigIds[1] == 6) {
            ped_pre_name_width(TTop, FidWidth);
        } else {
            ped_id_width(TTop, FidWidth);
        }
    }

    if (PidWidth != NULL) {
        switch(OrigIds[0]) {
        case 1:
        case 2:
            person_name_width(TTop, PidWidth);
            break;

        case 3:
        case 4:
            unique_id_width(TTop, PidWidth);
            break;

        case 5:
        default:
            person_id_width(TTop, PidWidth);
            break;

        case 6:
            person_pre_name_width(TTop, PidWidth);
            break;
        }
    }

    if (LocWidth != NULL) {
        locus_name_width(LTop, LocWidth);
    }

    if (AlleleWidth != NULL) {
        allele_number_width(LTop, AlleleWidth);
    }
}

void write_entry_record(FILE *fp, linkage_ped_top *Top, int ped, int per, int id,
                        int fid, int pid, int fpre, int ppre, int fnam, int pnam, int uid,
                        int oped, int oper)
{
    linkage_ped_tree *Ped   = &Top->Ped[ped];
    linkage_ped_rec *Entry  = &Ped->Entry[per];

    if (1) {
        fprintf(fp,  "%-*d ", fid , Ped->Num);     fprintf(fp, " %-*d ",  pid,  id);
    }

    if (Top->OrigIds == 1) {
        fprintf(fp, " %-*s ", fnam, Ped->Name);    fprintf(fp, " %-*s ",  pnam, Entry->OrigID);
    }

//  fprintf(fp, " %-*s ", fpre, Ped->PedPre);      fprintf(fp, " %-*s ",  ppre, Entry->PerPre);

    if (Top->UniqueIds == 1)
        fprintf(fp, " %-*s ", fnam, Ped->Name);    fprintf(fp, " %-*s ",  uid,  Entry->UniqueID);

    /* print the loop id */
    if ((Entry->loopbreakers != NULL) ||
        (Entry->OrigProband > 1)) {
	fprintf(fp, "      %2d", Entry->OrigProband);
    } else {
        fprintf(fp, "        ");
    }

    prID_ped(fp, ped, oped, Ped, "  ", " ");       prID_per(fp, oper, Entry, " ", "\n");
}

const char *prID_ped_type(int Id)
{
    switch (Id) {
        case 1:
        if (pedfile_type == POSTMAKEPED_PFT) {
            return ("Pedigree number (#1)");
        }

    case 2:
        if (pedfile_type == POSTMAKEPED_PFT) {
            return ("Ped: field (#2)");
        } else {
            return ("Premakeped pedigree number (#2).");
        }

    case 3:
        return ("Renumbered consecutively (#3)");

    case 4:
        /* This is only for nuclear pedigrees */
        return ("With extensions e.g. 1_2 etc. (#4)");

    case 5:
        /* This is a different type of an extension which is
           purely numeric */
        return ("With multipliers e.g 1002 etc. (#5)");

    case 6:
        return ("Original pre makeped pedigree (#6)");

    default:
        return ("prID_ped_type()");

    }
#if 0
    switch(Id) {
    case 1:  // Num
        return ("Post Makeped consecutive numbered (#1)");
    case 2:  // Name
        return ("Post Makeped Name(#2)");
    case 3:  // ped+1
        return ("Pre  Makeped consecutive numbered (#3)");

    case 4:  // Name (2)
        return ("Post Makeped Name (#4) (Nuclear fam)");
    case 5:  // Num
        return ("Post Makeped consecutive numbered (#5) (Nuclear fam)");

    case 6:  // PedPre
        return ("Pre  Makeped Name (#6)");

    default:  // Num
        return ("Post Makeped consecutive numbered (#?)");
    }
#endif
}

const char *prID_per_type(int Id)
{
    switch(Id) {
    case 1:
        return ("Individual id (#1)");

    case 2:
        return ("Per: field (#2)");

    case 3:
        return ("ID: field (#3)");

    case 4:
        return ("Unique id e.g 1_2, 1=ped, 2=ind (#4)");

    case 5:
        return ("Renumber consecutively in pedigree (#5)");

    case 6:
        return ("Original pre makeped person (#6)");

    default:
        return ("prID_per_type()");
    }
#if 0
    switch(Id) {
    case 1:  // OrigID
    case 2:
        return ("Post Makeped consecutive numbered (#1)");
    case 3:  // UniqueID
    case 4:
        return ("Pre  Makeped Unique ID (#3) [ped_per]");
    case 5:  // ID
        return ("Pre  Makeped consecutive numbered (#5)");
    case 6:  //PerPre
        return ("Pre  Makeped ID (#6)");
    default:  // ID
        return ("Post Makeped consecutive numbered (#?)");
    }
#endif
}

void prID_ped(FILE *fp, int ped, const char *format, linkage_ped_tree *pedp)
{
    if (OrigIds[1] == 2) {
        fprintf(fp, format, pedp->Name);
    } else if (OrigIds[1] == 3) {
        fprintf(fp, format, ped+1);
    } else if (OrigIds[1] == 4) {
        fprintf(fp, format, pedp->Name);
    } else if (OrigIds[1] == 6) {
        fprintf(fp, format, pedp->PedPre);
    } else {
        fprintf(fp, format, pedp->Num);
    }
}

void prID_per(FILE *fp, const char *format, linkage_ped_rec *tpme)
{
    if (OrigIds[0] == 1) {
        fprintf(fp, format, tpme->OrigID);
    } else if (OrigIds[0] == 2) {
        fprintf(fp, format, tpme->OrigID);
    } else if (OrigIds[0] == 3) {
        fprintf(fp, format, tpme->UniqueID);
    } else if (OrigIds[0] == 4) {
        fprintf(fp, format, tpme->UniqueID);
    } else if (OrigIds[0] == 6) {
        fprintf(fp, format, tpme->PerPre);
    } else {
        fprintf(fp, format, tpme->ID);
    }
}

void prID_fam(FILE *fp, const char *format, linkage_ped_rec *tpme, linkage_ped_rec *tpe)
{
    prID_per(fp, format, tpme);

    if (tpme->Father>0) {
        if (OrigIds[0] == 1) {
            fprintf(fp, format, tpe[tpme->Father-1].OrigID);
            fprintf(fp, format, tpe[tpme->Mother-1].OrigID);
        } else if (OrigIds[0] == 2) {
            fprintf(fp, format, tpe[tpme->Father-1].OrigID);
            fprintf(fp, format, tpe[tpme->Mother-1].OrigID);
        } else if (OrigIds[0] == 3) {
            fprintf(fp, format, tpe[tpme->Father-1].UniqueID);
            fprintf(fp, format, tpe[tpme->Mother-1].UniqueID);
        } else if (OrigIds[0] == 4) {
            fprintf(fp, format, tpe[tpme->Father-1].UniqueID);
            fprintf(fp, format, tpe[tpme->Mother-1].UniqueID);
        } else if (OrigIds[0] == 6) {
            fprintf(fp, format, tpe[tpme->Father-1].PerPre);
            fprintf(fp, format, tpe[tpme->Mother-1].PerPre);
        } else {
            fprintf(fp, format, tpe[tpme->Father-1].ID);
            fprintf(fp, format, tpe[tpme->Mother-1].ID);
        }
    }
}

void prID_rel(FILE *fp, const char *format, int key, linkage_ped_rec *tpe)
{
    if (key != 0) {
        if (OrigIds[0] == 1)
            fprintf(fp, format, tpe[key-1].OrigID);
        else if (OrigIds[0] == 2)
            fprintf(fp, format, tpe[key-1].OrigID);
        else if (OrigIds[0] == 3)
            fprintf(fp, format, tpe[key-1].UniqueID);
        else if (OrigIds[0] == 4)
            fprintf(fp, format, tpe[key-1].UniqueID);
        else if (OrigIds[0] == 6)
            fprintf(fp, format, tpe[key-1].PerPre);
        else
            fprintf(fp, format, tpe[key-1].ID);
    } else
        fprintf(fp, format, "0");
}

void prID_ped(FILE *fp, int ped, int w, linkage_ped_tree *pedp, const char *beg, const char *end)
{
    if (OrigIds[1] == 2) {
        fprintf(fp, "%s%*s%s", beg, w, pedp->Name, end);
    } else if (OrigIds[1] == 3) {
        fprintf(fp, "%s%*d%s", beg, w, ped+1, end);
    } else if (OrigIds[1] == 4) {
        fprintf(fp, "%s%*s%s", beg, w, pedp->Name, end);
    } else if (OrigIds[1] == 6) {
        fprintf(fp, "%s%*s%s", beg, w, pedp->PedPre, end);
    } else {
        fprintf(fp, "%s%*d%s", beg, w, pedp->Num, end);
    }
}

void prID_per(FILE *fp, int w, linkage_ped_rec *tpme, const char *beg, const char *end)
{
    if (OrigIds[0] == 1) {
        fprintf(fp, "%s%*s%s", beg, w, tpme->OrigID, end);
    } else if (OrigIds[0] == 2) {
        fprintf(fp, "%s%*s%s", beg, w, tpme->OrigID, end);
    } else if (OrigIds[0] == 3) {
        fprintf(fp, "%s%*s%s", beg, w, tpme->UniqueID, end);
    } else if (OrigIds[0] == 4) {
        fprintf(fp, "%s%*s%s", beg, w, tpme->UniqueID, end);
    } else if (OrigIds[0] == 6) {
        fprintf(fp, "%s%*s%s", beg, w, tpme->PerPre, end);
    } else {
        fprintf(fp, "%s%*d%s", beg, w, tpme->ID, end);
    }
}

void prID_fam(FILE *fp, int w, linkage_ped_rec *tpme, linkage_ped_rec *tpe, const char *beg, const char *end)
{
    prID_per(fp, w, tpme, beg, end);

    if (tpme->Father>0) {
        if (OrigIds[0] == 1) {
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Father-1].OrigID, end);
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Mother-1].OrigID, end);
        } else if (OrigIds[0] == 2) {
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Father-1].OrigID, end);
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Mother-1].OrigID, end);
        } else if (OrigIds[0] == 3) {
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Father-1].UniqueID, end);
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Mother-1].UniqueID, end);
        } else if (OrigIds[0] == 4) {
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Father-1].UniqueID, end);
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Mother-1].UniqueID, end);
        } else if (OrigIds[0] == 6) {
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Father-1].PerPre, end);
            fprintf(fp, "%s%*s%s", beg, w, tpe[tpme->Mother-1].PerPre, end);
        } else {
            fprintf(fp, "%s%*d%s", beg, w, tpe[tpme->Father-1].ID, end);
            fprintf(fp, "%s%*d%s", beg, w, tpe[tpme->Mother-1].ID, end);
        }
    }
}

void prID_rel(FILE *fp, int w, int key, linkage_ped_rec *tpe, const char *beg, const char *end)
{
    if (key != 0) {
        if (OrigIds[0] == 1)
            fprintf(fp, "%s%*s%s", beg, w, tpe[key-1].OrigID, end);
        else if (OrigIds[0] == 2)
            fprintf(fp, "%s%*s%s", beg, w, tpe[key-1].OrigID, end);
        else if (OrigIds[0] == 3)
            fprintf(fp, "%s%*s%s", beg, w, tpe[key-1].UniqueID, end);
        else if (OrigIds[0] == 4)
            fprintf(fp, "%s%*s%s", beg, w, tpe[key-1].UniqueID, end);
        else if (OrigIds[0] == 6)
            fprintf(fp, "%s%*s%s", beg, w, tpe[key-1].PerPre, end);
        else
            fprintf(fp, "%s%*d%s", beg, w, tpe[key-1].ID, end);
    } else
        fprintf(fp, "%s%*s%s", beg, w, "0", end);
}

/* write a mapping of the input pedigree, person, UniqueIDs to output ids*/
void write_key_file(char *ID_file, linkage_ped_top *Top)
{
    FILE *fp;
    int ped, per, lb;
    int fid, pid, uid, fnam, pnam, oped, oper, fpre, ppre;
    char format[20];
    char header1[MAX_NAMELEN], header2[MAX_NAMELEN];
    char headern[MAX_NAMELEN];
    int h1len;

    unique_id_width(Top, &uid);
    person_name_width(Top, &pnam);
    ped_name_width(Top, &fnam);

    person_id_width(Top, &pid);
    ped_id_width(Top, &fid);

    person_pre_name_width(Top, &fpre);
    ped_pre_name_width(Top, &ppre);

    uid = imax(uid, (int) strlen("ID:"));
    pid = imax(pid, (int) strlen("NPer"));
    fid = imax(fid, (int) strlen("NPed"));
    fnam = imax(fnam, (int) strlen("OPed"));
    pnam = imax(pnam, (int) strlen("OPer"));

    if (OrigIds[0] == 1 || OrigIds[0] == 2) {
        person_name_width(Top, &oper);
    } else if (OrigIds[0] == 3 || OrigIds[0] == 4) {
        unique_id_width(Top, &oper);
    } else {
        person_id_width(Top, &oper);
    }

    if (OrigIds[1] == 2 || OrigIds[1] == 4) {
        ped_name_width(Top, &oped);
    } else {
        ped_id_width(Top, &oped);
    }
    oper = imax(oper, (int) strlen("Person"));
    oped = imax(oped, (int) strlen("Pedigree"));

    if ((fp = fopen(ID_file, "w")) == NULL) {
        errorvf("Unable to open Mega2 ID file '%s' for writing.\n", ID_file);
        EXIT(FILE_WRITE_ERROR);
    }
    summary_time_stamp(mega2_input_files, fp, "");

    h1len=fid+pid+3;
    if (Top->OrigIds) {
        h1len += fnam + pnam + 4;
    }
    if (Top->UniqueIds) {
        h1len += fnam + uid + 4;
    }
    h1len += 9;

    fprintf(fp, "Output Pedigree: %s\n",   prID_ped_type(OrigIds[1]));
    fprintf(fp, "Output Person:   %s\n\n", prID_per_type(OrigIds[0]));

    sprintf(format, "%%-%ds", h1len);
    sprintf(header1, format, "INPUT");

    sprintf(format, "%%-%ds", pid + fid + 3);
    sprintf(headern, format, "Numeric");

    sprintf(format, "%%-%ds ", fid);
    sprintf(header2, format, "NPed");

    sprintf(format, " %%-%ds ", pid);
    grow(header2, format, "NPer");


    if (Top->OrigIds) {
        sprintf(format, " %%-%ds ", pnam + fnam + 2);
        grow(headern, format, "Original");

        sprintf(format, " %%-%ds ", fnam);
        grow(header2, format, "OPed");
        sprintf(format, " %%-%ds ", pnam);
        grow(header2, format, "OPer");
    }

    if (Top->UniqueIds) {
        sprintf(format, " %%-%ds ", pnam + uid + 4);
        grow(headern, format, "Unique");

        sprintf(format, " %%-%ds ", fnam);
        grow(header2, format, "Ped");
        sprintf(format, " %%-%ds ", uid);
        grow(header2, format, "ID:");
    }

//    h1len += (int) strlen(header1) + oped + oper + 4;

//    sprintf(format, "%%-%ds (%%d/%%d)", h1len);
#ifndef HIDESTATUS
    sprintf(format, "%%-%ds", (int) strlen(" OUTPUT")+3);
    grow(header1, format, " OUTPUT");
#else
    sprintf(format, "%%-%ds (%%d/%%d)", (int) strlen(" OUTPUT")+3);
    grow(header1, format, " OUTPUT", OrigIds[1], OrigIds[0]);
#endif

    strcat(header2, " Loop id ");

    sprintf(format, " %%-%ds ", oped);
    grow(header2, format, "Pedigree");
    sprintf(format, " %%-%ds ", oper);
    grow(header2, format, "Person");

    fprintf(fp, "%s\n", header1);
    fprintf(fp, "%s\n", headern);
    fprintf(fp, "%s\n", header2);
    for(ped=0; ped < Top->PedCnt; ped++) {
        for(per=0; per < Top->Ped[ped].EntryCnt; per++) {

            write_entry_record(fp, Top, ped, per, Top->Ped[ped].Entry[per].ID,
                               fid, pid, fpre, ppre, fnam, pnam, uid, -oped, -oper);

            if (Top->Ped[ped].Entry[per].loopbreakers == NULL) {
                continue;
            }
            lb=0;

            while(Top->Ped[ped].Entry[per].loopbreakers[lb] != -99) {
                /* copy the above information so many times */
                write_entry_record(fp, Top, ped, per, Top->Ped[ped].Entry[per].loopbreakers[lb],
                                   fid, pid, fpre, ppre, fnam, pnam, uid, -oped, -oper);
                lb++;
            }
        }
    }
    fclose(fp);
    return;
}

/* write a mapping of the nuclear pedigree and person to original pedigrees */
void write_nuc_ped_key_file(char *ID_file,
                            linkage_ped_top *Top, linkage_ped_top *NukeTop)

{
    FILE *fp;
    int exped, exper, ped, per;
    int fid, pid, uid, fnam, pnam, oped, oper, fpre, ppre;
    char format[20];
    char header1[MAX_NAMELEN], header2[MAX_NAMELEN];
    char headern[MAX_NAMELEN];
    int h1len;

    /* These fields are for the original pedigrees */

    unique_id_width(Top, &uid);
    person_name_width(Top, &pnam);
    ped_name_width(Top, &fnam);

    person_id_width(Top, &pid);
    ped_id_width(Top, &fid);

    person_pre_name_width(Top, &fpre);
    ped_pre_name_width(Top, &ppre);

    uid = imax(uid, (int) strlen("ID:"));
    pid = imax(pid, (int) strlen("NPer"));
    fid = imax(fid, (int) strlen("NPed"));
    fnam = imax(fnam, (int) strlen("OPed"));
    pnam = imax(pnam, (int) strlen("OPer"));

    field_widths(NukeTop, NULL, NULL, &oper, NULL, NULL);

    if (OrigIds[1] == 4) {
        ped_name_width(NukeTop, &oped);
    } else  {
        ped_id_width(NukeTop, &oped);
    }
    oper = imax(oper, (int) strlen("Person"));
    oped = imax(oped, (int) strlen("Pedigree"));

    if ((fp = fopen(ID_file, "a")) == NULL) {
        fprintf(fp, "Unable to open Mega2 keys file %s for writing.\n",
                ID_file);
        exit(FILE_WRITE_ERROR);
    }

    summary_time_stamp(mega2_input_files, fp, "");

    h1len=fid+pid+3;
    if (Top->OrigIds) {
        h1len += fnam + pnam + 4;
    }
    if (Top->UniqueIds) {
        h1len += fnam + uid + 4;
    }

    sprintf(format, "%%-%ds", h1len);
    sprintf(header1, format, "INPUT");

    sprintf(format, "%%-%ds", pid + fid + 3);
    sprintf(headern, format, "Numeric");

    sprintf(format, "%%-%ds ", fid);
    sprintf(header2, format, "NPed");


    sprintf(format, " %%-%ds ", pid);
    grow(header2, format, "NPer");

    if (Top->OrigIds) {
        sprintf(format, " %%-%ds ", pnam + fnam + 2);
        grow(headern, format, "Original");

        sprintf(format, " %%-%ds ", fnam);
        grow(header2, format, "OPed");
        sprintf(format, " %%-%ds ", pnam);
        grow(header2, format, "OPer");
    }

    if (Top->UniqueIds) {
        sprintf(format, " %%-%ds ", pnam + uid + 4);
        grow(headern, format, "Unique");

        sprintf(format, " %%-%ds ", fnam);
        grow(header2, format, "Ped");
        sprintf(format, " %%-%ds ", uid);
        grow(header2, format, "ID:");
    }
    grow(headern, "%s", "Nuclear Family");

#ifndef HIDESTATUS
    sprintf(format, "%%-%ds", (int) strlen(" OUTPUT")+3);
    grow(header1, format, " OUTPUT");
#else
    sprintf(format, "%%-%ds (%%d/%%d)", (int) strlen(" OUTPUT")+3);
    grow(header1, format, " OUTPUT", OrigIds[1], OrigIds[0]);
#endif
    strcat(header2, " Pedigree  Person");

    fprintf(fp, "%s\n", header1);
    fprintf(fp, "%s\n", headern);
    fprintf(fp, "%s\n", header2);

    for(ped=0; ped < NukeTop->PedCnt; ped++) {
        for(per=0; per < NukeTop->Ped[ped].EntryCnt; per++) {
            exped = NukeTop->Ped[ped].Entry[per].ext_ped_num - 1;
            exper = NukeTop->Ped[ped].Entry[per].ext_per_num - 1;
            /* print the linkage ids */
            if (1) {
                sprintf(format, "%%-%dd ", fid);
                fprintf(fp, format, Top->Ped[exped].Num);
                sprintf(format, " %%-%dd ", pid);
                fprintf(fp, format, Top->Ped[exped].Entry[exper].ID);
            }

            if (Top->OrigIds == 1) {
                sprintf(format, " %%-%ds ", fnam);
                fprintf(fp, format, Top->Ped[exped].Name);
                sprintf(format, " %%-%ds ", pnam);
                fprintf(fp, format, Top->Ped[exped].Entry[exper].OrigID);
            }

            if (Top->UniqueIds == 1) {
                sprintf(format, " %%-%ds ", fnam);
                fprintf(fp, format, Top->Ped[exped].Name);
                sprintf(format, " %%-%ds ", uid);
                fprintf(fp, format, Top->Ped[exped].Entry[exper].UniqueID);
            }

            /* print the output pedigree and person values */
/*
            if (OrigIds[1] == 2 || OrigIds[1] == 4) {
                sprintf(format, " %%-%ds ", oped);
                fprintf(fp, format, Top->Ped[ped].Name);
            }
            else if (OrigIds[1] == 3) {
                sprintf(format, " %%-%dd ", oped);
                fprintf(fp, format, ped+1);
*/
            if (OrigIds[1] == 4) {
                sprintf(format, " %%-%ds ", oped);
                fprintf(fp, format, NukeTop->Ped[ped].Name);
            } else {
                sprintf(format, " %%-%dd ", oped);
                fprintf(fp, format, NukeTop->Ped[ped].Num);
            }

            if (OrigIds[0] == 3 || OrigIds[0] == 4) {
                sprintf(format, " %%-%ds ", oper);
                fprintf(fp, format, NukeTop->Ped[ped].Entry[per].UniqueID);
            } else if (OrigIds[0] == 1) {
                sprintf(format, " %%-%ds ", oper);
                fprintf(fp, format, NukeTop->Ped[ped].Entry[per].OrigID);
            } else if (OrigIds[0] == 2) {
                sprintf(format, " %%-%ds ", oper);
                fprintf(fp, format, NukeTop->Ped[ped].Entry[per].PerPre);
            } else {
                sprintf(format, " %%-%dd ", oper);
                fprintf(fp, format, NukeTop->Ped[ped].Entry[per].ID);
            }

            /* write a newline  */
            fprintf(fp, "\n");

        }
    }
    fclose(fp);
    return;
}

//
// This part taken from  mendel files made into a separate function,
// should be called after field_widths are set
//
// The routine create_formats() was split up into two parts so that things like
// "FID_IID" can be generated. 
void create_formats_no_space(char *fformat, char *pformat)
{
    if (fformat != NULL)
        strcpy(fformat, (OrigIds[1] == 2 || OrigIds[1] == 4 || OrigIds[1] == 6 ? "%s" : "%d"));

    if (pformat != NULL)
        strcpy(pformat, ( (OrigIds[0] >= 1 && OrigIds[0] <= 4) || OrigIds[0] == 6 ? "%s" : "%d"));
}

// with a trailing space. This was the original...
void create_formats(const int famwid, const int perwid, char *fformat, char *pformat)
{
    if (fformat != NULL) {
        if (OrigIds[1] == 2 || OrigIds[1] == 4) {
            sprintf(fformat, "%%%ds ", famwid);
        } else if (OrigIds[1] == 6) {
            sprintf(fformat, "%%%ds ", famwid);
        } else {
            sprintf(fformat, "%%%dd ", famwid);
        }
    }

    if (pformat != NULL) {
        if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
            sprintf(pformat, "%%%ds ", perwid);
        } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
            sprintf(pformat, "%%%ds ", perwid);
        } else if (OrigIds[0] == 6) {
            sprintf(pformat, "%%%ds ", perwid);
        } else {
            sprintf(pformat, "%%%dd ", perwid);
        }
    }
}
