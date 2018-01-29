/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2018, University of Pittsburgh. All Rights Reserved.

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

#ifndef ALLUTIL_H
#define ALLUTIL_H

/* from fcmap.h
   Define the integer which keeps track of the number
   of tokens (format args and characters) that have
   been processed.
*/
#define FCMAP_NARGS fcmap_nargs

/*
  Define the status integer and bits and macros.
*/

#define FCMAP_STAT_TERM    1
#define FCMAP_STAT_TCHDEF  2
#define FCMAP_STAT_LCHDEF  4
#define FCMAP_TERM    (fcmap_stat & FCMAP_STAT_TERM)   /* terminate char or EOF  */
#define FCMAP_TCHDEF  (fcmap_stat & FCMAP_STAT_TCHDEF) /* terminate char defined */
#define FCMAP_LCHDEF  (fcmap_stat & FCMAP_STAT_LCHDEF) /* last character defined */

#define LCH   _fcmap_lch
#define TERM  _fcmap_term
#define NARGS fcmap_nargs
#define STAT  fcmap_stat

#define DEFTERM   fcmap_stat |= FCMAP_STAT_TERM
#define UNDEFTERM fcmap_stat = (fcmap_stat & ~FCMAP_STAT_TERM)
#define TERMDEF   FCMAP_TERM

#define DEFTCH    fcmap_stat |= FCMAP_STAT_TCHDEF
#define UNDEFTCH  fcmap_stat &= ~FCMAP_STAT_TCHDEF
#define TCHDEF    FCMAP_TCHDEF

#define DEFLCH    fcmap_stat |= FCMAP_STAT_LCHDEF
#define GETLCH    { LCH = fgetc(filep); DEFLCH; }
#define UNDEFLCH  fcmap_stat &= ~FCMAP_STAT_LCHDEF
#define LCHDEF    FCMAP_LCHDEF

#define LCH_IS_TERM  (TCHDEF && (LCH == TERM))
#define LCH_NOT_TERM (!TCHDEF || (LCH != TERM))

extern int NARGS, STAT;

/* from list.h
   list data structures */

/* Maybe it's gross for these to be macros, but wet the hay. */
#define free_list_entry(e) (free(e), e = NULL)
#define free_list(l) (free(l), l = NULL)
#define list_entry_data(e)  ((e)->Data)
#define next_list_entry(e)  ((e)->Next)
#define first_list_entry(l) ((l)->Head)
#define first_list_data(l)  ((l)->Head->Data)
#define last_list_entry(l)  ((l)->Tail)
#define last_list_data(l)   ((l)->Tail->Data)

#define LIST_INSERT_BEFORE_MODE  -1
#define LIST_INSERT_AFTER_MODE    1
#define LIST_DONT_INSERT_YET      0
#define LIST_START          -1
#define LIST_ENDED           1
#define LIST_OK              0

#define STOP_ITERATION 1

#define LIST_EMPTY(l) ((l)->Head == NULL)

/*  list data structures */
typedef struct _lrec_data {
    int PedID;
    int Proband;    /* may also contain loop information */
} lrec_data;

#define ALL_PEDS 0

typedef struct _fraction_interval {
    float Start, Increment, Stop;
} fraction_interval;


typedef struct _old_new {
    /* lwa */
    int oldID;
    int newID;

    int oldmem;
    int newmem;
    int oldped;
    int newped;
    struct _old_new *next;
} old_new;

typedef struct _affected_data{
    int AffectedCnt;
    int *Affected;    /* Affected[] is index in Entry[] */
} affected_data;

typedef struct _Untype_Info{
    int Pednumber;
    int Person;
} Untype_Info;

typedef struct Untype_List_Recs{
    Untype_Info *Data;
    struct Untype_List_Recs *Next;
} Untype_List_Rec;

typedef void *list_data;

typedef struct _list_entry {
    list_data Data;
    struct _list_entry *Next;
} list_entry;

typedef struct _list_rec {
    list_entry *Head;
    list_entry *Tail;
} list_rec;

typedef list_rec listhandle;
typedef listhandle *list_ptr;

typedef struct _program_data_rec {
    int Frac;
    listhandle *IntervalList;
} program_data_rec;

/* taken out all the extern definitions from below-Nandita */

#endif
