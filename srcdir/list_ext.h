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

#ifndef LIST_EXT_H
#define LIST_EXT_H

extern list_entry *append_to_list_tail(listhandle *List, list_data NewData);

extern void free_all_from_list(listhandle *List, void (*free_data_fun)(list_data EntryData));

extern void list_free_all(listhandle *List, void (*free_data_fun)(list_data EntryData));

extern list_entry *list_insert_when_told(listhandle *List, list_data NewData, int (*insert_when_true_fun)(list_data EntryData, list_data NewData, int flag));

extern void *list_iterate(listhandle *List, void *Foo, void *(*fun)(list_data EntryData, void *Foo));


extern list_rec *new_list(void);

extern list_data pop_first_list_entry(listhandle *List);


/*
extern list_entry *append_to_list_head(listhandle *List, list_data NewData);

extern list_entry *insert_between_list_entries(list_entry *Before, list_entry *After, list_data NewData);

extern list_entry *list_entry_free_all(list_entry *Entry, void (*free_data_fun)());

extern void list_iterate_simple(listhandle *List, int (*fun)());

extern void list_iterate_sinple_with_pop(listhandle *List, int (*fun)());

extern void *list_iterate_with_pop(listhandle *List, void *Foo, void *(*fun)());

extern list_entry *make_list_entry(list_entry *Next, list_data NewData);

extern list_entry *new_list_entry(void);

*/

#endif
