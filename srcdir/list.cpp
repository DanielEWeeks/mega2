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

#include <stdio.h>
#include <string.h>

#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
/*
     error_messages_ext.h:  my_malloc
*/



#ifdef USE_PROTOS
list_rec *new_list(void);
list_entry *new_list_entry(void);
list_entry *make_list_entry(list_entry *Next, list_data NewData);
list_entry *append_to_list_head(listhandle *List, list_data NewData);
list_entry *append_to_list_tail(listhandle *List, list_data NewData);
list_entry *insert_between_list_entries(list_entry *Before, list_entry *After, list_data NewData);
list_entry *list_insert_when_told(listhandle *List, list_data NewData, int (*insert_when_true_fun)(list_data EntryData, list_data NewData, int flag));
list_data pop_first_list_entry(listhandle *List);
list_entry *list_entry_free_all(list_entry *Entry, void (*free_data_fun)(list_data EntryData));
void list_free_all(listhandle *List, void (*free_data_fun)(list_data EntryData));
void free_all_from_list(listhandle *List, void (*free_data_fun)(list_data EntryData));
void *list_iterate(listhandle *List, void *Foo, void *(*fun)());
void list_iterate_simple(listhandle *List, int (*fun)(list_data EntryData));
void *list_iterate_with_pop(listhandle *List, void *Foo, void *(*fun)());
void list_iterate_sinple_with_pop(listhandle *List, int (*fun)(list_data EntryData));
#endif

/* Create a new list and set the head and tail pointers to NULL.
 */
list_rec *new_list(void) {
    list_rec *NewList = MALLOC(list_rec);
    if (NewList != NULL) {
        NewList->Head = NULL;
        NewList->Tail = NULL;
    }
    return NewList;
}


/* Create a new list entry and set its next and data pointers to
 * defaults.
 */
list_entry *new_list_entry(void) {
    list_entry *NewEntry = MALLOC(list_entry);
    if (NewEntry != NULL) {
        NewEntry->Next = NULL;
        NewEntry->Data = NULL;
    }
    return NewEntry;
}


/* Create a new list entry and set its next and data pointers to
 * the values given in the arguments.
 */
list_entry *make_list_entry(list_entry *Next, list_data NewData)
/*    list_entry *Next;
      list_data NewData; */
{
    list_entry *NewEntry = MALLOC(list_entry);
    if (NewEntry != NULL) {
        NewEntry->Next = Next;
        NewEntry->Data = NewData;
    }
    return NewEntry;
}


/* Make a new entry with data set to NewData and put it at the
 * head of the list.
 *
 * Return NULL if the list doesn't exist or can't create the new
 * entry, or else return a pointer to the new entry.
 */
list_entry *append_to_list_head(listhandle *List, list_data NewData)
/*   list *List;
     list_data NewData; */
{
    list_entry *Entry;
    if ((List == NULL) || ((Entry = make_list_entry(List->Head, NewData)) == NULL))
        return NULL;
    if (List->Tail == NULL) /* new list */
        List->Tail = Entry;
    List->Head = Entry;
    return Entry;
}


/* Same as above only do it at the end of the list.
 */
list_entry *append_to_list_tail(listhandle *List, list_data NewData)
/*   list *List;
     list_data NewData; */
{
    list_entry *Entry;
    if ((List == NULL) || ((Entry = make_list_entry(NULL, NewData)) == NULL))
        return NULL;
    if (List->Head == NULL) /* new list */
        List->Head = Entry;
    else List->Tail->Next = Entry;
    List->Tail = Entry;
    return Entry;
}


/* Insert a new entry with data set to NewData between two other entries.
 * If the entries are NULL this doesn't make much sense but let him do it.
 */
list_entry *insert_between_list_entries(list_entry *Before, list_entry *After, list_data NewData)
/* list_entry *Before;
   list_entry *After;
   list_data NewData; */
{
    list_entry *Entry = make_list_entry(After, NewData);
    if (Before != NULL) Before->Next = Entry;
    return Entry;
}


/* Go through the list and insert when you are told to.
 *
 * If insert_when_true_fun() returns 0, continue. If it returns
 * INSERT_AFTER_MODE, insert the new entry after the current entry.
 * If it returns INSERT_BEFORE_MODE (or anything else, actually) insert
 * before the current entry.
 *
 * insert_when_true_fun() should get all the information it needs from
 * it's single argument, Data.
 */
list_entry *list_insert_when_told(listhandle *List, list_data NewData, int (*insert_when_true_fun)(list_data EntryData, list_data NewData, int flag))
/*     list *List;
       list_data NewData;
       #ifdef USE_PROTOS
       int (*insert_when_true_fun)(list_data EntryData, list_data NewData, int flag);
       #else
       int (*insert_when_true_fun)( *//* list_data EntryData, list_data NewData, int flag *//* );
                                                                                               #endif
                                                                                            */
{
    list_entry *Entry, *PrevEntry = NULL, *NewEntry = NULL;
    int mode = 0;
    if (List == NULL) return NULL;
    if (List->Head == NULL) {
        if (insert_when_true_fun(NULL, NewData, LIST_START) != STOP_ITERATION)
            return append_to_list_head(List, NewData);
        return NULL;
    }
    for (Entry = List->Head; (Entry != NULL) && (mode == 0);
         PrevEntry = Entry, Entry = Entry->Next) {
        mode = insert_when_true_fun(Entry->Data, NewData, LIST_OK);
        if (mode) {
            if (mode == LIST_INSERT_AFTER_MODE) {
                if (Entry->Next == NULL)
                    NewEntry = append_to_list_tail(List, NewData);
                else
                    NewEntry = insert_between_list_entries(Entry, Entry->Next, NewData);
            } else if (mode == LIST_INSERT_BEFORE_MODE) {
                if (PrevEntry == NULL)
                    NewEntry = append_to_list_head(List, NewData);
                else
                    NewEntry = insert_between_list_entries(PrevEntry, Entry, NewData);
            } else return NULL;
        }
    }
    if (mode == 0) {
        if (insert_when_true_fun(NULL, NewData, LIST_ENDED))
            return append_to_list_tail(List, NewData);
        return NULL;
    }

    return NewEntry;
}


/* Free the first entry and return a pointer to the data.
 */
list_data pop_first_list_entry(listhandle *List)
/* list *List; */
{
    list_data OldData;
    list_entry *NewHead;
    if ((List == NULL) || (List->Head == NULL)) return NULL;
    OldData = List->Head->Data;
    NewHead = List->Head->Next;
    free(List->Head);
    List->Head = NewHead;
    if (List->Head == NULL) List->Tail = NULL;
    return OldData;
}


/* Given a function for freeing the data, free the entry and the data.
 */
/*  ORIGINAL VERSION:
    void list_entry_free_all(Entry, free_data_fun)
    list_entry *Entry;
    #ifdef USE_PROTOS
    void (*free_data_fun)(list_data EntryData);
    #else
    void (*free_data_fun)();
    #endif
    {
    if (free_data_fun != NULL) free_data_fun(Entry->Data);
    free(Entry);
    }
    END OF ORIGINAL VERSION */

list_entry *list_entry_free_all(list_entry *Entry, void (*free_data_fun)(list_data EntryData))
/*     list_entry *Entry;
       #ifdef USE_PROTOS
       void (*free_data_fun)(list_data EntryData);
       #else
       void (*free_data_fun)( *//* list_data EntryData */ /*);
                                                            #endif
                                                          */
{
    list_entry *temp;
    temp = Entry->Next;
    if (free_data_fun != NULL) free_data_fun(Entry->Data);
    free(Entry);
    return(temp);
}


/* Free the whole list. Must have a function to free the data.
 */
void list_free_all(listhandle *List, void (*free_data_fun)(list_data EntryData))
/*     list *List;
       #ifdef USE_PROTOS
       void (*free_data_fun)(list_data EntryData);
       #else
       void (*free_data_fun)( *//* list_data EntryData *//* );
                                                            #endif
                                                         */
{
    list_entry *Entry;
    if (List == NULL) return;
/* ORIGINAL CODE:  Error, since Entry will be removed before Entry->Next
   is accessed
   for (Entry = List->Head; Entry != NULL; Entry = Entry->Next) {
   list_entry_free_all(Entry, free_data_fun);
*/
    Entry = List->Head;
    while (Entry != NULL)
        Entry = list_entry_free_all(Entry, free_data_fun);
    free(List);
}


/* Free all entries and data and clear head
 */
void free_all_from_list(listhandle *List, void (*free_data_fun)(list_data EntryData))
{
    list_entry *Entry;
    if (List == NULL) return;
/* ORIGINAL CODE:  Error, since Entry will be removed before Entry->Next
   is accessed
   for (Entry = List->Head; Entry != NULL; Entry = Entry->Next) {
   list_entry_free_all(Entry, free_data_fun);
*/
    Entry = List->Head;
    while (Entry != NULL)
        Entry = list_entry_free_all(Entry, free_data_fun);
    List->Head = NULL;
    List->Tail = NULL;
}


/* Apply a function to every element of the list in turn. Support passing
 * of some user data structure between calls. Return the void* returned
 * by the fun() that quit.
 */
void *list_iterate(listhandle *List, void *Foo, void *(*fun)(list_data EntryData, void *Foo))
{
    list_entry *Entry;
    void *ret = NULL;

    if (List == NULL) return NULL;
    Entry = List->Head;
    while ((Entry != NULL) && ((ret = fun(Entry->Data, Foo)) == NULL))
        Entry = Entry->Next;
    return ret;
}


/* Same as above only don't bother passing around a data structure,
 * and don't return anything.
 */
void list_iterate_simple(listhandle *List, int (*fun)(list_data EntryData))
{
    list_entry *Entry;
    if (List == NULL) return;
    Entry = List->Head;
    while ((Entry != NULL) && (fun(Entry->Data) == 0))
        Entry = Entry->Next;
}


/* Apply a function to every element of the list in turn. Support passing
 * of some user data structure between calls. Destroy the list (but not
 * List itself).
 */
void *list_iterate_with_pop(listhandle *List, void *Foo, void *(*fun)(list_data EntryData, void *Foo))
{
    list_data Data;
    void *ret = NULL;

    if (List == NULL) return NULL;
    while ((Data = pop_first_list_entry(List), (!LIST_EMPTY(List)))
           && ((ret = fun(Data, Foo)) == NULL));
    return ret;
}


/* Same as above only don't bother passing around a data structure.
 */
void list_iterate_sinple_with_pop(listhandle *List, int (*fun)(list_data EntryData))
{
    list_data Data;
    if (List == NULL) return;
    while ((Data = pop_first_list_entry(List), (!LIST_EMPTY(List)))
           && (fun(Data) == 0));
}
