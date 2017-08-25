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

#ifndef ERRORNO_H
#define ERRORNO_H

#define SUCCESS 0

#define FILE_NOT_FOUND 1
#define FILE_READ_ERROR 2
#define FILE_WRITE_ERROR 3
#define INPUT_DATA_ERROR 4

#define RESERVED_WORD_ERROR 5
#define DATA_INCONSISTENCY 6
#define DATA_TYPE_ERROR 7
#define BATCH_FILE_ITEM_ERROR 8

#define EARLY_TERMINATION 9
#define OUTOF_BOUNDS_ERROR 10
#define MEMORY_ALLOC_ERROR 11
#define MEMORY_FREE_ERROR 12
#define MEMORY_REALLOC_ERROR 13
#define SYSTEM_ERROR 14
#define MEGA2_EXPIRED 15

#define LOOP_CONNECTION_ERROR 16
#define LOOP_BREAKING_ERROR 17
#define OUTPUT_FORMAT_ERROR 18

#define COMMAND_NOT_FOUND 19
#define CONNECTION_REFUSED 20
#define SERVER_NOT_FOUND 21

#define FILE_HEADER_ERROR 22
#define EXCEEDED_FIELD_MAX 23

#endif
