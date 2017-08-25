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

#ifndef ERROR_MESSAGES_EXT_H
#define ERROR_MESSAGES_EXT_H

extern void           close_logs(void);

extern void           empty_file(const char *filename, error_type etype);

extern void           err_or_warn(FILE **errfp, FILE **logfp);

extern void           errorf(const char *messg);
extern void           errorvf(const char *fmt, ...);

extern void           errsimf(const char *messg);

extern void           errf(const char *messg, const char *shout="", int showerrf = 1);
extern void           errvf(const char *fmt, ...);

extern void           mssgf(const char *messg, int show);
extern void           mssgf(const char *messg);
extern void           mssgvf(const char *fmt, ...);
extern void           msgvf(const char *fmt, ...);

extern void           * my_calloc(void * const ptr, const size_t nelem, const size_t elsize, const char * const file, const int line);

extern void           * my_malloc(void * const ptr, const size_t nelem, const size_t elsize, const char * const file, const int line);

extern void           * my_realloc(void * const ptr, const size_t nelem, const size_t elsize, const char * const file, const int line);

extern void           open_logs(void);

extern void           time_stamp_logs(void);

extern void           warnf(const char *messg);
extern void           warnvf(const char *fmt, ...);

extern void           dbgvf(const char *fmt, ...);

/*
extern void           input_file_error(void);

*/

#endif
