/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#ifndef CW_ROUTINES_EXT_H
#define CW_ROUTINES_EXT_H


#if (defined(SOLARIS) || defined(MINGW) || defined(_WIN))
extern char *strsep(char **stringp, const char *delim);
#endif


#if defined(__MWERKS__) || defined(SOLARIS) || defined(MINGW)
extern void cfree(void *ptr);
#endif

/*

#ifdef __MWERKS__
extern int access(const char*path, int amode);

extern int safe_tolower(int c);

extern int strcasecmp(char *s1, char *s2);

extern char *strdup(char *s1);

extern int strncasecmp(char *s1, char *s2, size_t n);

#endif

*/

#endif
