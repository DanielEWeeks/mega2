#   Mega2: Manipulation Environment for Genetic Analysis
#  
#   Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.
#  
#   Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
#   Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
#   and Daniel E. Weeks.

#   This file is part of the Mega2 program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation;
#   either version 3 of the License, or (at your option) any later
#   version.
#  
#   Mega2 is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#  
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#  
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
# 
# ===========================================================================
# Makefile for BGEN as a library for Mega2

# Put the library in the Mega2 srcdir above this...
# NOTE: -lz is needed when the library is linked...
ifeq ($(OSTYPE),windows)
UP = ..\\..\\
WINDIR=\\
#WINO = -Fo$(@D)/$*.o
WINO = -Fo$*.o
else
UP = ../../
WINDIR=/
WINO = -o $@
endif

DEFS = -D_mega2_ -DPACKAGE_NAME=\"sqlite\" -DPACKAGE_TARNAME=\"sqlite\" -DPACKAGE_VERSION=\"3.9.2\" -DPACKAGE_STRING=\"sqlite\ 3.9.2\" -DPACKAGE_BUGREPORT=\"http://www.sqlite.org\" -DPACKAGE_URL=\"\" -DPACKAGE=\"sqlite\" -DVERSION=\"3.9.2\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_DLFCN_H=1 -DLT_OBJDIR=\".libs/\" -DHAVE_FDATASYNC=1 -DHAVE_USLEEP=1 -DHAVE_GMTIME_R=1 -DHAVE_DECL_STRERROR_R=1 -DHAVE_STRERROR_R=1
# -DHAVE_READLINE=1
ifeq ($(OSTYPE),darwin)
DEFS += -DSQLITE_WITHOUT_ZONEMALLOC
endif

ifeq ($(SQLxDEFS),)
SQLxDEFS = -DHAVE_INTTYPES_H=1 -DHAVE_LOCALTIME_R=1
endif

THREADSAFE_FLAGS = -D_REENTRANT=1 -DSQLITE_THREADSAFE=1

AM_CFLAGS = -D_REENTRANT=1 -DSQLITE_THREADSAFE=1    -DSQLITE_ENABLE_FTS3 -DSQLITE_ENABLE_RTREE

OUT = $(UP)libsqlite3.a

#SRC = example/bgen_to_vcf.cpp src/MissingValue.cpp src/bgen.cpp src/zlib.cpp
BSRC = sqlite3.c
BOBJ = $(BSRC:.c=.o)
WOBJ = $(subst /,$(WINDIR),$(BOBJ))

# include directories '.' for other vcftools headers, '..' for Mega2 headers... and zlib headers
LOCFLAGS += -I.. $(DEFS) $(SQLxDEFS) $(THREADSAFE_FLAGS) $(AM_CFLAGS)

# CC, CFLAGS are defined in Mega2 Makedefs

.SUFFIXES: .c
.c.o:
	$(CC) $(INCLUDES) $(CFLAGS) $(LOCFLAGS) $(FLAGS) -c $< $(WINO)

all: $(OUT)


more: $(OUT) sqlite3

sqlite3: shell.o $(OUT)
	$(CC) -o sqlite3 $< $(OUT) 

sqlite3.o: sqlite3.c
#sqlite3.o: CFLAGS += $(L0xFLAGS)


$(OUT): $(BOBJ)
	$(AR) $(ARFLAGS) $(AROUT)$(OUT) $(BOBJ)

clean:
	$(RM) $(WOBJ) $(OUT) Makefile.bak
	$(RM) sqlite3 shell.o

distclean: clean

depend:
	makedepend -Y. -- $(CFLAGS) -- $(INCLUDES) $(BSRC)
# DO NOT DELETE
