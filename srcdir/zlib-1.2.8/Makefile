# Mega2 version of Makefile
# Makefile for zlib
# Copyright (C) 1995-2013 Jean-loup Gailly, Mark Adler
# For conditions of distribution and use, see copyright notice in zlib.h

# To compile and test, type:
#    ./configure; make test
# Normally configure builds both a static and a shared library.
# If you want to build just a static library, use: ./configure --static

# To use the asm code, type:
#    cp contrib/asm?86/match.S ./match.S
#    make LOC=-DASMV OBJA=match.o

# To install /usr/local/lib/libz.* and /usr/local/include/zlib.h, type:
#    make install
# To install in $HOME instead of /usr/local, use:
#    make install prefix=$HOME

INST = ../libz.a ../zlib.h ../zconf.h

all: ../libz.a

Makefile.Z: Makefile.Z.in
	./configure

../libz.a: libz.a
	cp  $(INST:../%=%) ..
	chmod 644 $(INST)

libz.a: Makefile.Z
	$(MAKE) $(MAKEFLAGS) -f Makefile.Z libz.a

clean: Makefile.Z
	$(MAKE) $(MAKEFLAGS) -f Makefile.Z clean
	$(RM) *.obj
	$(RM) *.lib *.exe *.exe *.pdb
	$(RM) $(INST)

distclean: clean
	$(RM) Makefile.Z
