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

go: all

ifeq      ($(OSTYPE),windows) #### WINDOWS ####


$(TARGET)::
	$(MAKE) $(MFLAGS) -f win32/Makefile.msc $(TARGET)

TARGET = zlib.lib
CP = copy

clean::
	$(MAKE) $(MFLAGS) -f win32/Makefile.msc clean

else ifeq ($(OSTYPE),mingw) #### MINGW ####


TARGET = libz.a
CP = cp

$(TARGET)::
	$(MAKE) $(MFLAGS) -f win32/Makefile.gcc $(TARGET)

clean::
	$(MAKE) $(MFLAGS) -f win32/Makefile.gcc clean

else                       #### DEFAULT ####

TARGET = libz.a
CP = cp

Makefile.zlib: Makefile.zlib.in
	./configure --static

$(TARGET):: Makefile.zlib
	$(MAKE) $(MFLAGS) -f Makefile.zlib $(TARGET)

clean::
ifneq ($(wildcard Makefile.zlib),)
	$(MAKE) $(MFLAGS) -f Makefile.zlib clean
endif
endif

INST = ../$(TARGET) ../zlib.h ../zconf.h

all: ../$(TARGET)

../$(TARGET): $(TARGET)
	$(foreach fl,$(INST:../%=%),$(CP) $(fl) ..;)
	chmod 644 $(INST)

clean::
	$(RM) *.obj
	$(RM) *.lib *.exe *.exe *.pdb
	$(RM) $(INST)

distclean: clean
	$(RM) Makefile.zlib
	$(CP) zconf.h.in zconf.h
