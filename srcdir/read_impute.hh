/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,
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

#ifndef READ_IMPUTE_HH
#define READ_IMPUTE_HH

class Marker {
public:
    Marker(Cstr name, Cstr chr, Cstr pos, Cstr A, Cstr B, bool read_info): name(name), chr(chr), pos(pos), A(A), B(B) {
        if (! read_info) info = 1.0;
        else info = 0.0;
    }

public:
    Str name;
    Str chr;
    Str pos;
    Str A;
    Str B;
    double info;
    double certainty;
};

typedef vectordb<Marker *> Vecmarkerp;
typedef vectordb<Marker *>::const_iterator Vecmarkerpp;

class ReadImputed {
public:
    ReadImputed(const char *imp, Cstr& info, Cstr& sam) :
        impute_file(imp), info_file(info), sample_file(sam){}
    ~ReadImputed() { (void) markers[0]; }

    void read_imputed_file();

    void read_info_file();

    void read_imputed_file_genotype(const char *imp);

    void read_sample_file();

public:
    Vecmarkerp markers;
    bool       read_info;
    Str        default_chrm;

private:
    const   char *impute_file;
    Str     info_file;
    Str     sample_file;
static
    Str     info_file_hdr;
static
    Str     sample_file_hdr;
    Globals G;
};

#endif
