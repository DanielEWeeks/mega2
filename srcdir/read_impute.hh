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

class ImpMarker {
public:
    ImpMarker(Cstr name, Cstr chr, Cstr pos, Cstr A, Cstr B, bool read_info): name(name), chr(chr), pos(pos), A(A), B(B) {
        if (! read_info)
            info = 1.0;
        else
            info = 0.0;
        skip = false;
    }

public:
    Str name;
    Str chr;
    Str pos;
    Str A;
    Str B;
    double info;
    double certainty;
    bool skip;
};

typedef vectordb<ImpMarker *> Vecmarkerp;
typedef vectordb<ImpMarker *>::const_iterator Vecmarkerpp;

typedef vectordb<Vecs> Vecvecs;
typedef vectordb<Vecs>::const_iterator Vecvecsp;
class Input;

class ReadImputed {
public:
    ReadImputed() {};
    ~ReadImputed() { (void) markers[0]; };

    void files(char *imp, Cstr& info, Cstr& sam) {
        impute_file = imp;
        info_file   = info;
        sample_file = sam;
    }

    void read_imputed_file();

    void read_info_file();

    void read_sample_file();

    void read_genotypes_file();


public:
    Vecmarkerp    markers;
    Vecvecs       people;
    bool          read_info;
    Str           default_chrm;
    double        info_threshold;
    Input_Impute *input;
private:
    const   char *impute_file;
    Str     info_file;
    Str     sample_file;
static
    Str     info_file_hdr;
static
    Str     sample_file_hdr;
    Vecs    sample_file_hdr1a;
    Vecs    sample_file_hdr1b;
    Vecs    sample_file_hdr2a;
    Vecs    sample_file_hdr2b;
};

#endif
