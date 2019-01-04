/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

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

#ifndef STR_UTILS_HH
#define STR_UTILS_HH

#include "types.hh"

/*
  read the next N token from a line
  parse with std::stringstream
*/

class Token {

public: 
    Token(int cnt=0, Cstr white=" \t\n\r"): sep(white), Cline(0), _cnt(cnt){
    }

    bool more(Str& token, int dbg=0);
    bool more(const char *& token);
    void rest(const char *& token) {
        token = &Cline[fo+1];
    }
    void rest(Str& token) {
        token = line.substr(of);
    }

    bool assign(Str& lhs, Str& rhs, int dbg=0);

    void getS(Vecs& vec, int cnt=0);

    void getD(Vecd& vec, int cnt=0);
    void getD(double vec[3], int cnt=0);
    typedef double d3[3];
    bool getDC(d3 &vec, int cnt=0);

    void set(Cstr& line);
    void set(char *line);

    bool   mo;

private:
    Str    sep;
    Str    line;
    char  *Cline;
//  int    npos;
    int    _cnt;
//  int    dbg;
    size_t of, fo;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
/*
  read the next N token from a line
  parse with std::stringstream
*/

#include <sstream>

class TokenSstream {

public: 
    TokenSstream(int cnt=0): _cnt(cnt) {
        ss.clear(std::ios::goodbit);
    }

    /**
      * \brief get a vec of {,_}cnt entries of type of vec from line using std::stringstream conversion
      *
      * \param vec
      * \return void
     */
    template <typename T>
    void get(std::vector<T> & vec, int cnt=0);

    void getS(Vecs& vec, int cnt);

    void getD(Vecd& vec, int cnt);

    void set(Cstr& line);

private:
    std::stringstream ss;
    int _cnt;
};

/**
 * \brief split line into fields using sep
 */
void split(Vecs &fields, Cstr& line, Cstr& sep=" \t\n\r", int cnt=0);
void split(Vecc &fields, char *line, const char *sep=" \t\n\r", int cnt=0);

void join(Vecs &vec, Str& ans, Cstr& sep);

Str& rtrim(Str& str, Cstr& trim=" \t\n\r");

Str& ltrim(Str& str, Cstr& trim=" \t\n\r");

class Globals {
public:
    Globals() {
        chrm_set.insert(chrmSet, chrmSet+31);
        chrm_renum.insert(chrmRenum, chrmRenum+31);
        chrm_int.insert(chrmInt, chrmInt+31);
    }

    Sets chrm_set;
    Mapss chrm_renum;
    Mapsi chrm_int;

private:
    static const Str chrmSet[31];
    static const Pairss chrmRenum[31];
    static const Pairsi chrmInt[31];
};

#endif
