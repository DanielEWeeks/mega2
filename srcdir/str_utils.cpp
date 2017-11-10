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
#include <stdlib.h>
#include <string.h>
#include "common.h"
#include "error_messages_ext.h"

#include <string>
#include <vector>
#include "str_utils.hh"

/*
  read the next N token from a line
  parse with stringstream
*/

using namespace std;

bool Token::more(Str& token, int dbg) {
    fo = line.find_first_of(sep, of);
//  if (dbg) cout << "of: " << of << "; fo: " << fo << endl;
    if (fo == std::string::npos) {
        fo = line.size();
//      if (dbg) cout << "of: " << of << "; fo: " << fo << endl;
//      if (dbg)cout << line.substr(of, fo-of) << endl;
        token = line.substr(of, fo-of);
        return false;
    } else {
//      if (dbg) cout << line.substr(of, fo-of) << endl;
        token = line.substr(of, fo-of);
        of = line.find_first_not_of(sep, fo);
        return true;
    }
}

bool Token::more(const char *&token) {
    fo = line.find_first_of(sep, of);
    if (fo == std::string::npos) {
        fo = line.size();
        token = &Cline[of];
        Cline[fo] = 0;
        return false;
    } else {
        token = &Cline[of];
        Cline[fo] = 0;

        of = line.find_first_not_of(sep, fo);
        if (of == std::string::npos)
            return false;
        else
            return true;
    }
}

//#include <iostream>
bool Token::assign(Str& lhs, Str& rhs, int dbg) {

    size_t mid = line.find_first_of("=", 0);
    if (mid == std::string::npos) {
        lhs = rhs = "";
        return false;
    }
    mid--;
    fo = line.find_last_not_of(sep, mid);
 
//  if (dbg) cout << "of: " << of << "; fo: " << fo << endl;
    if (fo == std::string::npos) {
        fo = mid;
//      if (dbg) cout << "of: " << of << "; fo: " << fo << endl;
    }
//  if (dbg) cout << line.substr(of, fo-of+1) << endl;
    lhs = line.substr(of, fo-of+1);

    mid += 2;
    of = line.find_first_not_of(sep, mid);
    if (of == std::string::npos)
        of = mid;
    fo = line.find_last_not_of(sep);
 
//  if (dbg) cout << "of: " << of << "; fo: " << fo << endl;
    if (fo < mid) {
        fo = line.size();
//      if (dbg) cout << "of: " << of << "; fo: " << fo << endl;
    }
//  if (dbg) cout << line.substr(of, fo-of+1) << endl;
    rhs = line.substr(of, fo-of+1);

    return true;
}

void Token::getS(Vecs& vec, int cnt) {
    string token;
    vec.clear();
    for (int i = 0, l = cnt ? cnt : _cnt; mo && i < l; i++) {
        mo = more(token);
        vec.push_back(token);
    }
}

void Token::getD(Vecd& vec, int cnt) {
    Str token;
    double dbl;
    char *fin;
    vec.clear();
    for (int i = 0, l = cnt ? cnt : _cnt; mo && i < l; i++) {
        mo = more(token);
        dbl = strtod(token.c_str(), &fin); //4.8
        if (*fin != 0) {
            dbl = 0.0;
            warnvf("Invalid floating point number: %s\n", token.c_str());
        }
        vec.push_back(dbl);
    }
}

void Token::getD(double vec[3], int cnt) {
    Str token;
    double dbl;
    char *fin;
    for (int i = 0, l = cnt ? cnt : _cnt; mo && i < l; i++) {
        mo = more(token);
        dbl = strtod(token.c_str(), &fin); //4.8
        if (*fin != 0) {
            dbl = 0.0;
            warnvf("Invalid floating point number: %s\n", C(token));
        }
        vec[i] = dbl;
    }
}

bool Token::getDC(d3 &vec, int cnt) {
    const char *token;
    char *fin;
    int i, l;
    for (i = 0, l = cnt ? cnt : _cnt; mo && i < l; i++) {
        mo = more(token);
        if (token[1] == 0) {
            if (*token == '1') {
                vec[i] = 1.0;
                continue;
            } else if (*token == '0') {
                vec[i] = 0.0;
                continue;
            }
        }
        vec[i] = strtod(token, &fin);
        if (*fin != 0) {
            vec[i] = 0.0;
            warnvf("Invalid floating point number: %s\n", token);
        }
    }
    return i == 3;
}

void Token::set(Cstr& lin) {
    this->line = lin;

    if (Cline) delete[] Cline;
    int sz = lin.size();
    Cline = new char[sz+1]; // to be safe

    lin.copy(Cline, sz, 0);
    Cline[sz] = 0;

    of = line.find_first_not_of(sep);
    if (of == std::string::npos) of = 0;
    fo = std::string::npos;
    mo = true;
}

void Token::set(char *lin) {
    this->line = string(lin);

    Cline = lin;

    of = line.find_first_not_of(sep);
    if (of == std::string::npos) of = 0;
    fo = std::string::npos;
    mo = true;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
/*
  read the next N token from a line
  parse with stringstream
*/

#include <sstream>

/**
  * \brief get a vec of {,_}cnt entries of type of vec from line using stringstream conversion
  *
  * \param vec
  * \return void
 */

template <typename T>
void TokenSstream::get(vector<T> & vec, int cnt) {
    ss.clear(ios::goodbit);
    vec.clear();
    T token;
    for (int i = 0, l = cnt ? cnt : _cnt; i < l; i++) {
        ss >> token;
        if (ss.eof() || ss.fail()) return;
        vec.push_back(token);
    }
}

void TokenSstream::getS(Vecs& vec, int cnt) {
    ss.clear(ios::goodbit);
    vec.clear();
    Str token;
    for (int i = 0, l = cnt ? cnt : _cnt; i < l; i++) {
        ss >> token;
        if (ss.eof() || ss.fail()) return;
        vec.push_back(token);
    }
}

void TokenSstream::getD(Vecd& vec, int cnt) {
    ss.clear(ios::goodbit);
    vec.clear();
    Str token;
    double dbl;
    for (int i = 0, l = cnt ? cnt : _cnt; i < l; i++) {
        ss >> token;
        if (ss.eof() || ss.fail()) return;
        dbl = atof(token.c_str());
        vec.push_back(dbl);
    }
}

void TokenSstream::set(Cstr& line) {
    ss.str(line);
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

/**
 * \brief split line into fields using sep
 */

void split(Vecs &fields, Cstr& line, Cstr& sep, int cnt) {
    int mult = (sep.find_first_of(" ") != std::string::npos) ? 1 : 0;
    size_t  fo = 0;
    size_t  of;
    if (mult) {
	of = line.find_first_not_of(sep, fo);
	if (of == std::string::npos) {
	    fields.push_back(line);
	    return;
	}
    } else
	of = fo;

    bool mo = true;
    int i = 1;
    while (mo) {
        fo = line.find_first_of(sep, of);
//      if (dbg) cout << "of: " << of << "; fo: " << fo << endl;
        if (fo == std::string::npos) {
            fo = line.size();
//          if (dbg) cout << "of: " << of << "; fo: " << fo << endl;
//          if (dbg)cout << line.substr(of, fo-of) << endl;
            fields.push_back(line.substr(of, fo-of));
            return;
        } else {
//          if (dbg) cout << line.substr(of, fo-of) << endl;
            fields.push_back(line.substr(of, fo-of));
            if (mult == 1) {
                of = line.find_first_not_of(sep, fo);
		if (of == std::string::npos) {
//		    fields.push_back("");
		    return;
		}
            } else
                of = fo + 1;

            if (cnt && i++ >= cnt) {
                fields.push_back(line.substr(of, line.size()-of+1));
                return;
            }
        }
    }
}

#if defined(MINGW) || defined(_WIN) || defined(SOLARIS)
static const char *index(const char *s, char c) {
    const char *sp;
    for (sp = s; *sp != 0 && *sp != c; sp++) ;

    return *sp ? sp : (const char *) 0;
}
#endif

void split(Vecc &fields, char *line, const char *sep, int cnt) {
    int mult = index(sep, ' ') != 0 ? 1 : 0;
    char *fo = line;
    char *of = fo;
    if (mult) {
	while ( *of && index(sep, *of) ) of++;
	if (of == 0) {
	    fields.push_back(line);
	    return;
	}
    }

    bool mo = true;
    int i = 1;
    while (mo) {
        fo = strpbrk(of, sep);
        if (fo == 0) {
            fields.push_back(of);
            return;
        } else {
            *fo = 0;
            fields.push_back(of);
	    of = fo + 1;
            if (mult == 1) {
		while ( *of && index(sep, *of) ) of++;
		if (of == 0) return;
	    }
            if (cnt && i++ >= cnt) {
                fields.push_back(of);
                return;
            }
        }
    }
}

void join(Vecs &vec, Str& ans, Cstr& sep) {
    ans = "";
    for (int i = 0, l = vec.size(); i < l; i++) {
        ans.append(vec[i]);
        if (i < l - 1)
            ans.append(sep);
    }
}

Str& rtrim(Str& str, Cstr& trim=" \t\n\r") {
    return str.erase(str.find_last_not_of(trim)+1);
}

Str& ltrim(Str& str, Cstr& trim=" \t\n\r") {
    return str.erase(0, str.find_first_not_of(trim));
}

const Str Globals::chrmSet[31] = {  "0",
         "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "23", "24", "25", "26",  "X",  "Y", "XY", "MT" };


const Pairss Globals::chrmRenum[31] = {
    make_pair("0", "0"),   make_pair("1", "1"),   make_pair("2", "2"),   make_pair("3", "3"),
    make_pair("4", "4"),   make_pair("5", "5"),   make_pair("6", "6"),   make_pair("7", "7"),
    make_pair("8", "8"),   make_pair("9", "9"),   make_pair("10", "10"), make_pair("11", "11"),
    make_pair("12", "12"), make_pair("13", "13"), make_pair("14", "14"), make_pair("15", "15"),
    make_pair("16", "16"), make_pair("17", "17"), make_pair("18", "18"), make_pair("19", "19"),
    make_pair("20", "20"), make_pair("21", "21"), make_pair("22", "22"), make_pair("23", "23"),
    make_pair("24", "24"), make_pair("25", "25"), make_pair("26", "26"), make_pair("X", "23"),
    make_pair("Y", "24"),  make_pair("XY", "25"), make_pair("MT", "26") };


const Pairsi Globals::chrmInt[31] = {
    make_pair("0", 0),   make_pair("1", 1),   make_pair("2", 2),   make_pair("3", 3),
    make_pair("4", 4),   make_pair("5", 5),   make_pair("6", 6),   make_pair("7", 7),
    make_pair("8", 8),   make_pair("9", 9),   make_pair("10", 10), make_pair("11", 11),
    make_pair("12", 12), make_pair("13", 13), make_pair("14", 14), make_pair("15", 15),
    make_pair("16", 16), make_pair("17", 17), make_pair("18", 18), make_pair("19", 19),
    make_pair("20", 20), make_pair("21", 21), make_pair("22", 22), make_pair("23", 23),
    make_pair("24", 24), make_pair("25", 25), make_pair("26", 26), make_pair("X", 23),
    make_pair("Y",  24), make_pair("XY", 25), make_pair("MT", 26) };
