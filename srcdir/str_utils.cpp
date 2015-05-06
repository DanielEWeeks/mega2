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

#include <stdlib.h>
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

void Token::set(Cstr& line) {
    this->line = line;
/*
    npos = this->line.size();
    if (this->line[npos] == '\n') {
        npos--;
        this->line.resize(npos);
    }
*/
    of = line.find_first_not_of(sep);
    if (of == std::string::npos) of = 0;
    fo = std::string::npos;
    mo = true;
}

void Token::set(const char * line) {
    set(string(line));
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
    size_t  fo = 0;
    bool mo = true;
    size_t  of = line.find_first_not_of(sep, fo);
    if (of == std::string::npos) {
        fields.push_back(line);
        return;
    }
    int  i  = 1;
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
            of = line.find_first_not_of(sep, fo);
            if (cnt && i++ >= cnt) {
                fields.push_back(line.substr(of, line.size()-of+1));
                return;
            }
            if (of == std::string::npos) return;
        }
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
