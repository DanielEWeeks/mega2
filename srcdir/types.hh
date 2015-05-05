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

////////////////////////////////////////////////////////////////////////===
*/

#ifndef TYPES_HH
#define TYPES_HH


#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>

typedef std::string                         Str;
typedef const std::string                   Cstr;

////

typedef std::set<std::string>               Sets;

template <typename T>
bool inline inSet(Cstr& str, T& set) {
    return (set.find(str) != set.cend());
}

////

typedef std::vector<std::string>            Vecs;
typedef std::vector<int>                    Veci;
typedef std::vector<double>                 Vecd;
typedef std::vector<std::string>::const_iterator            Vecsp;
typedef std::vector<int>::const_iterator                    Vecip;
typedef std::vector<double>::const_iterator                 Vecdp;

template <typename T>
class vectordb : public std::vector<T> {
public:
    T& operator[](std::size_t i);
};

template <typename T>
T& vectordb<T>::operator[](std::size_t i) { return (std::vector<T>::data()[i]); }

typedef vectordb<std::string> VecsDB;
typedef vectordb<int>         VeciDB;
typedef vectordb<double>      VecdDB;
typedef vectordb<std::string>::const_iterator            VecspDB;
typedef vectordb<int>::const_iterator                    VecipDB;
typedef vectordb<double>::const_iterator                 VecdpDB;

////

typedef std::list<std::string>              Lists;
typedef std::list<std::string>::const_iterator              Listsp;

template <typename T>
bool inline inList(Cstr& str, T& list) {
    return (list.find(str) != list.cend());
}

////

typedef std::map<const std::string, std::string>  Mapss;
typedef std::map<const std::string, int>          Mapsi;
typedef std::map<int, std::string>          Mapis;
typedef std::map<int, int>                  Mapii;
typedef std::map<const std::string, std::string>::const_iterator  Mapssp;
typedef std::map<const std::string, int>::const_iterator          Mapsip;
typedef std::map<int, std::string>::const_iterator                Mapisp;
typedef std::map<int, int>::const_iterator                        Mapiip;

template <typename T>
bool inline inMap(typename T::key_type& str, T& map) {
    return (map.find(str) != map.cend());
}

template <typename T>
bool map_get(T& map, typename T::key_type& lookup, typename T::mapped_type& ret) {
    typename T::const_iterator ip;
    ip = map.find(lookup);
    if (ip == map.cend())
        return false;
    else {
        ret = ip->second;
        return true;
    }
}

////

typedef std::pair<std::string, std::string> Pairss;
typedef std::pair<std::string, int>         Pairsi;

////

inline const char * C(Cstr& str) {return str.c_str();}
inline int C(int i) {return i;}
inline int C(double d) {return d;}


#endif
