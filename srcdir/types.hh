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

////////////////////////////////////////////////////////////////////////===
*/

#ifndef TYPES_HH
#define TYPES_HH

#include <string.h>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <limits>
#include <map>

typedef std::string                         Str;
typedef const std::string                   Cstr;

////

typedef std::set<std::string>               Sets;

template <typename T>
bool inline inSet(typename T::value_type& str, T& set) {
//  return (set.find(str) != set.cend()); // gcc has issues
    return (set.find(str) != set.end());
}

template <typename T>
bool inline inSet(Cstr& str, T& set) {
//  return (set.find(str) != set.cend()); // gcc has issues
    return (set.find(str) != set.end());
}

////

typedef std::vector<std::string>            Vecs;
typedef std::vector<const char *>           Vecc;
typedef std::vector<int>                    Veci;
typedef std::vector<double>                 Vecd;
typedef std::vector<std::string>::const_iterator            Vecsp;
typedef std::vector<const char *>::const_iterator           Veccp;
typedef std::vector<int>::const_iterator                    Vecip;
typedef std::vector<double>::const_iterator                 Vecdp;

template <class T>
class vectordb : public std::vector<T> {
public:
    T& operator[](unsigned long i);
    const T& operator[](unsigned long i) const;
};

template <class T>
T&      vectordb<T>::operator[](unsigned long i) { return (std::vector<T>::data()[i]); }

template <class T>
const T& vectordb<T>::operator[](unsigned long i) const { return (std::vector<T>::data()[i]); }

typedef vectordb<std::string>  VecsDB;
typedef vectordb<const char *> VeccDB;
typedef vectordb<int>          VeciDB;
typedef vectordb<double>       VecdDB;
typedef vectordb<std::string>::const_iterator            VecspDB;
typedef vectordb<const char *>::const_iterator           VeccpDB;
typedef vectordb<int>::const_iterator                    VecipDB;
typedef vectordb<double>::const_iterator                 VecdpDB;

////

typedef std::list<std::string>              Lists;
typedef std::list<std::string>::const_iterator              Listsp;

template <typename T>
bool inline inList(typename T::value_type& str, T& list) {
//  return (list.find(str) != list.cend()); // gcc has issues
    return (list.find(str) != list.end());
}

////

typedef std::map<const std::string, std::string>  Mapss;
typedef std::map<const std::string, int>          Mapsi;
class charsless {
public:
    bool operator() (const char *a, const char *b) const {
        return strcmp(a, b) < 0;
    }
};
class alleleless {
public:
    bool operator() (const char *a, const char *b) const {
        return a[1] && b[1] ? strcmp(a, b) < 0 : a[0] < b[0];
    }
};
typedef std::map<const char *, const char *, charsless>  Mapcc;
typedef std::map<const char *, int, charsless>           Mapci;
typedef std::map<int, std::string>                       Mapis;
typedef std::map<int, int>                               Mapii;
typedef std::map<const void *, int>                      Mapvi;
typedef std::map<const std::string, std::string>::const_iterator        Mapssp;
typedef std::map<const std::string, int>::const_iterator                Mapsip;
typedef std::map<const char *, const char *, charsless>::const_iterator Mapccp;
typedef std::map<const char *, int, charsless>::const_iterator          Mapcip;
typedef std::map<int, std::string>::const_iterator                      Mapisp;
typedef std::map<int, int>::const_iterator                              Mapiip;
typedef std::map<const void *, int>::const_iterator                     Mapvip;

template <typename T>
bool inline inMap(typename T::key_type& str, T& map) {
//  return (map.find(str) != map.cend()); // gcc has issues
    return (map.find(str) != map.end());
}

template <typename T>
bool map_get(T& map, typename T::key_type& lookup, typename T::mapped_type& ret) {
    typename T::const_iterator ip;
    ip = map.find(lookup);
//  if (ip == map.cend()) // gcc has issues
    if (ip == map.end())
        return false;
    else {
        ret = ip->second;
        return true;
    }
}

////

#ifdef __APPLE__
#include <unordered_map>
#define UM std::unordered_map
#define HH std::hash

#elif defined(_WIN)
#include <unordered_map>
#define UM std::unordered_map
#define HH std::hash

#else
#include <tr1/unordered_map>
#define UM std::tr1::unordered_map
#define HH std::tr1::hash
#endif

class charseq {
public:
    bool operator() (const char *a, const char *b) const {
        return strcmp(a, b) == 0;
    }
}
;
#if defined(__APPLE__) || defined(_WIN)
namespace std {
    template<>
    struct hash<const unsigned char *> {
    public:
        std::size_t operator()(const unsigned char *up) const {
            std::size_t r = 5381;
            int c;
            while ((c = *up++))
                r = ((r << 5) + r) + c;
            return r;
        }
    };
    template<>
    struct hash<const char *> {
    public:
        std::size_t operator()(const char * cp) const {
            std::size_t r = 5381;
            int c;
            while ((c = *cp++))
                r = ((r << 5) + r) + c;
            return r;
        }
    };
}
#else
namespace std {
  namespace tr1 {
    template<>
    struct hash<const unsigned char *> {
    public:
        std::size_t operator()(const unsigned char *up) const {
            std::size_t r = 5381;
            int c;
            while ((c = *up++))
                r = ((r << 5) + r) + c;
            return r;
        }
    };
    template<>
    struct hash<const char *> {
    public:
        std::size_t operator()(const char * cp) const {
            std::size_t r = 5381;
            int c;
            while ((c = *cp++))
                r = ((r << 5) + r) + c;
            return r;
        }
    };
  }
}
#endif

typedef UM<const std::string, std::string, HH<std::string> >        Hmapss;
typedef UM<const std::string, int, HH<std::string> >                Hmapsi;
typedef UM<const char *, const char *, HH<const char *>, charseq>   HMapcc;
typedef UM<const char *, int, HH<const char *>, charseq>            HMapci;
typedef UM<int, std::string, HH<int> >                              HMapis;
typedef UM<int, int, HH<int> >                                      HMapii;
typedef UM<const void *, int, HH<const void *> >                    HMapvi;
typedef UM<const std::string, std::string, HH<std::string> >::const_iterator   HMapssp;
typedef UM<const std::string, int, HH<std::string> >::const_iterator           HMapsip;
typedef UM<const char *, const char *, HH<const char *>, charseq>::const_iterator HMapccp;
typedef UM<const char *, int, HH<const char *>, charseq>::const_iterator       HMapcip;
typedef UM<int, std::string, HH<int> >::const_iterator                         HMapisp;
typedef UM<int, int, HH<int> >::const_iterator                                 HMapiip;
typedef UM<const void *, int, HH<const void *> >::const_iterator               HMapvip;

////

typedef std::pair<std::string, std::string> Pairss;
typedef std::pair<std::string, int>         Pairsi;
typedef std::pair<int, int>                 Pairii;
typedef std::pair<long long, long long>     Pairll;

////

inline const char * C(Cstr& str) {return str.c_str();}
inline int C(int i) {return i;}
inline double C(double d) {return d;}


#endif
