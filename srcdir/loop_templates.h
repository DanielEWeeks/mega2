/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2015 Robert Baron, Charles P. Kollar,
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

#ifndef LOOP_TEMPLATES_H
#define LOOP_TEMPLATES_H


struct FLOOPonce: public fileloop::once {
    const char *file_name;
    const char *file_type;
    const char *file_mode;

    FLOOPonce(linkage_ped_top *Top) : fileloop::once(Top), file_mode("w") { }
    FLOOPonce(linkage_ped_top *Top, const char *f_name, const char *f_mode) :
        fileloop::once(Top), file_name(f_name), file_mode(f_mode) { }

    void make_file() {
        if (*file_type) msgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

struct FLOOPboth: public fileloop::both {
    const char *file_name;
    const char *file_type;
    const char *file_mode;

    FLOOPboth(linkage_ped_top *Top) : fileloop::both(Top), file_mode("w") { }
    FLOOPboth(linkage_ped_top *Top, const char *f_name, const char *f_mode) :
        fileloop::both(Top), file_name(f_name), file_mode(f_mode) { }

    void make_file() {
        if (*file_type) mssgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

struct FLOOPchr: public fileloop::chr {
    const char *file_name;
    const char *file_type;
    const char *file_mode;

    FLOOPchr(linkage_ped_top *Top) : fileloop::chr(Top), file_mode("w") { }
    FLOOPchr(linkage_ped_top *Top, const char *f_name, const char *f_mode) :
        fileloop::chr(Top), file_name(f_name), file_mode(f_mode) { }

    void make_file() {
        if (*file_type) mssgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

struct FLOOPtrait: public fileloop::trait {
    const char *file_name;
    const char *file_type;
    const char *file_mode;

    FLOOPtrait(linkage_ped_top *Top) : fileloop::trait(Top), file_mode("w") { }
    FLOOPtrait(linkage_ped_top *Top, const char *f_name, const char *f_mode) :
        fileloop::trait(Top), file_name(f_name), file_mode(f_mode) { }

    void make_file() {
        if (*file_type) mssgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

struct DTshell: public sh_util {
    DTshell(linkage_ped_top *Top) : sh_util(Top) {}
    DTshell(linkage_ped_top *Top, fileloop::fileloop_data *fl) : sh_util(Top, fl) { }
    virtual ~DTshell() {}

    virtual void file_post() {
        chmod_X_file(path_);
    }
};

#endif
