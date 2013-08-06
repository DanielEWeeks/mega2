/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2013 Robert Baron, Charles P. Kollar,
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

#ifndef LOOP_H
#define LOOP_H

#include "file_ops.h"
#include "entry.h"
#include "sh_util.h"

#include "linkage_ext.h"

class loop {
////////////////////////////////////////////////////////////////
//                       Outer loops
////////////////////////////////////////////////////////////////
public:
    class once: public virtual entry {
    public:
        int   *_trp;
        bool _trait_affect;
        bool _trait_quant;

        once(linkage_ped_top *Top) : entry(Top) {
            _numchr     = 0;
            _chrom_loop = LoopOverChrm && main_chromocnt > 1;

            _trp          = global_trait_entries;
            _trait_affect = false;
            _trait_quant  = false;
            _tte          = 0;

            _opath      = &output_paths[0];
        }
        virtual ~once() {}

        void iterate();
        virtual void make_file() {}
    };

public:
    class outer: public virtual entry {
    public:
        int _trait_loop;
        int   *_trp;
        bool _trait_affect;
        bool _trait_quant;

        outer(linkage_ped_top *Top) : entry(Top) {
            _numchr       = 0;
            _chrom_loop   = LoopOverChrm && main_chromocnt > 1;

            _trait_loop   = LoopOverTrait && num_traits > 1;
            _trp          = global_trait_entries;
            _trait_affect = false;
            _trait_quant  = false;
            _tte          = 0;

            _opath        = _trait_loop ? &output_paths[1] : &output_paths[0];
        }
        virtual ~outer() {}

        void iterate();
        virtual void chr_start() {}
        virtual void trait_start() {}
        virtual void make_file() {}
        virtual void trait_end() {}
        virtual void chr_end() {}
    };

public:
    class chr: public virtual entry {
    public:
        int   *_trp;
        bool _trait_affect;
        bool _trait_quant;

        chr(linkage_ped_top *Top) : entry(Top) {
            _numchr     = 0;
            _chrom_loop = LoopOverChrm && main_chromocnt > 1;

            _trp          = global_trait_entries;
            _trait_affect = false;
            _trait_quant  = false;
            _tte          = 0;

            _opath      = &output_paths[0];
        }
        virtual ~chr() {}

        void iterate();
        virtual void chr_start() {}
        virtual void make_file() {}
        virtual void chr_end() {}
    };

public:
    class trait: public virtual entry {
    public:
        int _trait_loop;
        int   *_trp;
        bool _trait_affect;
        bool _trait_quant;

        trait(linkage_ped_top *Top) : entry(Top) {
            _numchr       = 0;
            _chrom_loop   = LoopOverChrm && main_chromocnt > 1;

            _trait_loop   = LoopOverTrait && num_traits > 1;
            _trp          = global_trait_entries;
            _tte          = 0;
            _trait_affect = false;
            _trait_quant  = false;

            _opath        = _trait_loop ? &output_paths[1] : &output_paths[0];
        }
        virtual ~trait() {}

        void iterate();
        virtual void trait_start() {}
        virtual void make_file() {}
        virtual void trait_end() {}
    };

////////////////////////////////////////////////////////////////
//                       inner loops
////////////////////////////////////////////////////////////////

public:
    class null: public virtual entry {
    public:

        null(linkage_ped_top *Top) : entry(Top) { }
       ~null() {}

        void run_loop(const char *dir, const char *fl_name, const char *mode="w");
        void run_loop(const char *fl_name) { run_loop(*_opath, fl_name, "w"); }
        virtual void inner() {}
    };

public:
    class ped_per: public virtual entry {
    public:
        ped_per(linkage_ped_top *Top) : entry(Top) { }
       ~ped_per() {}

        void run_loop(const char *dir, const char *fl_name, const char *mode="w");
        void run_loop(const char *fl_name) { run_loop(*_opath, fl_name, "w"); }
        virtual void ped_start() {};
        virtual void inner() {};
        virtual void ped_end() {};
    };

public:
    class ped_per_trait: public virtual entry {
        int *trp;
        bool trait_affect;
        bool trait_quant;
    public:
        ped_per_trait(linkage_ped_top *Top) : entry(Top) {
            trp          = global_trait_entries;
            trait_affect = false;
            trait_quant  = false;
        }

       ~ped_per_trait() {}

        void run_loop(const char *dir, const char *fl_name, const char *mode="w");
        void run_loop(const char *fl_name) { run_loop(*_opath, fl_name, "w"); }
        virtual void ped_start() {};
        virtual void per_start() {};
        virtual void inner() {};
        virtual void per_end() {};
        virtual void ped_end() {};
    };

public:
    class trait_ped_per: public virtual entry {
        int *trp;
        bool trait_affect;
        bool trait_quant;
    public:
        trait_ped_per(linkage_ped_top *Top) : entry(Top) {
            trp          = global_trait_entries;
            trait_affect = false;
            trait_quant  = false;
        }

       ~trait_ped_per() {}

        void run_loop(const char *dir, const char *fl_name, const char *mode="w");
        void run_loop(const char *fl_name) { run_loop(*_opath, fl_name, "w"); }
        virtual void trait_start() {};
        virtual void ped_start() {};
        virtual void inner() {};
        virtual void ped_end() {};
        virtual void trait_end() {};
    };

public:
    class ped_per_loci: public virtual entry {
    public:
        int _loci_allele_limit;

        ped_per_loci(linkage_ped_top *Top) : entry(Top) {
            _loci_allele_limit = 2;
        }
       ~ped_per_loci() {}

        void run_loop(const char *dir, const char *fl_name, const char *mode="w");
        void run_loop(const char *fl_name) { run_loop(*_opath, fl_name, "w"); }
        virtual void ped_start() {}
        virtual void per_start() {}
        virtual void inner() {}
        virtual void per_end() {}
        virtual void ped_end() {}
    };

public:
    class loci_ped_per: public virtual entry {
    public:
        int _loci_allele_limit;

        loci_ped_per(linkage_ped_top *Top) : entry(Top) {
            _loci_allele_limit = 2;
        }
       ~loci_ped_per() {}

        void run_loop(const char *dir, const char *fl_name, const char *mode="w");
        void run_loop(const char *fl_name) { run_loop(*_opath, fl_name, "w"); }
        virtual void loci_start() {}
        virtual void ped_start() {}
        virtual void inner() {}
        virtual void ped_end() {}
        virtual void loci_end() {}
    };

public:
    class loci: public virtual entry {
    public:
        int _loci_allele_limit;

        loci(linkage_ped_top *Top) : entry(Top) {
            _loci_allele_limit = 2;
        }
       ~loci() {}
        void delete_file(const char *file) { filep_delete(*_opath, file); }

        void run_loop(const char *dir, const char *fl_name, const char *mode="w");
        void run_loop(const char *fl_name) { run_loop(*_opath, fl_name, "w"); }
        virtual void inner() {}
    };

};
#endif
