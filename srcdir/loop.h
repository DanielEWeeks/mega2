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

#ifndef LOOP_H
#define LOOP_H

#include "file_ops.h"
#include "entry.h"

#include "linkage_ext.h"

#define vlpCLASS(derived,outer,inner) struct derived: public fileloop::outer, public dataloop::inner
#define  vlpCTOR(derived,outer,inner) derived(linkage_ped_top *Top) : person_locus_entry(Top), fileloop::outer(Top), dataloop::inner(Top)

#define lpCLASS(derived,base) struct derived: public dataloop::base
#define  lpCTOR(derived,base) derived(linkage_ped_top *Top, fileloop::fileloop_data *fl) : person_locus_entry(Top), dataloop::base(Top, fl)

namespace dataloop {
    class dataloop_data;
};

namespace fileloop {
////////////////////////////////////////////////////////////////
//                       Outer loops
////////////////////////////////////////////////////////////////
    class fileloop_data: public virtual person_locus_entry {
    public:
        fileloop_data(linkage_ped_top *Top) : person_locus_entry(Top), _dataloop(0) { }
        int    _trait_loop;
        int   *_trp;
        bool   _trait_affect;
        bool   _trait_quant;
        int    _chrom_loop;
        char **_opath;

/*
        linkage_locus_rec *_tte; // shared to dataloop
        int _numchr;             // shared to dataloop
        int _trait;              // shared to dataloop

        linkage_ped_top  *_Top;   // copy same as dataloop
        linkage_locus_top  *_LTop;// copy same as dataloop
*/

        dataloop::dataloop_data *_dataloop; // who called
    };

    class once: public fileloop_data {
    public:
        once(linkage_ped_top *Top) : person_locus_entry(Top), fileloop_data(Top) {
            _ttraitp     = 0;
            _numchr      = 0;

            _chrom_loop   = LoopOverChrm && main_chromocnt > 1;

            _trp          = global_trait_entries;
            _trait_loop   = false;
            _trait_affect = false;
            _trait_quant  = false;

            _opath      = &output_paths[0];
        }
        virtual ~once() {}

        void iterate();
        virtual void file_loop() {}
    };

    class both: public fileloop_data {
    public:
        both(linkage_ped_top *Top) : person_locus_entry(Top), fileloop_data(Top) {
            _ttraitp     = 0;
            _numchr      = 0;

            _chrom_loop   = LoopOverChrm && main_chromocnt > 1;

            _trait_loop   = LoopOverTrait && num_traits > 1;
            _trp          = global_trait_entries;
            _trait_affect = false;
            _trait_quant  = false;

            _opath        = _trait_loop ? &output_paths[1] : &output_paths[0];
        }
        virtual ~both() {}

        void iterate();
        virtual void chr_start() {}
        virtual void trait_start() {}
        virtual void file_loop() {}
        virtual void trait_end() {}
        virtual void chr_end() {}
    };

    class chr: public fileloop_data {
    public:
        chr(linkage_ped_top *Top) : person_locus_entry(Top), fileloop_data(Top) {
            _ttraitp     = 0;
            _numchr      = 0;

            _chrom_loop   = LoopOverChrm && main_chromocnt > 1;

            _trp          = global_trait_entries;
            _trait_loop   = false;
            _trait_affect = false;
            _trait_quant  = false;

            _opath      = &output_paths[0];
        }
        virtual ~chr() {}

        void iterate();
        virtual void chr_start() {}
        virtual void file_loop() {}
        virtual void chr_end() {}
    };

    class trait: public fileloop_data {
    public:
        trait(linkage_ped_top *Top) : person_locus_entry(Top), fileloop_data(Top) {
            _ttraitp     = 0;
            _numchr      = 0;

            _chrom_loop   = LoopOverChrm && main_chromocnt > 1;

            _trait_loop   = LoopOverTrait && num_traits > 1;
            _trp          = global_trait_entries;
            _trait_affect = false;
            _trait_quant  = false;

            _opath        = _trait_loop ? &output_paths[1] : &output_paths[0];
        }
        virtual ~trait() {}

        void iterate();
        virtual void trait_start() {}
        virtual void file_loop() {}
        virtual void trait_end() {}
    };
};

namespace dataloop {
////////////////////////////////////////////////////////////////
//                       inner loops
////////////////////////////////////////////////////////////////

    class dataloop_data: public virtual person_locus_entry {
    public:
        fileloop::fileloop_data *_fileloop;

        dataloop_data(linkage_ped_top *Top) : person_locus_entry(Top), _fileloop(0) { }
        dataloop_data(linkage_ped_top *Top, fileloop::fileloop_data *fl) : person_locus_entry(Top) {
            _fileloop = fl;
            fl->_dataloop = this;
        }

        virtual void data_loop(const char *dir, const char *fl_name, const char *mode) {};

        inline void data_init() {
            if (_fileloop) {
                this->_trait  = _fileloop->_trait;
                this->_numchr = _fileloop->_numchr;

                this->_ttraitp = _fileloop->_ttraitp;
            }
        }

    };

    class null: public dataloop_data {
    public:

        null(linkage_ped_top *Top) : person_locus_entry(Top), dataloop_data(Top) { }
        null(linkage_ped_top *Top, fileloop::fileloop_data *fl) : person_locus_entry(Top), dataloop_data(Top, fl) { }
       ~null() {}

        void data_loop(const char *dir, const char *fl_name, const char *mode);
        virtual void inner() {}
    };

    class ped_per: public dataloop_data {
    public:
        ped_per(linkage_ped_top *Top) : person_locus_entry(Top), dataloop_data(Top) { }
        ped_per(linkage_ped_top *Top, fileloop::fileloop_data *fl) : person_locus_entry(Top), dataloop_data(Top, fl) { }
       ~ped_per() {}

        void data_loop(const char *dir, const char *fl_name, const char *mode);
        virtual void ped_start() {};
        virtual void inner() {};
        virtual void ped_end() {};
    };

    class ped_per_trait: public dataloop_data {
        int *trp;
        bool trait_affect;
        bool trait_quant;
    public:
        ped_per_trait(linkage_ped_top *Top, fileloop::fileloop_data *fl) : person_locus_entry(Top), dataloop_data(Top, fl) {
            trp          = global_trait_entries;
            trait_affect = false;
            trait_quant  = false;
        }
        ped_per_trait(linkage_ped_top *Top) : person_locus_entry(Top), dataloop_data(Top) {
            trp          = global_trait_entries;
            trait_affect = false;
            trait_quant  = false;
        }
       ~ped_per_trait() {}

        void data_loop(const char *dir, const char *fl_name, const char *mode);
        virtual void ped_start() {};
        virtual void per_start() {};
        virtual void inner() {};
        virtual void per_end() {};
        virtual void ped_end() {};
    };

    class trait_ped_per: public dataloop_data {
        int *trp;
        bool trait_affect;
        bool trait_quant;
    public:
        trait_ped_per(linkage_ped_top *Top, fileloop::fileloop_data *fl) : person_locus_entry(Top), dataloop_data(Top, fl) {
            trp          = global_trait_entries;
            trait_affect = false;
            trait_quant  = false;
        }
        trait_ped_per(linkage_ped_top *Top) : person_locus_entry(Top), dataloop_data(Top) {
            trp          = global_trait_entries;
            trait_affect = false;
            trait_quant  = false;
        }
       ~trait_ped_per() {}

        void data_loop(const char *dir, const char *fl_name, const char *mode);
        virtual void trait_start() {};
        virtual void ped_start() {};
        virtual void inner() {};
        virtual void ped_end() {};
        virtual void trait_end() {};
    };

    class ped_per_loci: public dataloop_data {
    public:
        int _loci_allele_limit;

        ped_per_loci(linkage_ped_top *Top, fileloop::fileloop_data *fl) : person_locus_entry(Top), dataloop_data(Top, fl) {
            _loci_allele_limit = 0;
        }
        ped_per_loci(linkage_ped_top *Top) : person_locus_entry(Top), dataloop_data(Top) {
            _loci_allele_limit = 0;
        }
       ~ped_per_loci() {}

        void data_loop(const char *dir, const char *fl_name, const char *mode);
        virtual void ped_start() {}
        virtual void per_start() {}
        virtual void inner() {}
        virtual void per_end() {}
        virtual void ped_end() {}
    };

    class loci_ped_per: public dataloop_data {
    public:
        int _loci_allele_limit;

        loci_ped_per(linkage_ped_top *Top, fileloop::fileloop_data *fl) : person_locus_entry(Top), dataloop_data(Top, fl) {
            _loci_allele_limit = 0;
        }
        loci_ped_per(linkage_ped_top *Top) : person_locus_entry(Top), dataloop_data(Top) {
            _loci_allele_limit = 0;
        }
       ~loci_ped_per() {}

        void data_loop(const char *dir, const char *fl_name, const char *mode);
        virtual void loci_start() {}
        virtual void ped_start() {}
        virtual void inner() {}
        virtual void ped_end() {}
        virtual void loci_end() {}
    };

    class loci: public dataloop_data {
    public:
        int _loci_allele_limit;

        loci(linkage_ped_top *Top, fileloop::fileloop_data *fl) : person_locus_entry(Top), dataloop_data(Top, fl) {
            _loci_allele_limit = 0;
        }
        loci(linkage_ped_top *Top) : person_locus_entry(Top), dataloop_data(Top) {
            _loci_allele_limit = 0;
        }
       ~loci() {}

        void delete_file(const char *opath, const char *file) { filep_delete(opath, file); }

        void data_loop(const char *dir, const char *fl_name, const char *mode);
        virtual void inner() {}
    };

};

////////////////////////////////////////////////////////////////
//                   special outer loops
////////////////////////////////////////////////////////////////

struct FLOOPonce: public fileloop::once {
    const char *file_name;
    const char *file_type;
    const char *file_mode;

    FLOOPonce(linkage_ped_top *Top, const char *f_name, const char *f_mode) :
        person_locus_entry(Top), fileloop::once(Top), file_name(f_name), file_mode(f_mode) { }

    void file_loop() {
        if (*file_type) msgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

struct FLOOPboth: public fileloop::both {
    const char *file_name;
    const char *file_type;
    const char *file_mode;

    FLOOPboth(linkage_ped_top *Top, const char *f_name, const char *f_mode) :
        person_locus_entry(Top), fileloop::both(Top), file_name(f_name), file_mode(f_mode) { }

    void file_loop() {
        if (*file_type) mssgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

struct FLOOPchr: public fileloop::chr {
    const char *file_name;
    const char *file_type;
    const char *file_mode;

    FLOOPchr(linkage_ped_top *Top, const char *f_name, const char *f_mode) :
        person_locus_entry(Top), fileloop::chr(Top), file_name(f_name), file_mode(f_mode) { }

    void file_loop() {
        if (*file_type) mssgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

struct FLOOPtrait: public fileloop::trait {
    const char *file_name;
    const char *file_type;
    const char *file_mode;

    FLOOPtrait(linkage_ped_top *Top, const char *f_name, const char *f_mode) :
        person_locus_entry(Top), fileloop::trait(Top), file_name(f_name), file_mode(f_mode) { }

    void file_loop() {
        if (*file_type) mssgvf("%s%s/%s\n", file_type, *_opath, file_name);
        _dataloop->data_loop(*_opath, file_name, file_mode);
    }
};

#endif
