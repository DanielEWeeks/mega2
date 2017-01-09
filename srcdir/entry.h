/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#ifndef PERSON_LOCUS_ENTRY_H
#define PERSON_LOCUS_ENTRY_H

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "typedefs.h"

#include "create_summary_ext.h"
#include "output_routines_ext.h"

class person_locus_entry : public file_ops {
public:
    int _ped;
    int _per;
    int _locus;
    int _allele1, _allele2;
    linkage_ped_top  *_Top;
    linkage_locus_top  *_LTop;
    ext_linkage_locus_top *_EXLTop;
    linkage_ped_rec  *_tpersonp;
    linkage_ped_tree *_tpedtreep;
    linkage_locus_rec *_tlocusp;
    linkage_locus_rec *_ttraitp;
    int  _fwid;
    int  _pwid;
    int  _mwid;
    char _fformat[16];
    char _pformat[16];
    char _mformat[16];

    int _trait;

    int _numchr;
    int _chrom_loop;
public:
    person_locus_entry() {
        // It's really nice to initialize everything to '0' so that when
        // the code blows up, you know that the variable was never written
        // too as opposed to someone wrote some funny number there...
        init();
    }
    void init() {
        _ped = _per = _locus = _allele1 = _allele2 = 0;
        _Top    = (linkage_ped_top  *)NULL;
        _LTop   = (linkage_locus_top  *)NULL;
        _EXLTop = (ext_linkage_locus_top *)NULL;
        _tpersonp = (linkage_ped_rec  *)NULL;
        _tpedtreep = (linkage_ped_tree *)NULL;
        _tlocusp = (linkage_locus_rec *)NULL;
        _ttraitp = (linkage_locus_rec *)NULL;
        _fformat[0] = _pformat[0] = _mformat[0] = 0;
        _fwid = _pwid = _mwid = _trait = _numchr = _chrom_loop = 0;
    }
    person_locus_entry(linkage_ped_top *Top) {
        init();
        if (Top != (linkage_ped_top *)NULL) {
            this->_Top = Top;
            _LTop   = Top->LocusTop;
            _EXLTop = Top->EXLTop;
        }
    }
    virtual ~person_locus_entry() {
        if (_fln) {
            free(_fln);
            _fln = NULL;
        }
    }
    void load_formats_no_space(const int mwid);
    void load_formats(const int fwid, const int pwid, const int mwid);

    void setfln(const char *prefix, const char *sfx) 
    {
        char fa[FILENAME_LENGTH];
        sprintf(fa, "%s%s", prefix, sfx);
        _fln = strdup(fa);
    }
    void setfln(const char *prefix, const char *type, const char *sfx) 
    {
        char fa[FILENAME_LENGTH];
        sprintf(fa, "%s%s%s", prefix, type, sfx);
        _fln = strdup(fa);
    }
    void mkfln(char *buf, const char *prefix, const char *sfx) {
        extern void change_output_chr(char *buf, int _numchr);
        sprintf(buf, "%s%s", prefix, sfx);
        if (strchr(sfx+1, '.'))
            change_output_chr(buf, _numchr);
    }
    void mkfln(char *buf, const char *dir, const char *prefix, const char *sfx) {
        extern void change_output_chr(char *buf, int _numchr);
        sprintf(buf, "%s%s%s", dir, prefix, sfx);
        if (strchr(sfx+1, '.'))
            change_output_chr(buf+strlen(dir), _numchr);
    }

    void pr_id();
    void pr_uid();
    void pr_fam();
    void pr_per() { pr_per(_tpersonp); }
    void pr_per(linkage_ped_rec  *tpe);
    void pr_father();
    void pr_mother();
    void pr_parent(); // prints father then mother

    void pr_sex() { pr_sex(_tpersonp); }
    void pr_sex(linkage_ped_rec  *tpe);
    //natural language print sex
    void pr_sex_l(){ pr_sex_l(_tpersonp,"U");}
    void pr_sex_l(const char *unknown){ pr_sex_l(_tpersonp,unknown);}
    void pr_sex_l(linkage_ped_rec *tpe, const char *unkown);

    // returns -1 missing phenotype; 0 Control (unaffected); 1 Case (affected)
    int is_affected_pheno() { return is_affected_pheno(_tpersonp); }
    int is_affected_pheno(linkage_ped_rec  *tpe);

    int has_pheno() { return has_pheno(_tpersonp); }
    int has_pheno(linkage_ped_rec  *tpe);

    void pr_pheno() { pr_pheno(_tpersonp, 0); }
    void pr_pheno(const int affection_as_string) { pr_pheno(_tpersonp, affection_as_string); }
    void pr_pheno(linkage_ped_rec  *tpe) { pr_pheno(tpe, 0); }
    void pr_pheno(linkage_ped_rec  *tpe, const int affection_as_string);

    void pr_aff();
    void pr_quant();
    void pr_marker_name();

    // From loop.cpp it appears that _locus and _tle are always defined together...
    void pr_marker() { pr_marker(_locus, _tlocusp, _tpersonp); }
    void pr_marker(const int locus, linkage_locus_rec *tlocusp, linkage_ped_rec *tpersonp);

    void pr_marker_alleles();

    double get_genetic_distance(int *warnp);
    void pr_genetic_distance_warning(int warnp);
    int pr_genetic_distance(const char *format_string, const int display_warn);

    void pr_physical_distance(const char *format_string);

    const char *recode_name(const int allele, const char *zero, const char *other);
    linkage_locus_rec *get_quant_by_name(const char *name);
};

#endif
