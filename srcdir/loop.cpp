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

#include "loop.h"

void loop::once::iterate()
{
    int tr;
    linkage_locus_rec *tte;

    _tte = (linkage_locus_rec *)NULL; // None was selected
    for (tr=0; tr < num_traits; tr++) {

        _trait = *_trp++;
        if (_trait == -1) continue;
        if (_trait == -99) break;

        tte = &(_LTop->Locus[_trait]);
        if (_trait_affect && tte->Type != AFFECTION)
            continue;
        if (_trait_quant  && tte->Type != QUANT)
            continue;
        _tte = tte;
        break;
    }

    make_file(); // which must call run_loop() to access an 'inner' loop
}

// Sets both the trait and the chromosome... outer == chr_trait
void loop::outer::iterate()
{
    int i;
    int tr;
    int   *sv_trp;
    char **sv_opath;
    linkage_locus_rec *tte;
    
    for (i = 0; i < main_chromocnt; i++) {

        if (_chrom_loop) {
            _numchr=global_chromo_entries[i];
            _Top->analysis->replace_chr_number(file_names, _numchr);
        } else {
            if (main_chromocnt > 1) {
                _numchr=-1;
                _Top->analysis->replace_chr_number(file_names, 0); // all
            } else {
                _numchr=global_chromo_entries[i];
                _Top->analysis->replace_chr_number(file_names, _numchr);
            }
        }

        sv_trp   = _trp;
        sv_opath = _opath;

        chr_start();

        _tte = (linkage_locus_rec *)NULL; // None was selected
        for (tr=0; tr < num_traits; tr++) {

            _trait = *_trp++;
            if (_trait == -1) continue;
            if (_trait == -99) break;

            tte = &(_LTop->Locus[_trait]);
            if (_trait_affect && tte->Type != AFFECTION)
                continue;
            if (_trait_quant  && tte->Type != QUANT)
                continue;

            _tte = tte;

            // A trait (_tte) must be selected to get here...
            trait_start();

            make_file();

            trait_end();

            if (!_trait_loop) break;

            _opath++;

        }

        if (_tte == 0) { // NO traits
            _opath = &output_paths[0];
            make_file();
        }

        chr_end();

        if (!_chrom_loop) break;

        _trp   = sv_trp;
        _opath = sv_opath;
    }
}

void loop::chr::iterate()
{
    int i, tr;
    linkage_locus_rec *tte;

    _tte = (linkage_locus_rec *)NULL; // None was selected
    for (tr=0; tr < num_traits; tr++) {

        _trait = *_trp++;
        if (_trait == -1) continue;
        if (_trait == -99) break;

        tte = &(_LTop->Locus[_trait]);
        if (_trait_affect && tte->Type != AFFECTION)
            continue;
        if (_trait_quant  && tte->Type != QUANT)
            continue;
        _tte = tte;
        break;
    }

    for (i = 0; i < main_chromocnt; i++) {

        if (_chrom_loop) {
            _numchr=global_chromo_entries[i];
            _Top->analysis->replace_chr_number(file_names, _numchr);
        } else {
            if (main_chromocnt > 1) {
                _numchr=-1;
                _Top->analysis->replace_chr_number(file_names, 0);
            } else {
                _numchr=global_chromo_entries[i];
                _Top->analysis->replace_chr_number(file_names, _numchr);
            }
        }

        chr_start();

        make_file(); // which must call run_loop() to access an 'inner' loop

        chr_end();

        if (!_chrom_loop) break;

    }
}

void loop::trait::iterate()
{
    int tr;
    linkage_locus_rec *tte;

    for (tr=0; tr < num_traits; tr++) {

        _trait = *_trp++;
        if (_trait == -1) continue;
        if (_trait == -99) break;

        tte = &(_LTop->Locus[_trait]);
        if (_trait_affect && tte->Type != AFFECTION)
            continue;
        if (_trait_quant  && tte->Type != QUANT)
            continue;

        _tte = tte;

        // A trait (_tte) must be selected to get here...
        trait_start();

        make_file();

        trait_end();

        if (!_trait_loop) break;

        _opath++;

    }

    if (_tte == 0) { // NO traits
        _opath = &output_paths[0];
        make_file();
    }

}

////////////////////////////////////////////////////////////////
//                       inner loops
////////////////////////////////////////////////////////////////

void loop::null::run_loop(const char *dir, const char *fl_name, const char *mode)
{
    filep_open(dir, fl_name, mode);

    inner();

    filep_close();
}

void loop::ped_per::run_loop(const char *dir, const char *fl_name, const char *mode)
{

    filep_open(dir, fl_name, mode);

    // for each pedigree (consulting the linkage_ped_top structure)...
    for (_ped=0; _ped < _Top->PedCnt; _ped++) {

        // It there is some indication as to why this pedigree should not be included, don't...
        if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;
        // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...

        _tp = &(_Top->Ped[_ped]);
        ped_start();

        for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
            // record for the individual (a linkage_ped_rec)...
            _tpe = &(_tp->Entry[_per]);
            //_tpe = &(_Top->Ped[_ped].Entry[_per]);
            inner();
        }
        ped_end();
    }
    filep_close();
}

void loop::ped_per_trait::run_loop(const char *dir, const char *fl_name, const char *mode)
{
    int tr;
    int *retrp;

    filep_open(dir, fl_name, mode);

    // for each pedigree (consulting the linkage_ped_top structure)...
    for (_ped=0; _ped < _Top->PedCnt; _ped++) {

        // It there is some indication as to why this pedigree should not be included, don't...
        if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;
        // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...

        _tp = &(_Top->Ped[_ped]);
        ped_start();

        for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
            // record for the individual (a linkage_ped_rec)...
            _tpe = &(_tp->Entry[_per]);
            //_tpe = &(_Top->Ped[_ped].Entry[_per]);
            per_start();

            retrp = trp;

            // output all of the traits on the same line...
            for (tr = 0; tr < num_traits; tr++) {
                _trait = *retrp++;
                // If this is a marker, skip it...
	        if (_trait == -1) continue;
                if (_trait == -99) break;

                _tte = &(_LTop->Locus[_trait]);
                if (trait_affect && _tte->Type != AFFECTION)
                    continue;
                if (trait_quant  && _tte->Type != QUANT)
                    continue;

                inner();
            }
            per_end();
        }
        ped_end();
    }
    filep_close();
}

void loop::trait_ped_per::run_loop(const char *dir, const char *fl_name, const char *mode)
{
    int tr;
    int *retrp;

    filep_open(dir, fl_name, mode);

    retrp = trp;
    for (tr = 0; tr < num_traits; tr++) {
      _trait = *retrp++;

      // If this is a marker, skip it...
      if (_trait == -1) continue;

      if (_trait == -99) break;

      _tte = &(_LTop->Locus[_trait]);
      if (trait_affect && _tte->Type != AFFECTION)
	continue;
      if (trait_quant  && _tte->Type != QUANT)
	continue;

      trait_start();

      // for each pedigree (consulting the linkage_ped_top structure)...
      for (_ped=0; _ped < _Top->PedCnt; _ped++) {

        // It there is some indication as to why this pedigree should not be included, don't...
        if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;

        // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...
        _tp = &(_Top->Ped[_ped]);
	ped_start();

        for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
            // record for the individual (a linkage_ped_rec)...
            _tpe = &(_tp->Entry[_per]);
            //_tpe = &(_Top->Ped[_ped].Entry[_per]);
	    
	    inner();
	}
	ped_end();
      }
      trait_end();
    }
    filep_close();
}

void loop::ped_per_loci::run_loop(const char *dir, const char *fl_name, const char *mode)
{
    int m;

    // Loop through the Loci...
    markers_on_chromosome(_numchr);

    filep_open(dir, fl_name, mode);

    // for each pedigree (consulting the linkage_ped_top structure)...
    for (_ped=0; _ped < _Top->PedCnt; _ped++) {

        // It there is some indication as to why this pedigree should not be included, don't...
        if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;
        // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...

        _tp = &(_Top->Ped[_ped]);
        ped_start();

        for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
            // record for the individual (a linkage_ped_rec)...

            _tpe = &(_Top->Ped[_ped].Entry[_per]);

            per_start();

            for (m=0; m < NumChrLoci; m++) {
                // The actual Loci number associated with a SELECTED Loci...
                _locus = ChrLoci[m];
                _tle = &(_LTop->Locus[_locus]);
                
                // Ignore anything but a marker.
                // This should get us things of Locus[?].Type == {NUMBERED, BINARY, XLINKED, YLINKED}.
                if (_tle->Class == MARKER) {
                    
                    if (_loci_allele_limit && _tle->AlleleCnt > _loci_allele_limit) continue;
                    
                    // The markers for this individual...
                    _allele1 = _tpe->Data[_locus].Alleles.Allele_1;
                    _allele2 = _tpe->Data[_locus].Alleles.Allele_2;

                    inner();
                }
            }
            per_end();
        }
        ped_end();
    }
    filep_close();
}

void loop::loci_ped_per::run_loop(const char *dir, const char *fl_name, const char *mode)
{
    int m;

    // Loop through the SELECTED Loci (by count)...
    markers_on_chromosome(_numchr);

    filep_open(dir, fl_name, mode);

    for (m=0; m < NumChrLoci; m++) {
        // The actual Loci number associated with a SELECTED Loci...
        _locus = ChrLoci[m];
        _tle = &(_LTop->Locus[_locus]);

        // Ignore anything but a marker.
        // This should get us things of Locus[?].Type == {NUMBERED, BINARY, XLINKED, YLINKED}.
        if (_tle->Class == MARKER) {

            if (_loci_allele_limit && _tle->AlleleCnt > _loci_allele_limit) continue;

            loci_start();

            // for each pedigree (consulting the linkage_ped_top structure)...
            for (_ped=0; _ped < _Top->PedCnt; _ped++) {

                // It there is some indication as to why this pedigree should not be included, don't...
                if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;
                // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...

                _tp = &(_Top->Ped[_ped]);
                ped_start();

                for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
                    // record for the individual (a linkage_ped_rec)...
                    _tpe = &(_tp->Entry[_per]);
                    //_tpe = &(_Top->Ped[_ped].Entry[_per]);

                    _allele1 = _tpe->Data[_locus].Alleles.Allele_1;
                    _allele2 = _tpe->Data[_locus].Alleles.Allele_2;

                    inner();
                }
                ped_end();
            }
            loci_end();
        }
    }
    filep_close();
}

void loop::loci::run_loop(const char *dir, const char *fl_name, const char *mode)
{
    int m;

    // Loop through the SELECTED Loci (by count)...
    markers_on_chromosome(_numchr);

    filep_open(dir, fl_name, mode);

    for (m=0; m < NumChrLoci; m++) {
        // The actual Loci number associated with a SELECTED Loci...
        _locus = ChrLoci[m];
        _tle = &(_LTop->Locus[_locus]);

        // Ignore anything but a marker.
        // This should get us things of Locus[?].Type == {NUMBERED, BINARY, XLINKED, YLINKED}.
        if (_tle->Class == MARKER) {

            if (_loci_allele_limit && _tle->AlleleCnt > _loci_allele_limit) continue;

            inner();
        }
    }
    filep_close();
}
