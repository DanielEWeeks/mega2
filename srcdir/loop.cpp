/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2016 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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
#include "tod.hh"
#include "analysis.h"

#include "output_file_names_ext.h"

void fileloop::once::iterate()
{
    int tr;
    linkage_locus_rec *tte;

    _ttraitp = (linkage_locus_rec *)NULL; // None was selected
    for (tr=0; tr < num_traits; tr++) {

        _trait = *_trp++;
        if (_trait == -1) continue;
        if (_trait == -99) break;

        tte = &(_LTop->Locus[_trait]);
        if (_trait_affect && tte->Type != AFFECTION)
            continue;
        if (_trait_quant  && tte->Type != QUANT)
            continue;
        _ttraitp = tte;
        break;
    }

    file_loop(); // which must call data_loop() to access an 'inner' loop
}

// Sets both the trait and the chromosome... outer == chr_trait
void fileloop::both::iterate()
{
    int i;
    int tr;
    int   *sv_trp;
    char **sv_opath;
    linkage_locus_rec *tte;
    
    for (i = 0; i < main_chromocnt; i++) {

        if (_chrom_loop) {
            _numchr=global_chromo_entries[i];
            _Top->analysis->replace_chr_number(Outfile_Names, _numchr);
        } else {
            if (main_chromocnt > 1) {
                _numchr=-1;
                _Top->analysis->replace_chr_number(Outfile_Names, 0); // all
            } else {
                _numchr=global_chromo_entries[i];
                _Top->analysis->replace_chr_number(Outfile_Names, _numchr);
            }
        }
        if (_fln) change_output_chr(_fln, _numchr);

        sv_trp   = _trp;
        sv_opath = _opath;

        chr_start();

        _ttraitp = (linkage_locus_rec *)NULL; // None was selected
        for (tr=0; tr < num_traits; tr++) {

            _trait = *_trp++;
            if (_trait == -1) continue;
            if (_trait == -99) break;

            tte = &(_LTop->Locus[_trait]);
            if (_trait_affect && tte->Type != AFFECTION)
                continue;
            if (_trait_quant  && tte->Type != QUANT)
                continue;

            _ttraitp = tte;

            // A trait (_ttraitp) must be selected to get here...
            trait_start();

            file_loop();

            trait_end();

            if (!_trait_loop) break;

            _opath++;

        }

        if (_ttraitp == 0) { // NO traits
            _opath = &output_paths[0];
            file_loop();
        }

        chr_end();

        if (!_chrom_loop) break;

        _trp   = sv_trp;
        _opath = sv_opath;
    }
}

void fileloop::chr::iterate()
{
    int i, tr;
    linkage_locus_rec *tte;

    _ttraitp = (linkage_locus_rec *)NULL; // None was selected
    for (tr=0; tr < num_traits; tr++) {

        _trait = *_trp++;
        if (_trait == -1) continue;
        if (_trait == -99) break;

        tte = &(_LTop->Locus[_trait]);
        if (_trait_affect && tte->Type != AFFECTION)
            continue;
        if (_trait_quant  && tte->Type != QUANT)
            continue;
        _ttraitp = tte;
        break;
    }

    for (i = 0; i < main_chromocnt; i++) {

        if (_chrom_loop) {
            _numchr=global_chromo_entries[i];
            _Top->analysis->replace_chr_number(Outfile_Names, _numchr);
        } else {
            if (main_chromocnt > 1) {
                _numchr=-1;
                _Top->analysis->replace_chr_number(Outfile_Names, 0);
            } else {
                _numchr=global_chromo_entries[i];
                _Top->analysis->replace_chr_number(Outfile_Names, _numchr);
            }
        }
        if (_fln) change_output_chr(_fln, _numchr);

        chr_start();

        file_loop(); // which must call data_loop() to access an 'inner' loop

        chr_end();

        if (!_chrom_loop) break;

    }
}

void fileloop::trait::iterate()
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

        _ttraitp = tte;

        // A trait (_ttraitp) must be selected to get here...
        trait_start();

        file_loop();

        trait_end();

        if (!_trait_loop) break;

        _opath++;

    }

    if (_ttraitp == 0) { // NO traits
        _opath = &output_paths[0];
        file_loop();
    }

}

////////////////////////////////////////////////////////////////
//                       inner loops
////////////////////////////////////////////////////////////////

void dataloop::null::data_loop(const char *dir, const char *fl_name, const char *mode)
{
    data_init();

    filep_open(dir, fl_name, mode);

    inner();

    filep_close();
}

void dataloop::ped_per::data_loop(const char *dir, const char *fl_name, const char *mode)
{
    data_init();

    filep_open(dir, fl_name, mode);

    // for each pedigree (consulting the linkage_ped_top structure)...
    for (_ped=0; _ped < _Top->PedCnt; _ped++) {

        // It there is some indication as to why this pedigree should not be included, don't...
        if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;
        // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...

        _tpedtreep = &(_Top->Ped[_ped]);
        ped_start();

        for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
            // record for the individual (a linkage_ped_rec)...
            _tpersonp = &(_tpedtreep->Entry[_per]);
            //_tpersonp = &(_Top->Ped[_ped].Entry[_per]);
            inner();
        }
        ped_end();
    }
    filep_close();
}

void dataloop::ped_per_trait::data_loop(const char *dir, const char *fl_name, const char *mode)
{
    int tr;
    int *retrp;

    data_init();

    filep_open(dir, fl_name, mode);

    // for each pedigree (consulting the linkage_ped_top structure)...
    for (_ped=0; _ped < _Top->PedCnt; _ped++) {

        // It there is some indication as to why this pedigree should not be included, don't...
        if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;
        // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...

        _tpedtreep = &(_Top->Ped[_ped]);
        ped_start();

        for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
            // record for the individual (a linkage_ped_rec)...
            _tpersonp = &(_tpedtreep->Entry[_per]);
            //_tpersonp = &(_Top->Ped[_ped].Entry[_per]);
            per_start();

            retrp = trp;

            // output all of the traits on the same line...
            for (tr = 0; tr < num_traits; tr++) {
                _trait = *retrp++;
                // If this is a marker, skip it...
	        if (_trait == -1) continue;
                if (_trait == -99) break;

                _ttraitp = &(_LTop->Locus[_trait]);
                if (trait_affect && _ttraitp->Type != AFFECTION)
                    continue;
                if (trait_quant  && _ttraitp->Type != QUANT)
                    continue;

                inner();
            }
            per_end();
        }
        ped_end();
    }
    filep_close();
}

void dataloop::trait_ped_per::data_loop(const char *dir, const char *fl_name, const char *mode)
{
    int tr;
    int *retrp;

    data_init();

    filep_open(dir, fl_name, mode);

    retrp = trp;
    for (tr = 0; tr < num_traits; tr++) {
      _trait = *retrp++;

      // If this is a marker, skip it...
      if (_trait == -1) continue;

      if (_trait == -99) break;

      _ttraitp = &(_LTop->Locus[_trait]);
      if (trait_affect && _ttraitp->Type != AFFECTION)
	continue;
      if (trait_quant  && _ttraitp->Type != QUANT)
	continue;

      trait_start();

      // for each pedigree (consulting the linkage_ped_top structure)...
      for (_ped=0; _ped < _Top->PedCnt; _ped++) {

        // It there is some indication as to why this pedigree should not be included, don't...
        if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;

        // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...
        _tpedtreep = &(_Top->Ped[_ped]);
	ped_start();

        for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
            // record for the individual (a linkage_ped_rec)...
            _tpersonp = &(_tpedtreep->Entry[_per]);
            //_tpersonp = &(_Top->Ped[_ped].Entry[_per]);
	    
	    inner();
	}
	ped_end();
      }
      trait_end();
    }
    filep_close();
}

void dataloop::ped_per_loci::data_loop(const char *dir, const char *fl_name, const char *mode)
{
    int m;

    data_init();

    // Dataloop through the Loci...
    markers_on_chromosome(_numchr);

    filep_open(dir, fl_name, mode);

    // for each pedigree (consulting the linkage_ped_top structure)...
    for (_ped=0; _ped < _Top->PedCnt; _ped++) {

        // It there is some indication as to why this pedigree should not be included, don't...
        if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;
        // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...

        _tpedtreep = &(_Top->Ped[_ped]);
        ped_start();

        for (_per = 0; _per < _Top->Ped[_ped].EntryCnt; _per++) {
            // record for the individual (a linkage_ped_rec)...

            _tpersonp = &(_Top->Ped[_ped].Entry[_per]);

            per_start();

            for (m=0; m < NumChrLoci; m++) {
                // The actual Loci number associated with a SELECTED Loci...
                _locus = ChrLoci[m];
                _tlocusp = &(_LTop->Locus[_locus]);
                
                // Ignore anything but a marker.
                // This should get us things of Locus[?].Type == {NUMBERED, BINARY, XLINKED, YLINKED}.
                if (_tlocusp->Class == MARKER) {
                    
                    if (_loci_allele_limit && _tlocusp->AlleleCnt > _loci_allele_limit) continue;
                    
                    // The markers for this individual...
                    get_2alleles(_tpersonp->Marker, _locus, &_allele1, &_allele2);

                    inner();
                }
            }
            per_end();
        }
        ped_end();
    }
    filep_close();
}

void dataloop::loci_ped_per::data_loop(const char *dir, const char *fl_name, const char *mode)
{
    int m;

    data_init();

    // Loop through the SELECTED Loci (by count)...
    Tod tod_lpp_markers("get markers_on_chromosomes");
    markers_on_chromosome(_numchr);
    tod_lpp_markers();

    filep_open(dir, fl_name, mode);
 
    Tod tod_lpp_start(10);
    for (m=0; m < NumChrLoci; m++) {
        // The actual Loci number associated with a SELECTED Loci...
        _locus = ChrLoci[m];
        _tlocusp = &(_LTop->Locus[_locus]);

        // Ignore anything but a marker.
        // This should get us things of Locus[?].Type == {NUMBERED, BINARY, XLINKED, YLINKED}.
        if (_tlocusp->Class == MARKER) {

            if (_loci_allele_limit && _tlocusp->AlleleCnt > _loci_allele_limit) continue;

            loci_start();

            tod_lpp_start.reset();
            // for each pedigree (consulting the linkage_ped_top structure)...
            for (_ped=0; _ped < _Top->PedCnt; _ped++) {

                // It there is some indication as to why this pedigree should not be included, don't...
                if (UntypedPeds != NULL && UntypedPeds[_ped]) continue;
                // for each person/inividual in that pedigree (cousulting the linkage_ped_tree)...

                _tpedtreep = &(_Top->Ped[_ped]);
                ped_start();

                _tpersonp = _tpedtreep->Entry;
                for (_per = 0; _per < _tpedtreep->EntryCnt; _per++, _tpersonp++) {
                    // record for the individual (a linkage_ped_rec)...
                    //_tpersonp = &(_Top->Ped[_ped].Entry[_per]);

                    get_2alleles(_tpersonp->Marker, _locus, &_allele1, &_allele2);
                    inner();
                }
                ped_end();
            }
            loci_end();
            tod_lpp_start("loop over ped person");

        }
    }
    filep_close();
}

void dataloop::loci::data_loop(const char *dir, const char *fl_name, const char *mode)
{
    int m;

    data_init();

    // Loop through the SELECTED Loci (by count)...
    markers_on_chromosome(_numchr);

    filep_open(dir, fl_name, mode);

    for (m=0; m < NumChrLoci; m++) {
        // The actual Loci number associated with a SELECTED Loci...
        _locus = ChrLoci[m];
        _tlocusp = &(_LTop->Locus[_locus]);

        // Ignore anything but a marker.
        // This should get us things of Locus[?].Type == {NUMBERED, BINARY, XLINKED, YLINKED}.
        if (_tlocusp->Class == MARKER) {

            if (_loci_allele_limit && _tlocusp->AlleleCnt > _loci_allele_limit) continue;

            inner();
        }
    }
    filep_close();
}
