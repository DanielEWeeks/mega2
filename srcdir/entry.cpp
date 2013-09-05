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

#include "file_ops.h"
#include "user_input_ext.h"
#include "entry.h"

void entry::load_formats(const int fwid, const int pwid, const int mwid)
{
    create_formats(fwid, pwid, _fformat, _pformat);
    _fwid = fwid;
    _pwid = pwid;

    if (mwid != -1) {
        sprintf(_mformat, "%%%ds", mwid);
        _mwid = mwid;
    }

}

void entry::pr_id()
{
    pr_fam();
    pr_per();
}

void entry::pr_fam()
{
    // Family or pedigree (see comments for OrigIds in common.h)...
    if (_fwid == 0) {
        errorvf("Internal Error: load_formats() method not called.\n");
        EXIT(SYSTEM_ERROR);
    } else if (OrigIds[1] == 2 || OrigIds[1] == 4) {
        pr_printf(_fformat, _tp->Name);
    } else if (OrigIds[1] == 3) {
        pr_printf(_fformat, _ped+1);
    } else {
        pr_printf(_fformat, _tp->Num);
    }
}

void entry::pr_per(linkage_ped_rec  *tpe)
{
    // Always print the individual info. If there is a father, and a mother
    // then print that information in the same manner as the individual
    // NOTE: PLUNK seems to like to print the individual info as <pedigree>_<individual>
    // as opposed to <individual>.
    if (_fwid == 0) {
        errorvf("Internal Error: load_formats() method not called.\n");
        EXIT(SYSTEM_ERROR);
    } else if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
        pr_printf(_pformat, tpe->OrigID);
    } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
        pr_printf(_pformat, tpe->UniqueID);
    } else {
        pr_printf(_pformat, tpe->ID);
    }
}

void entry::pr_father()
{
    if (_fwid == 0) {
        errorvf("Internal Error: load_formats() method not called.\n");
        EXIT(SYSTEM_ERROR);
    } else if (_tpe->Father != 0) {
        if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Father-1].OrigID);
        } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Father-1].UniqueID);
        } else {
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Father-1].ID);
        }
    } else {
    // If there was no father then print 0 (e.g., none)...
    /* create the 0 strings */
        int p;
        for(p = 0; p < (_pwid-1); p++) pr_printf(" ");
        pr_printf("0 ");
    }
}

void entry::pr_mother()
{
    if (_fwid == 0) {
        errorvf("Internal Error: load_formats() method not called.\n");
        EXIT(SYSTEM_ERROR);
    } else if (_tpe->Father != 0) {
        if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Mother-1].OrigID);
        } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Mother-1].UniqueID);
        } else {
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Mother-1].ID);
        }
    } else {
    // If there was no father (or mother) then print 0 (e.g., none)...
    /* create the 0 strings */
        int p;
        for(p = 0; p < (_pwid-1); p++) pr_printf(" ");
        pr_printf("0 ");
    }
}

void entry::pr_parent()
{
    if (_fwid == 0) {
        errorvf("Internal Error: load_formats() method not called.\n");
        EXIT(SYSTEM_ERROR);
    } else if (_tpe->Father != 0) {
        if ((OrigIds[0] == 1) || (OrigIds[0] == 2))  {
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Father-1].OrigID);
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Mother-1].OrigID);
        } else if ((OrigIds[0] == 3) || (OrigIds[0] == 4)) {
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Father-1].UniqueID);
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Mother-1].UniqueID);
        } else {
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Father-1].ID);
            pr_printf(_pformat,
                    _tp->Entry[_tpe->Mother-1].ID);
        }
    } else {
    // If there was no father (or mother) then print 0 (e.g., none)...
    /* create the 0 strings */
        int p;
        for(p = 0; p < (_pwid-1); p++) pr_printf(" ");
        pr_printf("0 ");
        for(p = 0; p < (_pwid-1); p++) pr_printf(" ");
        pr_printf("0 ");
    }
}

void entry::pr_sex(linkage_ped_rec  *tpe)
{
    // The sex of the individual...
    pr_printf("%1d ", tpe->Sex);
}

/**
   @Brief Determine if the individual has phenotype information
   
   @return 1 if the user does; 0 if they do not
*/
int entry::has_pheno(linkage_ped_rec  *tpe)
{
    if (_tte == (linkage_locus_rec *)NULL) return 0;

    switch(_tte->Type) {
        case AFFECTION:
            int ase;
            if (_tte->Pheno->Props.Affection.ClassCnt == 1)
                ase = tpe->Pheno[_trait].Affection.Status;
            else
                ase = aff_status_entry(tpe->Pheno[_trait].Affection.Status,
                                       tpe->Pheno[_trait].Affection.Class,
                                       _tte); //&(_Top->LocusTop->Locus[_trait]))
            
            if (ase == 0) return 0; // missing phenotype
            break;
            
        case QUANT:
            if (fabs(tpe->Pheno[_trait].Quant - MissingQuant) <= EPSILON) return 0; // missing phenotype
            break;

        default:
            break;
    }

    return 1; // there is a phenotype
}

/**
   @returns -1 missing phenotype; 0 Control (unaffected); 1 Case (affected)
 */
int entry::is_affected_pheno(linkage_ped_rec  *tpe)
{
    if (_tte == (linkage_locus_rec *)NULL) return -1;
    
    if (_tte->Type == AFFECTION) {
            int ase;
            // linkage.h:linkage_pedrec_data is a union (Affection(2xint),
            // Quant(2xint), Alleles (2xint), RAlleles (2xchar*)
            if (_tte->Pheno->Props.Affection.ClassCnt == 1)
                ase = tpe->Pheno[_trait].Affection.Status;
            else
                ase = aff_status_entry(tpe->Pheno[_trait].Affection.Status,
                                       tpe->Pheno[_trait].Affection.Class,
                                       _tte); //&(_Top->LocusTop->Locus[_trait]))
            
	    return ase-1;
    }
    
    return -1;
}

void entry::pr_pheno(linkage_ped_rec  *tpe, const int affection_as_string)
{
    int ase;

    switch(_tte->Type) {
        case AFFECTION:
            // linkage.h:linkage_pedrec_data is a union (Affection(2xint),
            // Quant(2xint), Alleles (2xint), RAlleles (2xchar*)
            if (_tte->Pheno->Props.Affection.ClassCnt == 1)
                ase = tpe->Pheno[_trait].Affection.Status;
            else
                ase = aff_status_entry(tpe->Pheno[_trait].Affection.Status,
                                       tpe->Pheno[_trait].Affection.Class,
                                       _tte); //&(_Top->LocusTop->Locus[_trait]))
            
            if (affection_as_string == 1) {
                // Currently only used for Eigenstrat...
                // See EIG4.2/CONVERTF/README for Case/Control definitions...
                // NOTE: Eigenstrat does not allow missing phenotypes in the .ped file
                // NOTE: For Eigenstrat non-quantitative (affection status) phenotype,
                // values must be either 'Case' or 'Control'. There seems to be no concept of 'missing'.
                if (ase == 0) pr_printf("NA "); // missing phenotype
                else if (ase == 1) pr_printf("Control "); // unaffected
                else if (ase == 2) pr_printf("Case "); // affected
            } else {
                pr_printf("%1d ", ase);
            }
            break;
            
        case QUANT:
            if (fabs(tpe->Pheno[_trait].Quant - MissingQuant) <= EPSILON) {
                pr_printf(" %s ", Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
            } else {
                pr_printf("%10.5f ", tpe->Pheno[_trait].Quant);
            }
            break;
        default:
            break;
    }
}

void entry::pr_aff()
{
    int ase;

    switch(_tte->Type) {
    case AFFECTION:
        // linkage.h:linkage_pedrec_data is a union (Affection(2xint),
        // Quant(2xint), Alleles (2xint), RAlleles (2xchar*)
        if (_tte->Pheno->Props.Affection.ClassCnt == 1)
            ase = _tpe->Pheno[_trait].Affection.Status;
        else
            ase = aff_status_entry(_tpe->Pheno[_trait].Affection.Status,
                                   _tpe->Pheno[_trait].Affection.Class,
                                   _tte); //&(_Top->LocusTop->Locus[_trait]))
        pr_printf("%1d ", ase);
        break;
    default:
        break;
    }
}

void entry::pr_quant()
{
    // This if statement protects the switch statement from invalid loop indexes (which should
    // in itself be sufficient) and also limits the locus to traits.
    switch(_tte->Type) {
    case QUANT:
        if (fabs(_tpe->Pheno[_trait].Quant - MissingQuant) <= EPSILON) {
            pr_printf(" %s ", Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name);
        } else {
            pr_printf("%10.5f ", _tpe->Pheno[_trait].Quant);
        }
        break;
    default:
        break;
    }
}


void entry::pr_marker_name()
{
    // The Loci name (in this case a SNP ID)...
    if (_fwid == 0) {
        errorvf("Internal Error: load_formats() method not called.\n");
        EXIT(SYSTEM_ERROR);
    } else
      //pr_printf(_mformat, _tte->Name);
      pr_printf(_mformat, _LTop->Locus[_locus].Name);
}

void entry::pr_marker(linkage_ped_rec  *tpe, const int locus)
{
    int a1, a2;
    get_2alleles(tpe->Marker, locus, &a1, &a2);
    pr_printf("%d %d ", a1, a2);
}

/**
 @brief Print all of the alleles associated with the marker locus
 */
void entry::pr_marker_alleles()
{
    int i;
    for (i=0; i <_tle->AlleleCnt; i++) {
        const char *name = _tle->Allele[i].name;
        if (name == (const char *)NULL) {
            // This problem seems to appear when _LTop->PedRecDataType == Premakeped
            // This base problem needs to be fixed properly.

            pr_printf("%d ", i + 1); //input was a linkage locus file with freq specified.

        } else
            pr_printf("%s ", name);
    }
}

/**
   @Brief get the genetic distance

   @param[out] int *warnp; if !NULL values are...
      -1, no warning;
       0, no map specified;
       1 average used on sex;
       2 autosome but sex map;
       3 X only female map.
   @return the genetic distnace
*/
double entry::get_genetic_distance(int *warnp)
{
    double genetic_distance = 0.0;
    
    if (warnp != (int *)NULL) *warnp = -1; // No warning
    
    if (genetic_distance_index < 0) {
        if (warnp != (int *)NULL) *warnp = 0; // No genetic map specified
    } else {
        // The routines that build the EXLTop structure should be making sure that
        // sex maps can only be of the form <a>, <m,f>, <f>.
        if (genetic_distance_sex_type_map == SEX_AVERAGED_GDMT) {
            // we were told to only use the average map, or only an average map was specified...
            genetic_distance = _EXLTop->EXLocus[_locus].positions[genetic_distance_index];
            // But tell the user about it if used for a sex linked chromosome...
//          if (_tle->Type == XLINKED || _tle->Type == YLINKED) 
            if (_tle->Marker->chromosome == SEX_CHROMOSOME || _tle->Marker->chromosome == MALE_CHROMOSOME) {
                if (warnp != (int *)NULL) *warnp = 1; // Using average position for a sex map
            }
        } else if (genetic_distance_sex_type_map == SEX_SPECIFIC_GDMT) {
            // Male, Female are available, so pick the right one given the current locus type...
//          if (_tle->Type == XLINKED) // X or SEX_CHROMOSOME
            if (_tle->Marker->chromosome == SEX_CHROMOSOME) // X or SEX_CHROMOSOME
                genetic_distance = _EXLTop->EXLocus[_locus].pos_female[genetic_distance_index];
//          else if (_tle->Type == YLINKED) // Y or MALE_CHROMOSOME
            else if (_tle->Marker->chromosome == MALE_CHROMOSOME)  // Y or MALE_CHROMOSOME
                genetic_distance = _EXLTop->EXLocus[_locus].pos_male[genetic_distance_index];
            else {
                if (warnp != (int *)NULL) *warnp = 2; // Autosome when a sex map specified
            }
        } else if (genetic_distance_sex_type_map == FEMALE_GDMT) {
//          if (_tle->Type == XLINKED) // X or SEX_CHROMOSOME
            if (_tle->Marker->chromosome == SEX_CHROMOSOME) // X or SEX_CHROMOSOME
                genetic_distance = _EXLTop->EXLocus[_locus].pos_female[genetic_distance_index];
            else {
                if (warnp != (int *)NULL) *warnp = 3; // Female map so only X chromosome analysis is possible?
            }
        }
    }    

    return genetic_distance;
}

//
// print something usefull from the warning returned by 'get_genetic_distance(int *warnp)'
void entry::pr_genetic_distance_warning(int warnp)
{
    switch (warnp) {
    case -1:
        // No error
        break;
    case 0:
        warnf("No genetic map specified.");
        break;
    case 1:
        warnf("Since only a sex-averaged genetic map was used, these positions have been used");
        warnvf("for markers (%d, %s) on the X and Y chromosomes.\n", _tle->Marker->chromosome, _tle->Name);
        break;
    case 2:
        // autozomes will be set to missing is a result of the decision made to allow the user
        // to select a SEX_SPECIFIC_GDMT
        warnf("Since a sex specific genetic map was used, only analysis on the");
        warnf("X & Y chromosome is possible.");
        break;
    case 3:
        warnf("Since only a female genetic map was used, only analysis on the");
        warnf("X chromosome is possible.");
        break;
	}
}

/**
   @brief Print the genetic distance associated with the marker locus

   The genetic distance will be printed according to the genetic_distance_index and
   genetic_distance_sex_type_map specified by the user. Checks will be made to see if the map
   type matches the locus type. If it does not match, and 'display_warn' is non-zero, a
   warning will be displayed. The routine will return a value based on warning given.

   @param[in] display_warn  If != 0 warnings will be displayed when anomolies are detected
   @param[in] format_string If not null, string used to format the genetic distance output
   @return 0 if no warning was displayed; >0 indicating that a certain warning was displayed
 */
int entry::pr_genetic_distance(const char *format_string, const int display_warn)
{
    int warnp;
    double genetic_distance = get_genetic_distance(&warnp);

    if (display_warn != 0) pr_genetic_distance_warning(warnp);

#if defined(_WIN) || defined(MINGW)
    pr_printf((format_string != NULL) ? format_string : "%.6f ", genetic_distance);
#else
    pr_printf((format_string != NULL) ? format_string : "%.6lf ", genetic_distance);
#endif
    return warnp;
}

/**
   @brief Print the physical distance associated with the marker locus

   The physical distance will be printed according to the base_pair_position (map)
   specified by the user.

   @param[in] format_string If not null, string used to format the genetic distance output
   @return void
*/
void entry::pr_physical_distance(const char *format_string)
{
    double base_pair_position = 0.0;

    if (base_pair_position_index >= 0) {
      base_pair_position = _EXLTop->EXLocus[_locus].positions[base_pair_position_index];
    }
#if defined(_WIN) || defined(MINGW)
    pr_printf((format_string != NULL) ? format_string : "%.0f ", base_pair_position);
#else
    pr_printf((format_string != NULL) ? format_string : "%.0lf ", base_pair_position);
#endif
}

const char *entry::recode_name(const int allele, const char *zero, const char *other)
{
    linkage_allele_rec *lar = _LTop->Locus[_locus].Allele;
    if (allele == 0) 
        return zero;
    else if (allele == 1)
        return lar[0].name;
    else if (allele == 2)
        return lar[1].name;
    else 
        return other;
}
