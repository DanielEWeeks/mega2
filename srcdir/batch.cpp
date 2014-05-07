/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2014 Robert Baron, Nandita Mukhopadhyay,
  Lee Almasy, Mark Schroeder, William P. Mulvihill,
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

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "common.h"
#include "typedefs.h"
#include "batch.h"
#include "batch_input_ext.h"
#include "error_messages_ext.h"
#include "genetic_utils_ext.h"
#include "reorder_loci_ext.h"
#include "user_input_ext.h"
#include "plink_ext.h"
#include "utils_ext.h"
// for access()
#include "cw_routines_ext.h"


//
// The static (class variable) that holds the BatchFile singleton...
BatchFile *BatchFile::bfInstance = NULL;

/**
 It's hard to imagine that the standard library does not contain a trim function.
 This does a simple white space trim on a string.
 */
static void trim(string& str)
{
    string::size_type pos1 = str.find_first_not_of(" \t");
    string::size_type pos2 = str.find_last_not_of(" \t");
    str = str.substr(pos1 == string::npos ? 0 : pos1,
                     pos2 == string::npos ? str.length() - 1 : pos2 - pos1 + 1);
}

//
// A simple predicate for sorting batch items on their keyword only...
static bool BatchItemPointerKeywordCompare(BatchItem  *a, BatchItem *b) {
    return a->getKeyword() < b->getKeyword();
}

/**
 This initialization routine builds a list of all possible batch items that can
 be encountered in a batch file. This is necessary even if the item has never been
 read, because we need a place to store informaiton about the item in any event.
 Only BatchItem derived classes can (and should) be instantiated using the
 keyword/section pair. A default section is sepcified using the empty string.
 Default values can also be specified.
 
 See section 20.5 of the manual for details on the batch file itmes.
 
 NOTE: that this section has various rules for when an item is required
 and when it is not. This functionality has not been implemented to date.
 */
void BatchFile::initBatchItems() {
    // Details on batch file items...

    // These are required items in all cases...
    addBatchItem(new BatchItemFileInput("Input_Pedigree_File", ""));
    addBatchItem(new BatchItemFileInput("Input_Map_File", ""));
    
    // Rrequired item for linkage and Mega2 input format, but it is not used for either PLINK formats.
    addBatchItem(new BatchItemFileInput("Input_Locus_File", ""));
    
    // Required item for PLINK binary formats (e.g., the .BED file)
    // when BatchItemBoolean("bfile", "PLINK") is specified.
    addBatchItem(new BatchItemFileInput("Input_Binary_File", ""));
    
    // opitonal...
    addBatchItem(new BatchItemFileInput("Input_Omit_File", ""));
    
    // Optionally used to specify allele frequencies for Mega2 and PLINK formats.
    addBatchItem(new BatchItemFileInput("Input_Frequency_File", ""));
    
    // Optionally used to specify allele penetrances for Mega2 and PLINK formats.
    addBatchItem(new BatchItemFileInput("Input_Penetrance_File", ""));
    
    // Optionally used to specify additional phenotypes in PLINK format.
    addBatchItem(new BatchItemFileInput("Input_Phenotype_File", ""));

    // Required value of the 'Untyped Pedigree Option' menu...
    addBatchItem(new BatchItemIntGe0("Input_Untyped_Ped_Option", "", -1));
    
    // Optionally specifies taking a random genotyping error simulation step.
    addBatchItem(new BatchItemBoolean("Input_Do_Error_Sim", "", false));
    
    // Optionally determines the directory for creating the output files.
    // The current directory is the default value
    addBatchItem(new BatchItemDirectory("Output_Path", "", "."));

    // Required specifies the analysis option...
    addBatchItem(new BatchItemAnalysis("Analysis_Option", ""));
    // The check for validity of this sub option is performed later in the
    // checks section by calling the CLASS_ANALYAIS::sub_prog_name_to_sub_option() method
    // which should generate an error if the subanalysis is not valid.
    addBatchItem(new BatchItemString("Analysis_Sub_Option", ""));

    // The chromosome number selected for analysis from the Locus Reording Menu...
    // THe validity of the chromosome number is decided after the locus and map files are read
    // If unused or invalid Loci_Selected must be used.
    addBatchItem(new BatchItemIntGe0("Chromosome_Single", "", 1));
    
    // Determines the number of chromosomes to be considered from the "Chromosomes_Multiple"
    // list. It can be less than or equal to it's length.
    addBatchItem(new BatchItemIntGt0("Chromosomes_Multiple_Num", ""));
    
    // Optionally the list of chromosomes to be selected. Length of list should be at
    // least as great as "Chromosomes_Multiple_Num"...
    addBatchItem(new BatchItemIntVectorChromosome("Chromosomes_Multiple", ""));

    // If true then each chromosome gets it's own file...
    addBatchItem(new BatchItemBoolean("Loop_Over_Chromosomes", "", false));

    // Optionally the marker numbers to be selected for analysis terminated with 'e'
    // Markers can only be validated after marker data is read...
    addBatchItem(new BatchItemIntVectorEterm("Loci_Selected", ""));

    // Required and specifies the genetic map to be used for analysis...
    addBatchItem(new BatchItemInt("Value_Genetic_Distance_Index", "", -1));
    const int vgdstm[] = {0, 1, 2};
    const vector<int> ValueGeneticDistanceSexTypeMapLegalValues(vgdstm, vgdstm + (sizeof(vgdstm)/sizeof(int)));
    addBatchItem(new BatchItemIntOf("Value_Genetic_Distance_SexTypeMap", "", ValueGeneticDistanceSexTypeMapLegalValues, -1));
    
    addBatchItem(new BatchItemInt("Value_Base_Pair_Position_Index", "", -2)); // None (should this be -1?)
    
    // Value of the menu option of the Trait Reordering menu (positive integer).
    // Trait loci are numbered 1 - N.
    addBatchItem(new BatchItemIntGt0("Trait_Single", ""));
    
    // Value of the menu item 2 of the Trait Reordering menu where traits are selected
    // to be analyzed one at a time. Trait numbers should correspond to the order in
    // which they are listed int he locus file, starting from 1
    addBatchItem(new BatchItemIntVectorEterm("Traits_Loop_Over", ""));
    
    // Value of the menu item 3 of the Trait Reordering menu where the trait loci
    // are selected to be combined in the same output.
    addBatchItem(new BatchItemIntVectorEterm("Traits_Combine", ""));
    
    // Trait specific directories in which output files are created if the
    // "Traits_Loop_Over" option is defined.
    addBatchItem(new BatchItemStringVector("Trait_Subdirs", ""));
    
    // The missing quantitative phenotype value.
    addBatchItem(new BatchItemMissingQuant("Value_Missing_Quant_On_Input", "", 0.0));
    //addBatchItem(new BatchItemDeprecated("Value_Missing_Quant", ""));
    addBatchItem(new BatchItemMissingQuant("Value_Missing_Quant", ""));

    addBatchItem(new BatchItemString("Value_Missing_Quant_On_Output", "", "0"));
    
    // This is a list of status-class pairs which should be considered affected
    // each value should contain one entry per mulriple-liability trait.
    // The trait name is followed by the affection labels as follows
    // ValueAffecteds=TRAIT:2-1,2-3 trt2:2-*
    addBatchItem(new BatchItemStringVector("Value_Affecteds", ""));
    
    // Option defines a set of locus numbers at which random genotyping errors
    // will be introduced. This should be a subset of the markers chosen in
    // the locus reordering step.
    addBatchItem(new BatchItemIntVector("Error_Loci", ""));
    
    addBatchItem(new BatchItemIntVector("Error_Except_Loci", ""));
    
    addBatchItem(new BatchItemIntGe0("Error_Loci_Num", ""));

    // The value here is a single letter: U[niform], M[arker-specific], S[imWalk2].
    // batch_input.cpp:batchfile_init_Mega2BatchItems() uses 'X' as a default value.
    const string emlv[] = {"U", "M", "S"};
    const vector<string> errorModelLegalValues(emlv, emlv+(sizeof(emlv)/sizeof(string)));
    addBatchItem(new BatchItemStringOf("Error_Model", "", errorModelLegalValues, "X"));
    
    // Optionally specifies the error probability for error simulation.
    // can either be a single or a list of real values.
    addBatchItem(new BatchItemDoubleVector("Error_Probabilities", ""));

    // if true the Output File Names menu will be skipped for all options...
    addBatchItem(new BatchItemBoolean("Default_Outfile_Names", "", false));
    
    // if true genotypes are reset to unknown
    addBatchItem(new BatchItemBoolean("Default_Reset_Invalid", "", true));
    
    // if true default value parameters will be used for varoius menu options.
    addBatchItem(new BatchItemBoolean("Default_Other_Values", "", false));
    
    // if true do not exit on non-critical inconsistencies in the data.
    addBatchItem(new BatchItemBoolean("Default_Ignore_Nonfatal", "", false));
    
    // if true skip the R-plot parameters menu.
    addBatchItem(new BatchItemBoolean("Default_Rplot_Options", "", false));
    
    const int xlamlv[] = {0, 1, 2};
    const vector<int> xlinkedAnalysisModeLegalValues(xlamlv, xlamlv+(sizeof(xlamlv)/sizeof(int)));
    addBatchItem(new BatchItemIntOf("Xlinked_Analysis_Mode", "", xlinkedAnalysisModeLegalValues, 2));
    
    // trait loci that are output as covariates for the GeneHunter-format option
    addBatchItem(new BatchItemIntVector("Covariates_Selected", ""));
    
    // Values correspond to the choices in the individuals selection menu.
    // 1) Genotyped founders only
    // 2) Genotyped founders + a randomly chosen genotyped person from pedigrees without genotyped founders.
    // 3) Genotyped founders + genotyped individuals with unique alleles
    // *4) All genotyped individuals
    const int cglv[] = {1, 2, 3, 4};
    const vector<int> countGenotypesLegalValues(cglv, cglv+(sizeof(cglv)/sizeof(int)));
    addBatchItem(new BatchItemIntOf("Count_Genotypes", "", countGenotypesLegalValues, 4));
    
    // A list of numbers that refer to the list of statistis that are to be plotted.
    // The numbers match with the statistics menu displayed in the interactive mode.
    // The number (elements in the vector) of permitted statistics differs based on the analysis type:
    // TO_MERLIN (2), (5), (12).
    // Currently this is checked in R_output.cpp, but should be moved here.
    addBatchItem(new BatchItemIntVectorEterm("Rplot_Statistics", ""));
    
    // If true intividuals who are typed only at one allele are considered to be genotyped.
    addBatchItem(new BatchItemBoolean("Count_Halftyped", "", false));
    
    // Optional real value used to decide whether input allele-frequencies are different
    // from that observed in the data for any marker locus.
    addBatchItem(new BatchItemDoubleGe0("AlleleFreq_SquaredDev", "", 0.0));
    
    // Similar to "Count_Genotypes" but only considered withing the HWE estimation option.
    addBatchItem(new BatchItemInt("Count_HWE_genotypes", "", 2));
    
    // Optionally specifies the indicator for unknown allele and affection status values.
    addBatchItem(new BatchItemString("Value_Missing_Allele_Aff", ""));


    // PLINK command line option(s).
    /// See PLINK section for an example of replacing this line with a subsection.
    addBatchItem(new BatchItemString("PLINK", ""));
    //
    // These keywords are only valid in the PLINK section...
    
    addBatchItem(new BatchItemBoolean("file", "PLINK"));
    addBatchItem(new BatchItemBoolean("bfile", "PLINK"));
    
    // specify the genetic distance units...
    addBatchItem(new BatchItemBoolean("kosambi", "PLINK", true));
    addBatchItem(new BatchItemBoolean("haldane", "PLINK", false));
    
    // specify a name for the pedigree file trait column.
    addBatchItem(new BatchItemString("trait", "PLINK"));
    
    // is the trait quantitative or binary...
    addBatchItem(new BatchItemBoolean("quantitative", "PLINK"));
    addBatchItem(new BatchItemBoolean("affectionstatus", "PLINK"));
    
    addBatchItem(new BatchItemBoolean("no-fid", "PLINK", false));
    addBatchItem(new BatchItemBoolean("no-parentes", "PLINK", false));
    addBatchItem(new BatchItemBoolean("no-sex", "PLINK", false));
    addBatchItem(new BatchItemBoolean("no-pheno", "PLINK", false));
    addBatchItem(new BatchItemBoolean("map3", "PLINK", false));
    addBatchItem(new BatchItemBoolean("cM", "PLINK", false));
    addBatchItem(new BatchItemInt("missing-phenotype", "PLINK", -9));
    addBatchItem(new BatchItemBoolean("1", "PLINK", false));
    
    //
    // These keywords are only valid in the STRUCTURE section...
    
    //addBatchItem(new BatchItemString("PopDataPheno", "Structure"));
    addBatchItem(new BatchItemString("Structure.PopDataPheno", ""));

    
    // In v3.0 "Markers_Selected" was renamed to "Loci_Selected".
    addBatchItem(new BatchItemDeprecated("Markers_Selected", ""));
    
    addBatchItem(new BatchItemDeprecated("REMOutput_Map_Num", ""));
    addBatchItem(new BatchItemDeprecated("Output_Map_Num", ""));

    addBatchItem(new BatchItemDeprecated("Default_Ignore_Xlinked", ""));
    
    addBatchItem(new BatchItemDeprecated("Markers_Selected_Num", ""));
    addBatchItem(new BatchItemDeprecated("Loci_Selected_Num", ""));
    addBatchItem(new BatchItemDeprecated("REMTraits_Num", ""));
    
    addBatchItem(new BatchItemDeprecated("REMCovariates_Num", ""));
    addBatchItem(new BatchItemDeprecated("Covariates_Num", ""));

    
    // Finally we sort the entries for a binary search...
    sort(batchItems.begin(), batchItems.end(), BatchItemPointerKeywordCompare);
    sort(sections.begin(), sections.end());
}


/**
 This method will parse [read and process] the batch file.
 
 This is the guts of the batch file processing. It will get a line from the file
 and try to convert it into a batch file item if it is possible for it to do so.
 
 Differences between this version of batch file processing and the old version
 includes:
 1) Multiple occurrences of the keyword are flagged as an error, where is in the
 past a warning was generated and only the first definition was used.
 2) Keyword matches are done in a case insensitive manner with the entire word
 being compared and not just a few latters here and there.
 3) Acceptable boolean values are case insensitive and from the following set:
 for affirmative {YES, Y, TRUE, T, 1}, for negative {NO, N, FALSE, F, 0}.
 */
void BatchFile::parse(const string fileName) {
    string section = "";
    ifstream file(fileName.c_str());
    
    // The batch file should only be processed once...
    if (processed) {
        throw BatchFileException("The batch file '" + batchFileName + "' has already been processed.",
                                 BATCHFILE_FILE_HAS_BEEN_PROCESSED);
    }
    
    // Try to open the batch file for processing....
    if (!file.is_open()) {
        throw BatchFileException("The batch file '" + fileName + "' was not found.",
                                 BATCHFILE_FILE_NOT_FOUND);
    }
    
    // Record the name though it may not be used...
    batchFileName = fileName;
    
    // loop over the lines in the batch file trying to process them...
    while (!file.eof()) {
        string keyword, value, line;
        
        // The next line in the file...
        getline(file,line);
        trim(line);
        
        // The first line of the file may contain a version number indicator of the form...
        // #Version4.4
        // However, currently nothing is done with the version number of the batch file,
        // and so we won't bother looking for it...
        
        // skip blank or comment lines...
        if (line.length() == 0 || line.at(0) == COMMENT_CHAR) continue;
        
        // We are intersted in two forms of lines: those that contain keyword/value assignments
        // and those that begin new sections. Anything else that we find at this point is an error.
        size_t found = line.find_first_of("=");
        if (found != string::npos) {
            // This line should be of the form...
            // keyword = value # optional comment
            keyword = line.substr(0,found);
            trim(keyword);
            
            // Need to search for both the section and the keyword to find a valid match...
            vector<string> arg;
            arg.push_back(section);
            arg.push_back(keyword);
            
            vector<BatchItem *>::const_iterator itter =
            find_if(batchItems.begin(), batchItems.end(), bind2nd(BatchItemSectionKeywordEqual(), arg));
            // We should have been careful to instantiate a BatchItem for every section/keyword that we are
            // aware of. So not finding it in our list is a problem....
            if (itter == batchItems.end()) {
                if (section == "")
                    throw BatchFileException("The batch file keyword '" + keyword + "' is not valid.",
                                             BATCHFILE_BATCHITEM_INVALID_KEYWORD);
                else
                    throw BatchFileException("The batch file keyword '" + keyword + "' is not valid in section '" + section + "'.",
                                             BATCHFILE_BATCHITEM_INVALID_KEYWORD);
            }
            
            // So we know about the batch item keyword and section...
            BatchItem *item = *itter;
            
#ifndef DISABLE_EXCEPTION_IF_KEYWORD_DEPRECATED
            // Make sure that it is not deprecated...
            if (item->getIsDeprecated()) {
                if (section == "")
                    throw BatchFileException("Keyword '" + keyword + "' is deprecated.",
                                             BATCHFILE_DEPRECATED_ENTRY);
                else
                    throw BatchFileException("Section '" + section + "' with keyword '" + keyword + "' is deprecated.",
                                             BATCHFILE_DEPRECATED_ENTRY);
            }
#endif /* DISABLE_EXCEPTION_IF_KEYWORD_DEPRECATED */
            
            // Make sure that we havn't seen this one before...
            if (item->itemRead()){
                if (section == "")
                    throw BatchFileException("Duplicate entry for keyword '" + keyword + "'.",
                                             BATCHFILE_DUPLICATE_ENTRY);
                else
                    throw BatchFileException("Duplicate entry for section '" + section + "' with keyword '" + keyword + "'.",
                                             BATCHFILE_DUPLICATE_ENTRY);
            }
            
            // look for a comment on the line after the value, and remove if found...
            // IMPORTANT NOTE: we have no way of including the comment character (e.g., quoting it)
            // in a value. But this has not been a problem to this point.
            value = line.substr(found+1);
            found = value.find_first_of(COMMENT_CHAR);
            if (found != string::npos) value = value.substr(0,found);
            trim(value);
            
            // was there anything after the equal sign?
            if (value.length() == 0) {
                if (section == "")
                    throw BatchFileException("No value found in batch file for keyword '" + keyword + "'.",
                                             BATCHFILE_NO_VALUE_FOR_ENTRY);
                else
                    throw BatchFileException("No value found in batch file for section '" + section + "' and keyword '" + keyword + "'.",
                                             BATCHFILE_NO_VALUE_FOR_ENTRY);
            }
            
            // Finally we can set the value of the item...
            item->setValue(value);
            item->setWasReadFromBatchFile();
        } else if (line.at(0) == '[' && line.at(line.length()-1) == ']') {
            // This line should be of the form...
            // [ sectionName ]
            section = line.substr(1, line.length()-2);
            trim(section);
            // Make a record of having encountered this section
            if (find(sections.begin(), sections.end(), section) == sections.end()) sections.push_back(section);
        } else {
            throw BatchFileException("The batch file line '" + line + "' is not valid.",
                                     BATCHFILE_BATCHITEM_LINE_INVALID);
        }
    }
    
    file.close();
    
    processed = true;
}

/**
 Perform consistency checks over information that we have read from the batch file.
 */
void BatchFile::consistencyChecks() {

    //
    // Some batch items are required unconditionally...
    dynamic_cast<BatchItemAnalysis&>(BF()->getBatchItem("Analysis_Option")).setIsRequired();
    dynamic_cast<BatchItemFile&>(BF()->getBatchItem("Input_Pedigree_File")).setIsRequired();
    dynamic_cast<BatchItemFile&>(BF()->getBatchItem("Input_Map_File")).setIsRequired();
    // From user_input.cpp...
    dynamic_cast<BatchItemIntGe0&>(BF()->getBatchItem("Input_Untyped_Ped_Option")).setIsRequired();

    //
    // Warnings based on missing items...

#ifndef HIDESTATUS
    BatchItemFile f = dynamic_cast<BatchItemFile&>(BF()->getBatchItem("Input_Omit_File"));
    if (!f.itemRead())
      mssgvf("BatchItem '%s' not in batch file, assumed to be unspecified.\n",  f.getKeywordCstr());

    f = dynamic_cast<BatchItemFile&>(BF()->getBatchItem("Input_Frequency_File"));
    if (!f.itemRead())
      mssgvf("BatchItem '%s' not in batch file, assumed to be unspecified.\n",  f.getKeywordCstr());

    f = dynamic_cast<BatchItemFile&>(BF()->getBatchItem("Input_Penetrance_File"));
    if (!f.itemRead())
      mssgvf("BatchItem '%s' not in batch file, assumed to be unspecified.\n",  f.getKeywordCstr());

    f = dynamic_cast<BatchItemFile&>(BF()->getBatchItem("Input_Phenotype_File"));
    if (!f.itemRead())
      mssgvf("BatchItem '%s' not in batch file, assumed to be unspecified.\n",  f.getKeywordCstr());

    // Input_Binary_File is handled when the Analysis_Option PLINK is checked for...
#endif /* HIDESTATUS */

    //
    // Some batch items are requied conditionally...

    // From user_input.cpp...
#ifndef HIDESTATUS
    BatchItemString outputPath = dynamic_cast<BatchItemString&>(BF()->getBatchItem("Output_Path"));
    if (!outputPath.itemRead()) {
      mssgvf("Keyword %s not in batch file, using default '%s' (current directory).\n",
	     outputPath.getKeywordCstr(),
	     outputPath.getDefaultValueCstr()
	     );
    }
#endif /* HIDESTATUS */

    // From user_input.cpp...
    BatchItemDoubleGe0 squaredDev = dynamic_cast<BatchItemDoubleGe0&>(BF()->getBatchItem("AlleleFreq_SquaredDev"));
    if (!squaredDev.itemRead()) {
      mssgvf("Keyword %s not in batch file, using default = 'no limit'.\n", squaredDev.getKeywordCstr());
    }
    
    // If the user defined the Chromosome_Multiple vector, they are required to specify Chromosomes_Multiple_Num, and
    // it must be less than or equal to the number of elements in Chromocome_Multiple...
    BatchItemIntVector chromosomesMultiple = dynamic_cast<BatchItemIntVector&>(BF()->getBatchItem("Chromosomes_Multiple"));
    if (chromosomesMultiple.itemRead()) {
        BatchItemIntGt0 chromosomesMultipleNum = dynamic_cast<BatchItemIntGt0&>(BF()->getBatchItem("Chromosomes_Multiple_Num"));
        chromosomesMultipleNum.setIsRequired();
        
        if (!chromosomesMultipleNum.itemRead()) {
            errorf("The batch file item Chromosomes_Multiple requires the specification of the");
            errorf("batch file item Chromosomes_Multiple_Num.");
            EXIT(BATCH_FILE_ITEM_ERROR);
        }
        if ((int)chromosomesMultiple.size() < chromosomesMultipleNum.getValue()) {
            errorf("The batch file item Chromosomes_Multiple_Num must be less than or equal to the");
            errorf("length of the number of items in Chromocomes_Multiple");
            EXIT(BATCH_FILE_ITEM_ERROR);
        }
    }
    
    if (BF()->itemRead("Chromosome_Single") ||
        (BF()->itemRead("Chromosomes_Multiple_Num") && BF()->itemRead("Chromosomes_Multiple"))) {
        
        setBatchReorder();
        mssgf("Markers, chromosome(s) read in from batch file.");
    } else {
        warnf("Locus selections not specified in batch file.");
    }
    
#if 0
    if (dynamic_cast<BatchItemStringVector&>(BF()->getBatchItem("Trait_Subdirs")).itemRead())
      dynamic_cast<BatchItemIntVectorEterm&>(BF()->getBatchItem("Traits_Loop_Over")).setIsRequired();

    if (dynamic_cast<BatchItemDoubleVector&>(BF()->getBatchItem("Error_Probabilities")).itemRead())
      dynamic_cast<BatchItemStringOf&>(BF()->getBatchItem("Error_Model")).setIsRequired();
#endif
    
    // Error_Loci or Error_Except_Loci, both need Error_Loci_Num to be specified...
    if (dynamic_cast<BatchItemIntVector&>(BF()->getBatchItem("Error_Loci")).itemRead())
      dynamic_cast<BatchItemIntGe0&>(BF()->getBatchItem("Error_Loci_Num")).setIsRequired();
    if (dynamic_cast<BatchItemIntVector&>(BF()->getBatchItem("Error_Except_Loci")).itemRead())
      dynamic_cast<BatchItemIntGe0&>(BF()->getBatchItem("Error_Loci_Num")).setIsRequired();
    
    // These are Required for Linkage and Mega2 input format...
    //f = (dynamic_cast<BatchItemFile&>(BF()->getBatchItem("Input_Locus_File")));

    // Throw an exception of a required item was not read from the batch file...
    for (vector<BatchItem *>::iterator it = batchItems.begin(); it != batchItems.end(); ++it) {
        if ((*it)->getIsRequired() && !(*it)->itemRead()) {
            if ((*it)->getSection() == "")
                throw BatchFileException("The batch file item '" + (*it)->getKeyword() + "' is a required item.",
                                         BATCHFILE_BATCHITEM_REQUIRED);
            else
                throw BatchFileException("The batch file item '" + (*it)->getKeyword() + "' is a required item in section '" + (*it)->getSection() + "'.",
                                         BATCHFILE_BATCHITEM_REQUIRED);
        }
    }
    
    // THERE ARE A NUMBER OF CHECKS BASED ON THE ANALYSIS TYPE IN read_files.cpp:ext_linkage_locus_top() that could be placed here.
    
    BatchItemAnalysis bi_analysis = dynamic_cast<BatchItemAnalysis&>(BF()->getBatchItem("Analysis_Option"));
    CLASS_ANALYSIS *analysis = bi_analysis.getValue();
    BatchItemString bi_analysisSubOption = dynamic_cast<BatchItemString&>(BF()->getBatchItem("Analysis_Sub_Option"));

    if (!analysis->is_enabled()) {
	  throw BatchFileException("The 'Analysis_Option' ('" + string(analysis->_name) + "' is disabled.",
				   BATCHFILE_ANALYSIS_OPTION_DISABLED);
    }
    
    if (bi_analysisSubOption.itemRead() && analysis->has_sub_options()) {
        string sub_analysis = bi_analysisSubOption.getValue();
        analysis->sub_prog_name_to_sub_option((char *)sub_analysis.c_str(), &analysis);
        // Because the above could return a diffenent analysis object than you called it with we need to reset the class!!!!
        bi_analysis.setValue(analysis);
        analysis->_subname = (char *)sub_analysis.c_str();
    }
    
    BatchItemString bi_mqoo = dynamic_cast<BatchItemString&>(BF()->getBatchItem("Value_Missing_Quant_On_Output"));
    if (bi_mqoo.itemRead()) {
        if (!analysis->output_quant_can_define_missing_value()) {
            warnvf("'Analysis_Option' = '%s' does not allow the definition of 'Value_Missing_Quant_On_Output'.\n",
                   analysis->_name);
        }

        if (analysis->output_quant_must_be_numeric() && ! bi_mqoo.is_numeric()) {
            errorvf("'Analysis_Option' = '%s' requires a numeric value for 'Value_Missing_Quant_On_Output'.\n",
                    analysis->_name);
            EXIT(BATCH_FILE_ITEM_ERROR);
        }
    }

    // This needs to be removed at some point and Value_Missing_Quant marked as deprecated in the definition section...
    BatchItemMissingQuant mq = dynamic_cast<BatchItemMissingQuant&>(BF()->getBatchItem("Value_Missing_Quant"));
    if (mq.itemRead()) {
       BatchItemMissingQuant mqoi = dynamic_cast<BatchItemMissingQuant&>(BF()->getBatchItem("Value_Missing_Quant_On_Input"));
       if (mqoi.itemRead()) {
            errorf("Value_Missing_Quant is deprecated and used with Value_Missing_Quant_On_Input.");
            EXIT(BATCH_FILE_ITEM_ERROR);
       }
       warnf("Value_Missing_Quant is deprecated, use Value_Missing_Quant_On_Input.");
       mq.resetWasReadFromBatchFile();
       mqoi.setValue(mq.getValue());
       mqoi.setWasReadFromBatchFile();
    }

    BatchItemString bi_PLINK = dynamic_cast<BatchItemString&>(BF()->getBatchItem("PLINK"));
    if (bi_PLINK.itemRead()) {
        // Determine if there are any bad arguments in the "PLINK" BatchItem...
        string commandLineArgs = bi_PLINK.getValue();
        if (!PLINK_args((char *)commandLineArgs.c_str(), 1)) {
            throw BatchFileException("A bad parameter was detected in the PLINK batch file keyword string.",
                                     BATCHFILE_BATCHITEM_PLINK_BAD_PARAMETER);
        }
        
        // Make sure that the PLINK batch parameter '--missing-phenotype n' matches with the BatchItem 'Value_Missing_Quant'.
        if (PLINK.missing_pheno == 1) {
	    //BatchItemMissingQuant mq = dynamic_cast<BatchItemMissingQuant&>(BF()->getBatchItem("Value_Missing_Quant_On_Input"));
            if (mq.itemRead()) {
                if (PLINK.pheno_value != mq.getValue()) {
                    errorvf("When both batch file parameters 'Value_Missing_Quant_On_Input' (%f) and\n",
                            mq.getValue());
                    errorvf("'PLINK --missing-phenotype' (%f) are specified, they must be equal. \n",
                            PLINK.pheno_value);
                    EXIT(BATCH_FILE_ITEM_ERROR);
                }
            } else {
                // for consistency set the Value_Missing_Quant...
                mq.setValue(PLINK.pheno_value);
            }
        }
        
        // If the PLINK binary input file format is specified on the 'PLINK' BatchItem, then
        // the user must also have specified the binary input file...
        if (PLINK.plink == binary_PED_format && !BF()->getBatchItem("Input_Binary_File").itemRead()) {
            throw BatchFileException("A 'Input_Binary_File' must be specified when the 'PLINK' BatchItem contains a '-bfile' flag.",
                                     BATCHFILE_BATCHITEM_PLINK_BAD_PARAMETER);
        }
    }

    // Analysis_Option has a dependency on Analysis_Sub_Option unless it's PLINK, and then we fudge it
    // for backward compatability...
    if (TO_PLINK->_name == analysis->_name) {
        // If the default analysis is not specified, then make it 'lgen'...
        if (!bi_analysisSubOption.itemRead() && analysis->_suboption == -1) analysis->_suboption = 1; // 'lgen
    }

#ifndef DISABLE_EXCEPTION_ANALYSIS_DEPENDS_ON_SUBANALYSIS
    bi_analysis.checkDependencyOn(bi_analysisSubOption);
#endif /* DISABLE_EXCEPTION_ANALYSIS_DEPENDS_ON_SUBANALYSIS */
    
    // If the user defines Trait sub directories, the number of traits must be set using the trait-selection item...
    BF()->getBatchItem("Trait_Subdirs").checkDependencyOn(BF()->getBatchItem("Traits_Loop_Over"));
    
    BF()->getBatchItem("Error_Loci").checkDependencyOn(BF()->getBatchItem("Error_Loci_Num"));
    BF()->getBatchItem("Error_Except_Loci").checkDependencyOn(BF()->getBatchItem("Error_Loci_Num"));
    
    if (BF()->itemRead("Value_Affecteds")) {
        // Depends on one of these keywords to be defined...
        if (!(BF()->itemRead("Trait_Single") || BF()->itemRead("Traits_Loop_Over") || BF()->itemRead("Traits_Combine"))) {
            throw BatchFileException("BatchItem 'Value_Affecteds' depends on one of the following being defined: Trait_Single, Traits_Loop_Over, or Traits_Combine",
                                     BATCHFILE_BATCHITEM_DEPENDENCY_NOT_MET);
        }
    }

#if 0    
    // These are related and checking should be done here...
    BatchItemStringOf bi_errorModel = dynamic_cast<BatchItemStringOf&>(BF()->getBatchItem("Error_Model"));
    BF()->getBatchItem("Error_Probabilities").checkDependencyOn(bi_errorModel);
    string eMv = bi_errorModel.getValue();
    if (eMv == "U" || eMv == "u") {
    } else  if (eMv == "S" || eMv == "s") {
    }
#endif /* 0 */
    
    // Now that we have taken care of all of the conditions for which we have no menus...
    // take care of those for which we do...
    
    // If anything is missing present the user with a menu...
    if (BF()->itemRead("Input_Locus_File") || BF()->itemRead("PLINK")) {
        setBatchInputFiles();
        mssgf("Input filenames and missing value indicator read in from batch file.");
    }
    
    if (BF()->itemRead("Trait_Single") || BF()->itemRead("Traits_Loop_Over") || BF()->itemRead("Traits_Combine")) {
        setBatchTraitSelection();
        mssgf("Trait selection(s) read in from batch file.");
    } else {
        warnf("Trait selections not specified in batch file.");
    }
    
    if ((BF()->itemRead("Error_Loci") || BF()->itemRead("Error_Except_Loci")) && BF()->itemRead("Error_Loci_Num")) {
        setBatchErrorParameters();
        mssgf("Error parameters read in from batch file.");
    }
}

/**
 Search the dictionary of BatchItems and return the one associated with the keyword.
 */
#if 0
BatchItem& BatchFile::getBatchItem(const string keyword) {
    size_t s = batchItems.size();
    for (size_t i = 0; i < s; i++) {
        BatchItem *bi = batchItems[i];
        if (bi->getKeyword() == keyword) return *bi;
    }
    throw BatchFileException("No batch item found for keyword '" + keyword + "'",
                             BATCHFILE_BATCHITEM_NOT_FOUND);
}
#endif /* 0 */
BatchItem& BatchFile::getBatchItem(const string keyword) {
    vector<BatchItem *>::iterator itter =
    find_if(batchItems.begin(), batchItems.end(), bind2nd(BatchItemKeywordEqual(), keyword));
    if (itter != batchItems.end()) {
        return **itter;
    } else {
        throw BatchFileException("No batch item found for keyword '" + keyword + "'",
                                 BATCHFILE_BATCHITEM_NOT_FOUND);
    }
}
BatchItem& BatchFile::getBatchItem(const string keyword, const string section) {
    vector<string> arg;
    arg.push_back(section);
    arg.push_back(keyword);
    
    vector<BatchItem *>::const_iterator itter =
    find_if(batchItems.begin(), batchItems.end(), bind2nd(BatchItemSectionKeywordEqual(), arg));
    if (itter != batchItems.end()) {
        return **itter;
    } else {
        throw BatchFileException("No batch item found for section '" + section + "' keyword '" + keyword + "'",
                                 BATCHFILE_BATCHITEM_NOT_FOUND);
    }
}

/**
 Search the dictionary of BatchItems to determine if there is one associated with the keyword.
 If one is found, ask it if it was read from the batch file.
 */
bool BatchFile::itemRead(const string keyword) {
    vector<BatchItem *>::const_iterator itter =
        find_if(batchItems.begin(), batchItems.end(), bind2nd(BatchItemKeywordEqual(), keyword));
    return itter != batchItems.end() && (*itter)->itemRead();
}

bool BatchFile::itemRead(const string keyword, const string section) {
    vector<string> arg;
    arg.push_back(section);
    arg.push_back(keyword);
    
    vector<BatchItem *>::const_iterator itter =
    find_if(batchItems.begin(), batchItems.end(), bind2nd(BatchItemSectionKeywordEqual(), arg));
    return itter != batchItems.end();
}


/**
 Write a batch file that contains the current values.
 If a file corresponding to 'fileName' exists it will be over written.
 */
void BatchFile::overWrite(const string fileName) {
    ofstream bf;
    bf.open(fileName.c_str());
    
    bf << COMMENT_CHAR <<  "Version4.5\n";
    
    time_t   time_val;
#ifdef HIDEDATE
    struct tm tm = { 0, 14, 2, 29, 8-1, 1997-1900, 0, 0, 1};
    time_val = mktime(&tm);
#else
    time(&time_val);
#endif
    struct tm *time_stt = localtime (&time_val);
    bf << COMMENT_CHAR << "          " << asctime(time_stt);
    
    bf << COMMENT_CHAR << " \n";
    bf << COMMENT_CHAR << " A batch file contains keyword value pairs one per line. A comment is permitted after\n";
    bf << COMMENT_CHAR << " the '" << COMMENT_CHAR << "' character following a keyword value pair, or on a line without a keyword value pair.\n";
    bf << COMMENT_CHAR << " The comment character does not need to be the first character on a line.\n";
    bf << COMMENT_CHAR << "     #  For example this is a comment line\n";
    bf << COMMENT_CHAR << " foo = 17 # Comment following a keyword/value pair saying that seventeen is the value of foo\n";
    bf << COMMENT_CHAR << " \n";
    bf << COMMENT_CHAR << " There is a notion of sections (named and generic) as there is in a Windows like\n";
    bf << COMMENT_CHAR << " .INI file. Section demarcations appear on a line by themselves and consist of an optional\n";
    bf << COMMENT_CHAR << " section name enclosed in square brackets. For a named section, the name is not optional and\n";
    bf << COMMENT_CHAR << " is the section name. A generic section consists of empty brackets without a name. For backward \n";
    bf << COMMENT_CHAR << " compatibility the batch file begins with an implied generic section\n";
    bf << COMMENT_CHAR << " All keyword value pairs define following a section marker belong to that section.\n";
    bf << COMMENT_CHAR << " For example....\n";
    bf << COMMENT_CHAR << " Input_Locus_File = mylocus.txt # defined in the generic section\n";
    bf << COMMENT_CHAR << " [ PLINK ]\n";
    bf << COMMENT_CHAR << " cm = yes # a PLINK section keyword\n";
    bf << COMMENT_CHAR << " \n";
    bf << COMMENT_CHAR << " Currently implemented keywords:\n";
    
    for (vector<string>::iterator section_it = sections.begin(); section_it != sections.end(); ++section_it) {
        if (*section_it != "") bf << COMMENT_CHAR << " [" << *section_it << "]\n";
        for (vector<BatchItem *>::iterator bi_it = batchItems.begin(); bi_it != batchItems.end(); ++bi_it) {
            if (!(*bi_it)->getIsDeprecated() && (*bi_it)->getSection() == *section_it) {
                bf << COMMENT_CHAR << "   " << (*bi_it)->getKeyword() << "\n";
                //if ((*bi_it)->getHasDefaultValue()) bf << "=" << (*bi_it)->getDefaultValue() << "\n";
            }
        }
    }
    
    bf << COMMENT_CHAR << " \n";
    bf << COMMENT_CHAR << " Restrictions on usage:\n";
    bf << COMMENT_CHAR << "   Use either Chromosome_Single or \n";
    bf << COMMENT_CHAR << "      Chromosome_Single and Loci_Selected or\n";
    bf << COMMENT_CHAR << "      Chromosomes_Multiple and Chromsomes_Multiple_Num\n";
    bf << COMMENT_CHAR << " \n";
    bf << COMMENT_CHAR << "   Use either Trait_Single or\n";
    bf << COMMENT_CHAR << "              Traits_Loop_Over or\n";
    bf << COMMENT_CHAR << "              Traits_Combine.\n";
    bf << COMMENT_CHAR << "   Keyword    Trait_Subdirs is allowed only if\n";
    bf << COMMENT_CHAR << "   keyword    Traits_Loop_Over is defined.\n";
    bf << COMMENT_CHAR << " \n";
    bf << COMMENT_CHAR << "      Error_Loci and Error_Loci_Num or\n";
    bf << COMMENT_CHAR << "      Error_Except_Loci and Error_Loci_Num.\n";
    bf << COMMENT_CHAR << " \n";
    bf << COMMENT_CHAR << " Default settings :\n";
    bf << COMMENT_CHAR << " Default_Reset_Invalid: \n";
    bf << COMMENT_CHAR << "   \"yes\"= set inconsistent genotypes to 0 and continue.\n";
    bf << COMMENT_CHAR << "   \"no\"= continue without setting inconsistent genotypes to 0.\n";
    bf << COMMENT_CHAR << " Default_Ignore_Nonfatal:\n";
    bf << COMMENT_CHAR << "   Don't pause for other non-fatal errors in input:\n";
    bf << COMMENT_CHAR << " Default_Other_Values:\n";
    bf << COMMENT_CHAR << "    Use Mega2's default values instead of asking user\n";
    bf << COMMENT_CHAR << "   inside analysis-option menus.\n";
    bf << COMMENT_CHAR << "   Xlinked_Analysis_Mode:\n";
    bf << COMMENT_CHAR << " Xlinked_Analysis_Mode:\n";
    bf << COMMENT_CHAR << "   Set x-linked based on chromosome number(human).\n";
    bf << COMMENT_CHAR << "   treat all markers as autosomal,\n";
    bf << COMMENT_CHAR << " Count_Genotypes and Count_Halftyped\n";
    bf << COMMENT_CHAR << "   These correspond to the Select Individuals menu\n";
    bf << COMMENT_CHAR << "   within the allele-recoding step.\n";
    bf << COMMENT_CHAR << " \n";
    
    // output by sections...
    for (vector<string>::iterator section_it = sections.begin(); section_it != sections.end(); ++section_it) {
        if (*section_it != "") bf << "\n[" << *section_it << "]\n";
        for (vector<BatchItem *>::iterator bi_it = batchItems.begin(); bi_it != batchItems.end(); ++bi_it) {
            if (!(*bi_it)->getIsDeprecated() && (*bi_it)->getSection() == *section_it) {
                bf << (*bi_it)->getKeyword() << "=";
                bf << **bi_it << "\n";
            }
        }
    }
    
    bf.close();
}

/**
 Append to an existing batch file any values that have been replaced.
 */
void BatchFile::appendReplaced(const string fileName) {
    ofstream bf(fileName.c_str(), fstream::out | fstream::app);
    bool firstTimeComment = true;
    
    // output by sections...
    for (vector<string>::iterator section_it = sections.begin(); section_it != sections.end(); ++section_it) {
        for (vector<BatchItem *>::iterator bi_it = batchItems.begin(); bi_it != batchItems.end(); ++bi_it) {
            if (!(*bi_it)->getIsDeprecated() &&
                (*bi_it)->getSection() == *section_it &&
                (*bi_it)->getValueWasReplaced()) {
                if (firstTimeComment == true) {
                    firstTimeComment = false;
                    bf << " \n" << COMMENT_CHAR << " \n";
                    bf << COMMENT_CHAR << " The following batch items have been appended by Mega2.\n" << " \n";
                }
                bf << "[" << *section_it << "]\n";
                bf << (*bi_it)->getKeyword() << "=";
                bf << **bi_it << "\n";
            }
        }
    }
    
    bf.close();
    
}

void BatchItemDirectory::valueCheck() {
#ifndef DISABLE_FILE_CHECKS
  if (!(access(getValueCstr(), W_OK) == 0 || is_dir((char *)getValueCstr()) == 0))
    throw BatchFileException("The batch file keyword '" + getKeyword() + "' does not specify a writable directory.",
			     BATCHFILE_BATCHITEM_INVALID_DIRECTORY);
#endif /* DISABLE_FILE_CHECKS */
}

void BatchFile::CopyToMega2BatchItems() {
    int i;
    for (i=0; i<NUM_KEYS; i++) {
        char *keyword = Mega2BatchItems[i].keyword;
        int value_type = Mega2BatchItems[i].value_type;
        BatchItem& bi = BF()->getBatchItem(keyword);
        if (!bi.itemRead()) continue;
        // Do not process deprecated items...
        if (dynamic_cast<BatchItemDeprecated *>(&BF()->getBatchItem(keyword)) != NULL) continue;
        Mega2BatchItems[i].item_read=1;
        switch (value_type) {
            case INT:
                Mega2BatchItems[i].value.option = dynamic_cast<BatchItemInt&>(bi).getValue();
                break;
            case INT_LIST:
                Mega2BatchItems[i].value.mult_opts = dynamic_cast<BatchItemIntVector&>(bi).getValueClist();
                break;
            case FLOAT:
                Mega2BatchItems[i].value.fvalue = dynamic_cast<BatchItemDouble&>(bi).getValue();
                break;
            case FLOAT_LIST:
                Mega2BatchItems[i].value.mult_fvalues = dynamic_cast<BatchItemDoubleVector&>(bi).getValueClist();
                break;
            case YORN:
                Mega2BatchItems[i].value.copt = dynamic_cast<BatchItemBoolean&>(bi).getCharValue();
                break;
            case STRING:
                if (i == Analysis_Option) {
                    CLASS_ANALYSIS *anal = dynamic_cast<BatchItemAnalysis&>(bi).getValue();
                    strcpy(Mega2BatchItems[i].value.name, anal->_name);
                } else if (i == Loci_Selected || i == Traits_Loop_Over || i == Traits_Combine || i == Rplot_Statistics) {
                    strcpy(Mega2BatchItems[i].value.name, dynamic_cast<BatchItemIntVectorEterm&>(bi).getValueCstr());
                } else if (i == Covariates_Selected) {
                    strcpy(Mega2BatchItems[i].value.name, dynamic_cast<BatchItemIntVector&>(bi).getValueCstr());
                } else {
                    strcpy(Mega2BatchItems[i].value.name, dynamic_cast<BatchItemString&>(bi).getValueCstr());
                }
                break;
            case NAME_LIST:
                Mega2BatchItems[i].value.mult_names = (char **)dynamic_cast<BatchItemStringVector&>(bi).getValueClist();
                break;
            case CHAR:
                const char *str = dynamic_cast<BatchItemString&>(bi).getValueCstr();
                Mega2BatchItems[i].value.copt = str[0];
                break;
        }
    }
}


#ifdef NEW_BATCH

void batchfile_process(char *batch_file_name, analysis_type *analysis)
{
   try {
       BatchFile::BF()->ProcessBatchFile(batch_file_name);
       BatchFile::BF()->CopyToMega2BatchItems();
       BatchItem& bi = BatchFile::BF()->getBatchItem("Analysis_Option");
       *analysis = dynamic_cast<BatchItemAnalysis&>(bi).getValue();
       check_batch_items(); // This sets all of the  'batch??????' globals...
   } catch (BatchFileException &e) {
       // Here we give the user a textual description of the exception and exit with the appropriate   value...
     errorf("In processing Batch File:");
     errorf(e.what());
     EXIT(BATCH_FILE_ITEM_ERROR);
     //exit(e.getType());
   }
}

//
// Write the item associated with the Mega2BatchItem to the batch file.
void batchf(int item) {
    FILE *batchfp;
    int first_time = ((access(Mega2Batch, F_OK) != 0)? 1: 0);

    const char *keyword = Mega2BatchItems[item].keyword;
    BatchItem& bi = BatchFile::BF()->getBatchItem(keyword);

    if (first_time) create_batchfile();

    if ((batchfp = fopen(Mega2Batch, "a")) == NULL) {
        errorvf("Could not open Mega2 batch file (no write permission?).\n");
        EXIT(FILE_WRITE_ERROR);
    }

    if (first_time) {
        fprintf(batchfp, "%cVersion4.4\n", COMMENT_CHAR);
        fprintf(batchfp, "%c          ", COMMENT_CHAR);
        write_time(batchfp);
        batch_file_doc(batchfp);
    }

    // Everything but BatchItemDeprecated has a .getValueCstr() method
    // that should do the right thing, and if it's deprecated we shouldn't
    // be writing it's value to the batch file!
    fprintf(batchfp, "%s=%s\n", bi.getKeywordCstr(), bi.getValueCstr());

    fclose(batchfp);
}

#endif /* NEW_BATCH */
