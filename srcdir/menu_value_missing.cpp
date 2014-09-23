/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2014 Robert Baron, Charles P. Kollar,
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


#include <stdio.h>
#include <string.h>
#include <errno.h>

#include "common.h"
#include "typedefs.h"

#include "batch_input.h"

#include "batch_input_ext.h"
#include "error_messages_ext.h"
#include "fcmap_ext.h"
#include "plink_ext.h"
#include "user_input_ext.h"
#include "utils_ext.h"

/*
     error_messages_ext.h:  errorf mssgf my_calloc warnf
              utils_ext.h:  EXIT draw_line get_line
*/

char Str_Missing_Quant_On_Input[50]   = "";
char Str_Missing_Affect_on_Input[50]  = "";
char Str_Missing_Quant_On_Output[50]  = "";
char Str_Missing_Affect_On_Output[50] = "";
char Str_read[50] = "";

// If the Value_Missing_Quant_On_Output is either not specified or was not permitted
// to have been specified in the batch file (e.g., output_quant_can_define_missing_value()
// returns false), then check to see if there is a default value for that analysis mode (e.g.,
// output_quant_default_value() != NULL). If so, then make it look as if
// Value_Missing_Quant_On_Output was read as that default.
static void fix_Value_Missing_Quant_On_Input(analysis_type *analysis, struct itl *itp, const char *value);
static void fix_Value_Missing_Quant_On_Output(analysis_type *analysis, struct itl *itp, const char *value);
static void fix_Value_Missing_Affect_On_Input(analysis_type *analysis, struct itl *itp, const char *value);
static void fix_Value_Missing_Affect_On_Output(analysis_type *analysis, struct itl *itp, const char *value);
static int fix_Value_Missing_check_allow(analysis_type *analysis, int vmidx);
static int fix_Value_Missing_check_numeric(analysis_type *analysis, struct itl *itp, const char *value);
static int fix_Value_Missing(analysis_type *analysis, int vmidx, int init_str);
static void fix_Value_Missing_default(analysis_type *analysis, int vmidx);
static int fix_Value_Missing_all(analysis_type *analysis, int init_str);
static void Value_Missing_menu(analysis_type *analysis);

struct itl {
    int it;
    int inherit;
    char *str;
    const char *name;
    const char *put;
    void (*store)(analysis_type *analysis, struct itl *itp, const char *value);
    int source;
    int set;
};

static const char *source_name[] = { "menu input", "in batch", "plink --missing-pheno", "plink default", "analysis rule", "inherited", "cannot change", "output default", "linkage default", "mega2 default"};

struct itl Value_Missing[]  =  {
    {Value_Missing_Quant_On_Input,   0,                             Str_Missing_Quant_On_Input,
     "Quantitative", "Input",   fix_Value_Missing_Quant_On_Input, 0, 0},

    {Value_Missing_Quant_On_Output,  Value_Missing_Quant_On_Input,  Str_Missing_Quant_On_Output,
     "Quantitative", "Output",  fix_Value_Missing_Quant_On_Output, 0, 0},

    {Value_Missing_Affect_On_Input,  Value_Missing_Quant_On_Input,  Str_Missing_Affect_on_Input,
     "Affection", "Input",         fix_Value_Missing_Affect_On_Input, 0, 0},

    {Value_Missing_Affect_On_Output, Value_Missing_Affect_On_Input, Str_Missing_Affect_On_Output,
     "Affection", "Output",        fix_Value_Missing_Affect_On_Output, 0, 0},

    {0, 0, 0, 
     0, 0, (void (*)(analysis_type *analysis, struct itl *itp, const char *value)) 0, 0}
};


/* *******************************************************************************
                               Interactive Menu
******************************************************************************* */

static void Value_Missing_menu(analysis_type *analysis)
{
//  return;
    int all_i = 1, qin_i = 2, qout_i = 3, ain_i = 4, aout_i = 5;
    int ans, idx, not_done = 1;
    int ttl;
    int trans[6];
    const char *show;

    fix_Value_Missing_default(analysis, 0);
    fix_Value_Missing_default(analysis, 1);
    fix_Value_Missing_default(analysis, 2);
    fix_Value_Missing_default(analysis, 3);
    while (not_done) {
        ttl = 0;
        fflush(stdout);
        draw_line();
        printf("              Mega2 %s Missing Value menu:\n", Mega2Version);
        draw_line();
        printf("If it is necessary specify the value to indicate that a trait is missing,\n");
        printf("both for input to Mega2 and/or output from Mega2.\n");
        printf("Note: Output entries that are marked with a \"#\" can not be changed.\n\n");
        printf("0) Done with this menu - please proceed\n");
        if (1) {
            printf("%2d) Specify default for ALL missing traits:                         %s\n",
                   ++ttl, "");
            trans[ttl] = all_i;
        }

        if (1) {
            printf("%2d) Specify missing value for Quantitative traits read:    %7s [%s]\n",
                   ++ttl, Value_Missing[qin_i-2].str, source_name[Value_Missing[qin_i-2].source]);
            trans[ttl] = qin_i;
        }

        if (1) {
            printf("%2d) Specify missing value for Affection status read:       %7s [%s]\n",
                   ++ttl, Value_Missing[ain_i-2].str, source_name[Value_Missing[ain_i-2].source]);
            trans[ttl] = ain_i;
        }

        show = Value_Missing[qout_i-2].str;
        if (!*show) show = (*analysis)->output_quant_default_value();
        if (!show) show = "";
        if (fix_Value_Missing_check_allow(analysis, 1)) {
            printf("%2d) Specify missing value for Quantitative traits written: %7s [%s]\n",
                   ++ttl, show, source_name[Value_Missing[qout_i-2].source]);
            trans[ttl] = qout_i;
        } else {
            Value_Missing[qout_i-2].source = 6;
            printf(" #) Fixed missing value for Quantitative traits written:   %7s [%s]\n",
                   show, source_name[Value_Missing[qout_i-2].source]);
        }

        show = Value_Missing[aout_i-2].str;
        if (!*show) show = (*analysis)->output_affect_default_value();
        if (!show) show = "";
        if (fix_Value_Missing_check_allow(analysis, 3)) {
            printf("%2d) Specify missing value for Affection status written:    %7s [%s]\n",
                   ++ttl, show, source_name[Value_Missing[aout_i-2].source]);
            trans[ttl] = aout_i;
        } else {
            Value_Missing[aout_i-2].source = 6;
            printf(" #) Fixed missing value for Affection status written:      %7s [%s]\n",
                   show, source_name[Value_Missing[aout_i-2].source]);
        }

        printf("Select from options 0-%d > ", ttl);

        fcmap(stdin, "%d", &ans); newline;
        if (ans == 0) {
//          asm("int $3");
            not_done = fix_Value_Missing_all(analysis, 0);
            continue;
        } else if (ans > ttl || ans < 1) {
            printf("The allowed values are:");
            for (idx = 0; idx <= ttl; idx++)
                printf(" %d%s", idx, idx == ttl? "." : ",");
            printf("\n");
            continue;
        }
        idx = trans[ans];
        if (idx != all_i) {
            printf("Please enter Missing %s %s Value (or enter \"clear\") > ",
                   Value_Missing[idx - 2].name, Value_Missing[idx - 2].put);
            Value_Missing[idx - 2].set = 1;
            fcmap(stdin, "%s", Str_read);
            if (!strcmp(Str_read, "clear")) {
                Value_Missing[idx - 2].str[0] = 0;
                Value_Missing[idx - 2].set = 0;
            } else if (! fix_Value_Missing_check_numeric(analysis, &Value_Missing[idx-2], Str_read)) {
                strcpy(Value_Missing[idx - 2].str, Str_read);
                Value_Missing[idx-2].source = 0;
            }
        } else {
            printf("Please enter Default Missing Value (or enter \"clear\") > ");
            fcmap(stdin, "%s", Str_read);
            Value_Missing[0].set = 1;

            if (! fix_Value_Missing_check_numeric(analysis, &Value_Missing[0], Str_read)) {
                strcpy(Value_Missing[0].str, Str_read);
                Value_Missing[0].source = 0;
            }

            if (!strcmp(Value_Missing[0].str, "clear")) {
                Value_Missing[0].str[0] = 0;
                Value_Missing[0].set = 0;
            }

            if (fix_Value_Missing_check_allow(analysis, 1)) {
                if (! fix_Value_Missing_check_numeric(analysis, &Value_Missing[1], Str_read)) {
                    strcpy(Value_Missing[1].str, Str_read);
                    Value_Missing[1].set = Value_Missing[0].set;
                    Value_Missing[1].source = 0;
                }
            }

            if (! fix_Value_Missing_check_numeric(analysis, &Value_Missing[2], Str_read)) {
                strcpy(Value_Missing[2].str, Str_read);
                Value_Missing[2].set = Value_Missing[0].set;
                Value_Missing[2].source = 0;
            }

            if (fix_Value_Missing_check_allow(analysis, 3)) {
                if (! fix_Value_Missing_check_numeric(analysis, &Value_Missing[3], Str_read)) {
                    strcpy(Value_Missing[3].str, Str_read);
                    Value_Missing[3].set = Value_Missing[0].set;
                    Value_Missing[3].source = 0;
                }
            }
        }
    }
}

/* *******************************************************************************
                                   BATCH FILE
******************************************************************************* */

void fix_Value_Missing_inherit_Quant_2_Affect() {
    Value_Missing[2  /*Value_Missing[_Affect_On_Input]*/].inherit = Value_Missing_Quant_On_Input;
    Value_Missing[3 /*Value_Missing[_Affect_On_Output]*/].inherit = Value_Missing_Quant_On_Output;
}

static int fix_Value_Missing_check_allow(analysis_type *analysis, int vmidx) {
    struct itl *itp = &Value_Missing[vmidx];
    int allow = 1;
    if (itp->it == Value_Missing_Quant_On_Input)
        allow = 1;
    else if (itp->it == Value_Missing_Affect_On_Input)
        allow = 1;
    else if (itp->it == Value_Missing_Quant_On_Output)
        allow = (*analysis)->output_quant_can_define_missing_value();
    else if (itp->it == Value_Missing_Affect_On_Output)
        allow = (*analysis)->output_affect_can_define_missing_value();
    return allow;
}

static int fix_Value_Missing_check_numeric(analysis_type *analysis, struct itl *itp, const char *value) {
    int qnum = 0;
    int num = 0;
    if (itp->it == Value_Missing_Quant_On_Input) {
        if (Input_Format != in_format_mega2 || strcasecmp(value, "na"))
            qnum = 1;
    } else if (itp->it == Value_Missing_Affect_On_Input) {
        if (Input_Format != in_format_mega2 || strcasecmp(value, "na")) {
            if (PLINK.plink || PLINK.xcf)
                qnum = 1;
            else
                num = 1;
        }
    } else if (itp->it == Value_Missing_Quant_On_Output) {
        qnum = (*analysis)->output_quant_must_be_numeric();
    } else if (itp->it == Value_Missing_Affect_On_Output) {
        if ((*analysis)->output_affect_must_be_numeric()) {
            if (*analysis == TO_PLINK)
                qnum = 1;
            else
                num = 1;
        }
    }

    // Since the user did not specify the default value, this is an internal check.
    if (num || qnum) {
        char *end;
        if (qnum)
            (void) strtod(value, &end);
        else if (num)
            (void) strtol(value, &end, 10);
        if (strlen(value) == 0 || strlen(end) != 0 || errno == ERANGE) {
            // Conversion of the entire string was not successful or some other error...
            errorvf("%s %s missing value \"%s\" must be a numeric string representing %s number.\n",
                    itp->name, itp->put, value, qnum ? "a float" : "an integer");
            return 1;
        }
    }
    return 0;
}

/*
            switch (Input_Format) {
            case in_format_mega2:
                break;
            case in_format_linkage:
            case in_format_extended_linkage:
                break;
            case in_format_binary_PED:
            case in_format_PED:
            case in_format_binary_VCF:
            case in_format_compressed_VCF:
            case in_format_VCF:
                break;
            default:
                break;
            }
*/

static void fix_Value_Missing_default(analysis_type *analysis, int vmidx)
{
    int ret = 0;
    struct itl *itp = &Value_Missing[vmidx];
    const char *df = 0;
    int allow;
    char dfp[32];

    allow = fix_Value_Missing_check_allow(analysis, vmidx);

    if (itp->it == Value_Missing_Quant_On_Input) {
        switch (Input_Format) {
        case in_format_mega2:
            itp->set = 1;
            itp->source = 9;
            df = "NA";
            break;
        case in_format_linkage:
        case in_format_extended_linkage:
            itp->set = 1;
            itp->source = 8;
            df = "-999";
            break;
        case in_format_binary_PED:
        case in_format_PED:
        case in_format_binary_VCF:
        case in_format_compressed_VCF:
        case in_format_VCF:
            itp->set = 1;
            if (PLINK.missing_pheno) {
                itp->source = 2;
                sprintf(dfp, "%g", PLINK.pheno_value);
                df = dfp;
            } else {
                itp->source = 3;
                df = "-9";
            }
            break;
        default:
            break;
        }
    } else if (itp->it == Value_Missing_Affect_On_Input)
        switch (Input_Format) {
        case in_format_mega2:
            itp->set = 1;
            itp->source = 9;
            df = "NA";
            break;
        case in_format_linkage:
        case in_format_extended_linkage:
            itp->set = 1;
            itp->source = 8;
            df = "0";
            break;
        case in_format_binary_PED:
        case in_format_PED:
        case in_format_binary_VCF:
        case in_format_compressed_VCF:
        case in_format_VCF:
            itp->set = 1;
            if (PLINK.missing_pheno) {
                itp->source = 2;
                sprintf(dfp, "%g", PLINK.pheno_value);
                df = dfp;
            } else {
                itp->source = 3;
                df = "-9";
            }
            break;
        default:
            break;
        }
    else if (itp->it == Value_Missing_Quant_On_Output) {
        itp->set = allow;
        df = (*analysis)->output_quant_default_value();
        if (df == (const char *)NULL) {
            itp->source = 7;
            df = Mega2BatchItems[itp->it].value.name;
        } else
            itp->source = 4;
    } else if (itp->it == Value_Missing_Affect_On_Output) {
        itp->set = allow;
        df = (*analysis)->output_affect_default_value();
        if (df == (const char *)NULL) {
            itp->source = 7;
            df = Mega2BatchItems[itp->it].value.name;
        } else
            itp->source = 4;
    }

    if (df != (const char *)NULL) {
        ret = fix_Value_Missing_check_numeric(analysis, itp, df);
        if (! ret) {
//              itp->store(analysis, itp, df);
            strcpy(itp->str, df);
        }
    }
}

static int fix_Value_Missing(analysis_type *analysis, int vmidx, int init_str)
{
    int ret = 0;
    struct itl *itp = &Value_Missing[vmidx];
    int allow = fix_Value_Missing_check_allow(analysis, vmidx);
//  if (ITEM_READ(itp->it))
    if (itp->set)
    {
        if (! allow) {
            warnvf("'Analysis_Option' = '%s' does not allow the definition of Missing %s %s Value.\n",
                   ProgName, itp->name, itp->put);
            return 0; // error; but not before
        }
        ret = fix_Value_Missing_check_numeric(analysis, itp, itp->str);
        if (! ret &&  !init_str) {
            itp->store(analysis, itp, itp->str);
            if (! batchINPUTFILES) {
                Mega2BatchItems[itp->it].item_read = 1;
                batchf(itp->it);
            }
        }
        return ret;
    } else {
        const char *df = 0;
        char dfp[32];
        if (itp->it == Value_Missing_Quant_On_Input) {
            if (PLINK.plink || PLINK.xcf) {
                if (PLINK.missing_pheno) {
                    if (!init_str) itp->source = 2;
                    sprintf(dfp, "%g", PLINK.pheno_value);
                    df = dfp;
                } else {
                    if (!init_str) itp->source = 3;
                    df = "-9";
                }
            } else
                df = 0;
        } else if (itp->it == Value_Missing_Affect_On_Input)
/*
            if (PLINK.plink || PLINK.xcf) {
                if (PLINK.missing_pheno) {
                    sprintf(dfp, "%g", PLINK.pheno_value);
                    df = dfp;
                    if (!init_str) itp->source = 2;
                } else {
                    df = "-9";
                    if (!init_str) itp->source = 3;
                }
            } else
*/
                df = 0;
        else if (itp->it == Value_Missing_Quant_On_Output) {
            df = (*analysis)->output_quant_default_value();
            if (!init_str && df) itp->source = 4;
        } else if (itp->it == Value_Missing_Affect_On_Output) {
            df = (*analysis)->output_affect_default_value();
            if (!init_str && df) itp->source = 4;
        }
 
        if (df != (const char *)NULL) {
            ret = fix_Value_Missing_check_numeric(analysis, itp, df);
            if (! ret) {
                if (!init_str)
                    itp->store(analysis, itp, df);
                strcpy(itp->str, df);
            }
        } else if (allow) {
// Note:
//            Mega2BatchItems[itp->it].item_read = 1;
// is not set.  Because we are guessing value.
            if (!init_str) itp->source = 5;
            int inherit = itp->inherit;
            if (inherit) {
                switch(itp->it) {
                case Value_Missing_Affect_On_Input:
                    switch(inherit) {
                    case Value_Missing_Quant_On_Input:
                        if (!init_str)
                            strcpy(Mega2BatchItems[itp->it].value.name,
                                   Value_Missing[0 /*Value_Missing[_Quant_On_Input]*/].str);
                        strcpy(Value_Missing[2 /*Value_Missing[_Affect_On_Input]*/].str,
                               Value_Missing[0 /*Value_Missing[_Quant_On_Input]*/].str);
                        break;
                    default:
                        break;
                    }
                    break;
                case Value_Missing_Quant_On_Output:
                    switch(inherit) {
                    case Value_Missing_Quant_On_Input:
                        if (!init_str)
                            strcpy(Mega2BatchItems[itp->it].value.name,
                                   Value_Missing[0 /*Value_Missing[_Quant_On_Input]*/].str);
                        strcpy(Value_Missing[1 /*Value_Missing[_Quant_On_Output]*/].str,
                               Value_Missing[0 /*Value_Missing[_Quant_On_Input]*/].str);
                        break;
                    default:
                        break;
                    }
                    break;
                case Value_Missing_Affect_On_Output:
                    switch(inherit) {
                    case Value_Missing_Quant_On_Input:
                        if (!init_str)
                            strcpy(Mega2BatchItems[itp->it].value.name,
                                   Value_Missing[0 /*Value_Missing[_Quant_On_Input]*/].str);
                        strcpy(Value_Missing[3 /*Value_Missing[_Quant_On_Output]*/].str,
                               Value_Missing[0 /*Value_Missing[_Quant_On_Input]*/].str);

                        break;
                    case Value_Missing_Quant_On_Output:
                        if (!init_str)
                            strcpy(Mega2BatchItems[itp->it].value.name,
                                   Value_Missing[1 /*Value_Missing[_Quant_On_Output]*/].str);
                        strcpy(Value_Missing[3 /*Value_Missing[_Quant_On_Output]*/].str,
                               Value_Missing[1 /*Value_Missing[_Quant_On_Input]*/].str);
                        break;
                    case Value_Missing_Affect_On_Input:
                        if (!init_str)
                            strcpy(Mega2BatchItems[itp->it].value.name,
                                   Value_Missing[2 /*Value_Missing[_Affect_On_Input]*/].str);
                        strcpy(Value_Missing[3 /*Value_Missing[_Quant_On_Output]*/].str,
                               Value_Missing[2 /*Value_Missing[_Quant_On_Input]*/].str);
                        break;
                    default:
                        break;
                    }
                    break;
                default:
                    break;
                }
            }
        } else {
            if (!init_str)  itp->source = 6;
        }
    }
    return 0;
}

static void fix_Value_Missing_Quant_On_Input(analysis_type *analysis, struct itl *itp, const char *value) {
    // QMISSING is the internal numeric value of NA.
    if (strcasecmp(value, "NA") == 0) {
        Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].value.fvalue = QMISSING;
    } else {
        char *endptr;
        // Interesting, but MissingQuant was not assigned here in the past!
        // NOTE: the use of 'strtok' here assumes that there is no comment following the value
        MissingQuant = strtod(value, &endptr);
        Mega2BatchItems[/* 17 */ Value_Missing_Quant_On_Input].value.fvalue = MissingQuant;
    }
}

static void fix_Value_Missing_Quant_On_Output(analysis_type *analysis, struct itl *itp, const char *value) {
    strcpy(Mega2BatchItems[/* 49 */ Value_Missing_Quant_On_Output].value.name, value);
}

static void fix_Value_Missing_Affect_On_Input(analysis_type *analysis, struct itl *itp, const char *value) {
//    if (strcasecmp(value, "NA") == 0)
//        Mega2BatchItems[itp->it].value.option = 0;
//    else
        strcpy(Mega2BatchItems[itp->it].value.name, value);
}

static void fix_Value_Missing_Affect_On_Output(analysis_type *analysis, struct itl *itp, const char *value) {
    strcpy(Mega2BatchItems[/* 49 */ Value_Missing_Affect_On_Output].value.name, value);
}

static int fix_Value_Missing_all(analysis_type *analysis, int init_str) {
    int ret = 0;
    ret = fix_Value_Missing(analysis, 0, init_str);
    if (ret) return 1;

    ret = fix_Value_Missing(analysis, 1, init_str);
    if (ret) return 1;

    ret = fix_Value_Missing(analysis, 2, init_str);
    if (ret) return 1;

    ret = fix_Value_Missing(analysis, 3, init_str);
    if (ret) return 1;
    return ret;
}

/* *******************************************************************************
                                   Entry Switch
******************************************************************************* */

void Value_Missing_get(analysis_type *analysis)
{
    struct itl *itp = &Value_Missing[0], *itpn;

    if (batchINPUTFILES) {
        for (int i = 0; i < 4; i++) {
            itpn = itp + i;
            itpn->source = 1;
            itpn->set = Mega2BatchItems[itpn->it].item_read;
        }

        if (fix_Value_Missing_all(analysis, 0))
            EXIT(DATA_INCONSISTENCY);
    } else {
        Value_Missing_menu(analysis);
    }

    itpn = itp + 0;
    if (Mega2BatchItems[itpn->it].value.fvalue == QMISSING)
        msgvf("%s  %s Missing Value \"%s\" [\"%s\" %s]\n",
              itpn->name, itpn->put, "AN", itpn->str, 
              source_name[itpn->source]);
    else
        msgvf("%s  %s Missing Value  %g  [\"%s\" %s]\n",
              itpn->name, itpn->put, Mega2BatchItems[itpn->it].value.fvalue, itpn->str, 
              source_name[itpn->source]);

    itpn = itp + 2;
    msgvf("%s     %s Missing Value \"%s\" [\"%s\" %s]\n",
          itpn->name, itpn->put, Mega2BatchItems[itpn->it].value.name, itpn->str, 
          source_name[itpn->source]);

    itpn = itp + 1;
    msgvf("%s %s Missing Value \"%s\" [\"%s\" %s]\n",
          itpn->name, itpn->put, Mega2BatchItems[itpn->it].value.name, itpn->str, 
          source_name[itpn->source]);

    itpn = itp + 3;
    msgvf("%s    %s Missing Value \"%s\" [\"%s\" %s]\n",
          itpn->name, itpn->put, Mega2BatchItems[itpn->it].value.name, itpn->str, 
          source_name[itpn->source]);

}
