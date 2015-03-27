/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2015 Robert Baron, Charles P. Kollar,
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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02120-1301, USA.

  For further information contact:
      Daniel E. Weeks
      e-mail: weeks@pitt.edu

===========================================================================
*/

#ifndef UTILS_EXT_H
#define UTILS_EXT_H

extern char           sumdir[];

extern void           invalid_value_field(int batch_item);
extern void           unknown_prog(char *prog_name);

extern int            hashtable_strcmp(const void *prev, const void *next);

extern void           Exit(int arg, const char *file, const int line, const char *err);

extern void           backup_file(char *flname);

#ifdef EXPIRE
extern void           check_expiration(void);
#endif

extern void           chmod_X_file(const char *flname);

extern void           chomp(char *str);

extern void           copy_file(char *flname1, char *flname2);

extern void           delete_file(const char *flname);

extern double         dmax(double d1, double d2);

extern void           draw_line(void);

extern void           fcopy(const char *file1, const char *file2, const char *write_mode);

extern void           get_line(FILE *fp, char *retline);

extern int            global_ou_files(int default_opt);

extern void           goodbye(int exit);

extern void           hello(FILE *fp);

extern int            imax(int i1, int i2);

extern int            imin(int i1, int i2);

extern void           init_file(char *file_name);

extern int            is_dir(char *file_name);

extern int            linecount(FILE *fp);

extern int            columncount(FILE *fp);

extern void           log_line(void (*log_func)(const char *messg));

extern void           makedir(char *dirname);

extern void           mega2_opts(int argc, char **argv);

extern const char*    mklogdir(void);

extern void           move_file(char *flname1, char *flname2);

extern void           move_logs(char *new_location);

extern const char*    perl_pgm(const char *pl);

extern int            press_return(void);

extern void           print_mega2_help(void);

extern void           print_mega2_version(void);

extern void           print_only (const char *messg);

extern int            randomindex(int N);

extern double         randomnum(void);

extern void           script_time_stamp(FILE *script_fp);

extern void           seed_random(void);

extern int            show_system_cmd(const char *str, const char *file, const int line);

extern void           summary_time_stamp(char **input_files, FILE *fp,  const char *count_option);

extern char *         strtail(const char *name, const size_t n);

extern void executable_in_path_does_not_exist_csh(FILE *fp, const char *prog, const char *msg);

extern void fprintf_env_checkset_csh (FILE *fp, const char *var, const char *default_value);

extern void fprintf_status_check_csh (FILE *fp, const char *prog, const int exit_on_error);

extern int  have_trait_b (linkage_locus_top *LTop, bool affect = true, bool quant = true);

/*

extern void           TCL_help_information(void);

extern double         choose(int n, int k);

extern double         factorial(int n);

extern int            icomp(int i, int j);

extern void           missing_arg(char *arg);

extern void           missing_param(char *arg);

extern void           string_to_upper(char *strg);

*/
#endif
