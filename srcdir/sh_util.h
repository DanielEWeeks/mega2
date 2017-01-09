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

#ifndef SH_UTIL_H
#define SH_UTIL_H

namespace dataloop {

    class sh_util: public virtual person_locus_entry {
    public:

        sh_util(linkage_ped_top  *Top) : person_locus_entry(Top) { }
        virtual ~sh_util() {}
        void sh_shell_type() {
            pr_printf("#!/bin/csh -f\n");
        }
        void sh_id() {
    #ifndef HIDEFILE
            /* print some identification */
            pr_printf("# C-shell file name: %s\n", path_);
    #endif /* HIDEFILE */
        }
        void sh_rm(const char *filename) {
    #ifdef RM_F
            pr_printf("rm -f %s\n", filename);
    #else /* RM_F */
            pr_printf("if (-e %s) then\n", filename);
            pr_printf("  rm %s\n", filename);
            pr_printf("endif\n");
            pr_printf("\n");
    #endif /* RM_F */
        }
        void sh_ln(const char *target, const char *newname) {
            pr_printf("if (-l %s) then\n", newname);
            pr_printf("  rm %s\n", newname);
            pr_printf("endif\n");
            pr_printf("ln -s %s %s\n", target, newname);
            pr_printf("\n");
        }
        void sh_need_data(const char *pgm, const char *filename) {
            pr_printf("if ((-z %s)) then\n", filename);
            pr_printf("   echo ERROR: %s has no data in %s\n", pgm, filename);
            pr_printf("   exit(-1)\n");
            pr_printf("endif\n");
        }
        void sh_run(const char *pgm, const char *cmd) {
            pr_printf("echo Running the %s program with the command:\n", pgm);
            pr_printf("echo %s", cmd);
            pr_printf("echo\n");
            pr_printf(cmd, "");
        }
        void sh_run(const char *pgm, const char *cmd, const char *out) {
            pr_printf("echo Running the %s program with the command:\n", pgm);
            pr_printf("echo %s \">&\" %s\n", cmd, out);
            pr_printf("echo The output from %s is being saved to %s.\n", pgm, out);
            pr_printf("echo\n");
            pr_printf("%s >& %s\n", cmd, out);
        }
        void sh_echo(const char *cmd) {
            pr_printf("echo %s\n", cmd);
            pr_printf("%s\n", cmd);
        }
        void sh_echo(const char *cmd, const char *out) {
            pr_printf("echo %s \">&\" %s\n", cmd, out);
            pr_printf("%s >& %s\n", cmd, out);
        }
        void sh_cat(const char *cmd, const char *target) {
            pr_printf("if (-e %s.old) then\n", target);
            pr_printf("  echo rm %s.old\n", target);
            pr_printf("  rm %s.old\n", target);
            pr_printf("endif\n");

            pr_printf("if (-e %s) then\n", target);
            pr_printf("  echo mv %s %s.old\n", target, target);
            pr_printf("  mv %s %s.old\n", target, target);
            pr_printf("endif\n");

            pr_printf("echo \"%s > %s\"\n", cmd, target);
            pr_printf("echo \"# DO\" NOT EDIT THIS \"FILE;\" IT IS AUTOMATICALLY GENERATED \"by:\" > %s\n", target);
            pr_printf("echo \"# \" %s >> %s\n", cmd, target);
            pr_printf("echo \"# YOU\" WANT TO EDIT ONE OF THE FILES BEING \"cat'ed\" >> %s\n", target);
            pr_printf("echo >>%s\n", target);
            pr_printf("%s >> %s\n", cmd, target);
        }
        void sh_stdin(const char *pgm, const char *cmd, const char *file) {
            pr_printf("echo Running the %s program with the command:\n", pgm);
            pr_printf("echo %s <%s\n", cmd, file);
            pr_printf("%s <%s\n", cmd, file);
        }
        void sh_inline(const char *pgm, const char *cmd, ...) {
            va_list ap;
            pr_printf("echo Running the %s program with the command:\n", pgm);
            va_start(ap, cmd);
            pr_printf("echo \"       \" \"");
            vfprintf(_filep, cmd, ap);
            pr_printf("\"\n");
            va_end(ap);

            va_start(ap, cmd);
            vfprintf(_filep, cmd, ap);
            fprintf(_filep, "\n");
            va_end(ap);
        }
        void sh_main() {
            sh_shell_type();
            sh_id();
            script_time_stamp(_filep);
            pr_nl();
            pr_printf("set ARGS=");
            pr_nl();
        }

        // NOTE: This is the same code that is found in utils.cpp:fprintf_env_checkset_csh().
        // But because fprintf is coopted, it is not possible to call this routine here...
        void sh_env_checkset_csh(const char *var, const char *default_value) {
            pr_printf("if ( ! $?%s) then\n", var);
            pr_printf("  set %s='%s'\n", var, default_value);
            pr_printf("endif\n");
        }

        void sh_sh(sh_util *sub_shell) {
            pr_nl();
            if (strcmp(sub_shell->dir_, ".")) {
                pr_printf("echo pushd %s\n", sub_shell->dir_);
                pr_printf("pushd %s\n", sub_shell->dir_);
            }

            pr_printf("echo csh %s $ARGS\n", sub_shell->file_);
            pr_printf("csh %s $ARGS\n", sub_shell->file_);

            if (strcmp(sub_shell->dir_, "."))
                pr_printf("popd\n");
            pr_nl();
        }
        void sh_status(const char *pgm, const char *filename) {
            pr_printf("if (!(-e %s)) then\n", filename);
            pr_printf("   echo ERROR: %s failed to create the %s\n", pgm, filename);
            pr_printf("   exit(-1)\n");
            pr_printf("endif\n");
        }
        void sh_save_output(const char *pgm, const char *filename, const char *outfilename) {
            pr_printf("if (-e %s) then\n", filename);
            pr_printf("  mv %s %s\n", filename, outfilename);
            pr_printf("  echo Renamed %s to %s\n", filename, outfilename);
            pr_printf("else\n");
            pr_printf("  echo %s failed to run successfully\n", pgm);
            pr_printf("  exit(-1)\n");
            pr_printf("endif\n");
        }
        void sh_show(const char *pgm, const char *filename) {
            pr_printf("@ xx = -Z %s\n", filename);
            pr_printf("echo stdout for %s in %s is $xx bytes\n", pgm, filename);
            pr_printf("echo\n");
            pr_printf("echo tail -100 %s\n", filename);
            pr_printf("tail -100 %s\n", filename);
        }

        virtual void inner() {}
        virtual void data_loop(const char *dir, const char *fl_name, const char *mode)
        {
            filep_open(dir, fl_name, mode);
            inner();
            filep_close();
        }
    };

    class sh_exec: public sh_util {
    public:
        sh_exec(linkage_ped_top  *Top) : sh_util(Top) { }
        virtual void file_post() {
            chmod_X_file(path_);
        }
    };
};
#endif /* SH_UTIL_H */
