/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 2012-2014 Robert Baron, Charles P. Kollar,
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

#ifndef FILE_OPS_H
#define FILE_OPS_H

#include <stdio.h>
#include <string.h>
#include <zlib.h>
#include "common.h"
#include "typedefs.h"

#include "error_messages_ext.h"
#include "utils_ext.h"

class file_ops {
protected:
    char path_[2*FILENAME_LENGTH];
    const char *dir_;
    const char *file_;
public:
    char **_opath;
    FILE *_filep;
    // In Version 1.2.6 gzFile is (gzFile_s *), and in 1.2.5 it is (voidp) which is a (Byte *).
    // So, it really is a pointer and not an inline structure, and so we can compare it to NULL.
    // Also note that the documentation for 'gzFile gzopen()' states that it returns NULL on failure.
    gzFile _gzfile;

    file_ops() : _filep(NULL), _gzfile(NULL)  { }

    virtual ~file_ops() {
        // It is not a problem to try to close the file multiple times...
        filep_close();
    }
    
    /**
       @brief Write a file header
     
       The user will specialize this method in order to write a file header
       to the newly opened file @see filep_open
       @return void
    */
    virtual void file_header() {};
    
    /**
       @brief Open the file for write or append
     
       The file will be written to the 'file' name and 'directory' provided.
     
       If the extension of the file is '.gz' it will be processed as a gzip file
       if gzip is supported on this platform, otherwise an error will be issued.
     
       @param[in] dir       string associated with directory prefix
       @param[in] file      string associated with the file name
       @param[in] mode      string "w"rite, "a"ppend
    */
    void filep_open(const char *dir, const char *file, const char *mode) {
        sprintf(path_, "%s/%s", dir, file);
        int _is_gzfile = (strlen(file) > 3 && strcasecmp(".gz", &file[strlen(file)-3]) == 0) ? 1 : 0;

        // The current implementation of this object is for output only, so make sure that we are never reading...
        if ((strchr(mode, 'w') == (char *)NULL) && (strchr(mode, 'a') == (char *)NULL)) {
            draw_line();
            errorvf("INTERNAL: Must open file '%s' for write/append only.\n", path_);
            EXIT(FILE_WRITE_ERROR);
        }
        
        this->dir_   = dir;
        this->file_ = &path_[strlen(dir)+1];
        
        // Open the file...
        if (_is_gzfile) {
            char mode_gz[FILENAME_LENGTH], *mode_ptr = (char *)mode;
            // mode must be "wb" or "ab" for gzipped file...
            if (strchr(mode, 'b') == (char *)NULL) {
                strcpy(mode_gz, mode);
                strcat(mode_gz, "b");
                mode_ptr = mode_gz;
            }
            // Returns NULL on failure...
            _gzfile = gzopen(path_, mode_ptr);
        } else
            _filep = fopen(path_, mode);

        // Was there an error???
        if ((_is_gzfile ? (void *)_gzfile : (void *)_filep) == (void *)NULL) {
            draw_line();
            errorvf("Unable to open file '%s' for writing\n", path_);
            EXIT(FILE_WRITE_ERROR);
        }
        
        // @see file_header
        file_header();
    }
    
    const char *file_path() { return path_; }

    /**
       @brief Print a formatted string to the file.

       The file must have been opened (@see filep_open) must be called before calling this method.
     */
    void pr_printf(const char *fmt, ...) {
        va_list ap;
        va_start(ap, fmt);
        if ((void *)_gzfile != (void *)NULL) {
#define BUFFER_SIZE            20480
            char buffer[BUFFER_SIZE];
            int rtn;
            rtn = vsnprintf(buffer, BUFFER_SIZE-1, fmt, ap);
            if (rtn >= BUFFER_SIZE) {
                errorvf("INTERNAL: Buffer too small for zipped file processing '%s'\n", path_);
                EXIT(FILE_WRITE_ERROR);
            }
            gzputs(_gzfile, buffer);
        } else
            vfprintf(_filep, fmt, ap);
        va_end(ap);
    }

    /**
       @brief Write a null terminated string to the file.

       The file must have been opened (@see filep_open) must be called before calling this method.
     */
    void pr_puts(const char *s) {
        if ((void *)_gzfile != (void *)NULL)
            gzputs(_gzfile, s);
        else
            fputs(s, _filep);
    }

    /**
       @brief Write a new line character to the file.

       NOTE: We should consider doing the right thing when writing a file in windows (e.g., "\r\n").
    */
    void pr_nl() { pr_printf("\n"); }

    void file_flush() {
      if ((void *)_gzfile != (void *)NULL)
	gzflush(_gzfile, 0);
      else
       if (_filep != (FILE *)NULL)
	fflush(_filep);
    }

    /**
       @brief Write to file before it is closed

       @see filep_close
    */
    virtual void file_trailer() {};

    /**
       @brief Perform actions on file after it exists.

       @see filep_close
     */
    virtual void file_post() {};

    /**
       @brief Close the file

       Close the file. Use @see file_trailer to write information to the end of
       the file before closing it. Use @see file_post to modify aspects of the file
       after it has been closed (e.g., permissions).

       It is not an error to close a file multiple times, but bad form.
    */
    void filep_close() {
        if ((void *)_gzfile != (void *)NULL) {
            file_trailer();
            gzclose(_gzfile);
            _gzfile = NULL;
            file_post();
        } else if (_filep != (FILE *)NULL) {
            file_trailer();
            fclose(_filep);
            _filep = NULL;
            file_post();
        }
    }

    /**
       @brief Delete the file
     */
    void filep_delete(const char *dir, const char *file) {
        sprintf(path_, "%s/%s", dir, file);
        delete_file(path_);
    }
};

#endif /* FILE_OPS_H */
