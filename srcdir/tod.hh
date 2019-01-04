/*
  Mega2: Manipulation Environment for Genetic Analysis.

  Copyright 1999-2019, University of Pittsburgh. All Rights Reserved.

  Contributors to Mega2: Robert Baron, Justin R. Stickel, Charles P. Kollar,
  Nandita Mukhopadhyay, Lee Almasy, Mark Schroeder, William P. Mulvihill,
  and Daniel E. Weeks.

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

#ifndef TOD_HH
#define TOD_HH

#ifdef TOD
#include <sys/time.h>
extern void msgvf(const char *fmt, ...);
#endif

class Tod {
public:

#ifndef TOD
    Tod() {}
    Tod(int lim) {}
    Tod(const char *str, int lim=0) {}
    void reset(void) {}
    double operator()(void) {return 0.0;}
    void   operator()(const char *) {}
    int ok(void) {return 0;}

#else

    Tod():lim(0),str(0) {
        reset();
    }

    Tod(int lim):lim(lim), cnt(0), str(0) {
        reset();
    }

    Tod(const char *str, int lim=0): lim(lim), cnt(0), str(str) {
        reset();
    }

   ~Tod()     { }

    void reset(void) {
        gettimeofday(&tv_base, (void *) 0);
    }

    void operator()(const int xx) {
        reset();
    }

    double operator()(void) {
        double td = tdiff();
        if (str != NULL) {
            if (lim == 0 || cnt++ < lim) 
                msgvf("%s %s: %.6f\n", td > 1.0 ? "@@@@": "@@", str, td);
            
            return td;
        }
        else
            return td;
    }

    void operator()(const char *xx) {
        double td = tdiff();
        if (lim == 0 || cnt++ < lim)
            msgvf("%s %s: %.6f\n", td > 30 ? "@@@@": "@@", xx, td);
    }

    int ok(void) {
        return (lim == 0 || cnt < lim) ? 1 : 0;
    }

private:
    int lim, cnt;
    const char *str;
    struct timeval tv_base, tmp, diff;

    double tdiff(void) {
        gettimeofday(&tmp, (void *) 0);
        timersub(&tmp, &tv_base, &diff);
        return diff.tv_sec + diff.tv_usec/1000000.0;
    }

#endif

};

#endif

