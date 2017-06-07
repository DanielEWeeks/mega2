/*
  Mega2: Manipulation Environment for Genetic Analysis
  Copyright (C) 1999-2017 Robert Baron, Justin R. Stickel, Charles P. Kollar,
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

#ifndef TOD_HH
#define TOD_HH

#include <sys/time.h>
extern void msgvf(const char *fmt, ...);

class Todd {
public:

    Todd():lim(0),str(0) {
        reset();
    }

    Todd(int lim):lim(lim), cnt(0), str(0) {
        reset();
    }

    Todd(const char *str, int lim=0): lim(lim), cnt(0), str(str) {
        reset();
    }

   ~Todd()     { }

    void reset(void) {
        gettimeofday(&tv_base, /* (void *) */ 0);
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
        gettimeofday(&tmp, /* (void *) */ 0);
//      timersub(&tmp, &tv_base, &diff);
        diff.tv_sec  = tmp.tv_sec  - tv_base.tv_sec;
        diff.tv_usec = tmp.tv_usec - tv_base.tv_usec;
        if (tv_base.tv_usec < 0) {
            tv_base.tv_sec--;
            tv_base.tv_usec += 1000000;
        }
        return diff.tv_sec + diff.tv_usec/1000000.0;
    }

};

class Tod: public Todd {
public:
#ifdef TOD
    Tod() : Todd() {}
    Tod(int lim) : Todd(lim) {}
    Tod(const char *str, int lim=0) : Todd(str, lim) {}

#else

    Tod() {}
    Tod(int lim) {}
    Tod(const char *str, int lim=0) {}

    void reset(void) {}
    void operator()(const int xx) {}
    double operator()(void) {return 0.0;}
    void   operator()(const char *) {}
    int ok(void) {return 0;}

#endif

   ~Tod()     { }


};

#endif

