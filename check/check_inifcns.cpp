/** @file check_inifcns.cpp
 *
 *  This test routine applies assorted tests on initially known higher level
 *  functions. */

/*
 *  GiNaC Copyright (C) 1999-2000 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "checks.h"

/* Some tests on the sine trigonometric function. */
static unsigned inifcns_check_sin(void)
{
    unsigned result = 0;
    bool errorflag = false;
    
    // sin(n*Pi) == 0?
    errorflag = false;
    for (int n=-10; n<=10; ++n) {
        if (sin(n*Pi).eval() != numeric(0) ||
            !sin(n*Pi).eval().info(info_flags::integer))
            errorflag = true;
    }
    if (errorflag) {
        // we don't count each of those errors
        clog << "sin(n*Pi) with integer n does not always return exact 0"
             << endl;
        ++result;
    }
    
    // sin((n+1/2)*Pi) == {+|-}1?
    errorflag = false;
    for (int n=-10; n<=10; ++n) {
        if (!sin((n+numeric(1,2))*Pi).eval().info(info_flags::integer) ||
            !(sin((n+numeric(1,2))*Pi).eval() == numeric(1) ||
              sin((n+numeric(1,2))*Pi).eval() == numeric(-1)))
            errorflag = true;
    }
    if (errorflag) {
        clog << "sin((n+1/2)*Pi) with integer n does not always return exact {+|-}1"
             << endl;
        ++result;
    }
    
    // compare sin((q*Pi).evalf()) with sin(q*Pi).eval().evalf() at various
    // points.  E.g. if sin(Pi/10) returns something symbolic this should be
    // equal to sqrt(5)/4-1/4.  This routine will spot programming mistakes
    // of this kind:
    errorflag = false;
    ex argument;
    numeric epsilon(double(1e-8));
    for (int n=-340; n<=340; ++n) {
        argument = n*Pi/60;
        if (abs(sin(evalf(argument))-evalf(sin(argument)))>epsilon) {
            clog << "sin(" << argument << ") returns "
                 << sin(argument) << endl;
            errorflag = true;
        }
    }
    if (errorflag) {
        ++result;
    }
    
    return result;
}

/* Simple tests on the cosine trigonometric function. */
static unsigned inifcns_check_cos(void)
{
    unsigned result = 0;
    bool errorflag;
    
    // cos((n+1/2)*Pi) == 0?
    errorflag = false;
    for (int n=-10; n<=10; ++n) {
        if (cos((n+numeric(1,2))*Pi).eval() != numeric(0) ||
            !cos((n+numeric(1,2))*Pi).eval().info(info_flags::integer))
            errorflag = true;
    }
    if (errorflag) {
        clog << "cos((n+1/2)*Pi) with integer n does not always return exact 0"
             << endl;
        ++result;
    }
    
    // cos(n*Pi) == 0?
    errorflag = false;
    for (int n=-10; n<=10; ++n) {
        if (!cos(n*Pi).eval().info(info_flags::integer) ||
            !(cos(n*Pi).eval() == numeric(1) ||
              cos(n*Pi).eval() == numeric(-1)))
            errorflag = true;
    }
    if (errorflag) {
        clog << "cos(n*Pi) with integer n does not always return exact {+|-}1"
             << endl;
        ++result;
    }
    
    // compare cos((q*Pi).evalf()) with cos(q*Pi).eval().evalf() at various
    // points.  E.g. if cos(Pi/12) returns something symbolic this should be
    // equal to 1/4*(1+1/3*sqrt(3))*sqrt(6).  This routine will spot
    // programming mistakes of this kind:
    errorflag = false;
    ex argument;
    numeric epsilon(double(1e-8));
    for (int n=-340; n<=340; ++n) {
        argument = n*Pi/60;
        if (abs(cos(evalf(argument))-evalf(cos(argument)))>epsilon) {
            clog << "cos(" << argument << ") returns "
                 << cos(argument) << endl;
            errorflag = true;
        }
    }
    if (errorflag) {
        ++result;
    }
    
    return result;
}

/* Simple tests on the tangent trigonometric function. */
static unsigned inifcns_check_tan(void)
{
    unsigned result = 0;
    bool errorflag;
    
    // compare tan((q*Pi).evalf()) with tan(q*Pi).eval().evalf() at various
    // points.  E.g. if tan(Pi/12) returns something symbolic this should be
    // equal to 2-sqrt(3).  This routine will spot programming mistakes of 
    // this kind:
    errorflag = false;
    ex argument;
    numeric epsilon(double(1e-8));
    for (int n=-340; n<=340; ++n) {
        if (!(n%30) && (n%60))  // skip poles
            ++n;
        argument = n*Pi/60;
        if (abs(tan(evalf(argument))-evalf(tan(argument)))>epsilon) {
            clog << "tan(" << argument << ") returns "
                 << tan(argument) << endl;
            errorflag = true;
        }
    }
    if (errorflag) {
        ++result;
    }
    
    return result;
}

unsigned check_inifcns(void)
{
    unsigned result = 0;

    cout << "checking consistency of symbolic functions" << flush;
    clog << "---------consistency of symbolic functions:" << endl;
    
    result += inifcns_check_sin();  cout << '.' << flush;
    result += inifcns_check_cos();  cout << '.' << flush;
    result += inifcns_check_tan();  cout << '.' << flush;
    
    if (!result) {
        cout << " passed " << endl;
        clog << "(no output)" << endl;
    } else {
        cout << " failed " << endl;
    }
    
    return result;
}