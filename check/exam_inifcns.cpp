/** @file exam_inifcns.cpp
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

#include "exams.h"

/* Assorted tests on other transcendental functions. */
static unsigned inifcns_consist_trans(void)
{
    unsigned result = 0;
    symbol x("x");
    ex chk;
    
    chk = asin(1)-acos(0);
    if (!chk.is_zero()) {
        clog << "asin(1)-acos(0) erroneously returned " << chk
             << " instead of 0" << endl;
        ++result;
    }
    
    // arbitrary check of type sin(f(x)):
    chk = pow(sin(acos(x)),2) + pow(sin(asin(x)),2)
        - (1+pow(x,2))*pow(sin(atan(x)),2);
    if (chk != 1-pow(x,2)) {
        clog << "sin(acos(x))^2 + sin(asin(x))^2 - (1+x^2)*sin(atan(x))^2 "
             << "erroneously returned " << chk << " instead of 1-x^2" << endl;
        ++result;
    }
    
    // arbitrary check of type cos(f(x)):
    chk = pow(cos(acos(x)),2) + pow(cos(asin(x)),2)
        - (1+pow(x,2))*pow(cos(atan(x)),2);
    if (!chk.is_zero()) {
        clog << "cos(acos(x))^2 + cos(asin(x))^2 - (1+x^2)*cos(atan(x))^2 "
             << "erroneously returned " << chk << " instead of 0" << endl;
        ++result;
    }
    
    // arbitrary check of type tan(f(x)):
    chk = tan(acos(x))*tan(asin(x)) - tan(atan(x));
    if (chk != 1-x) {
        clog << "tan(acos(x))*tan(asin(x)) - tan(atan(x)) "
             << "erroneously returned " << chk << " instead of -x+1" << endl;
        ++result;
    }
    
    // arbitrary check of type sinh(f(x)):
    chk = -pow(sinh(acosh(x)),2).expand()*pow(sinh(atanh(x)),2)
        - pow(sinh(asinh(x)),2);
    if (!chk.is_zero()) {
        clog << "expand(-(sinh(acosh(x)))^2)*(sinh(atanh(x))^2) - sinh(asinh(x))^2 "
             << "erroneously returned " << chk << " instead of 0" << endl;
        ++result;
    }
    
    // arbitrary check of type cosh(f(x)):
    chk = (pow(cosh(asinh(x)),2) - 2*pow(cosh(acosh(x)),2))
        * pow(cosh(atanh(x)),2);
    if (chk != 1) {
        clog << "(cosh(asinh(x))^2 - 2*cosh(acosh(x))^2) * cosh(atanh(x))^2 "
             << "erroneously returned " << chk << " instead of 1" << endl;
        ++result;
    }
    
    // arbitrary check of type tanh(f(x)):
    chk = (pow(tanh(asinh(x)),-2) - pow(tanh(acosh(x)),2)).expand()
        * pow(tanh(atanh(x)),2);
    if (chk != 2) {
        clog << "expand(tanh(acosh(x))^2 - tanh(asinh(x))^(-2)) * tanh(atanh(x))^2 "
             << "erroneously returned " << chk << " instead of 2" << endl;
        ++result;
    }
    
    return result;
}

/* Simple tests on the tgamma function.  We stuff in arguments where the results
 * exists in closed form and check if it's ok. */
static unsigned inifcns_consist_gamma(void)
{
    unsigned result = 0;
    ex e;
    
    e = tgamma(ex(1));
    for (int i=2; i<8; ++i)
        e += tgamma(ex(i));
    if (e != numeric(874)) {
        clog << "tgamma(1)+...+tgamma(7) erroneously returned "
             << e << " instead of 874" << endl;
        ++result;
    }
    
    e = tgamma(ex(1));
    for (int i=2; i<8; ++i)
        e *= tgamma(ex(i));    
    if (e != numeric(24883200)) {
        clog << "tgamma(1)*...*tgamma(7) erroneously returned "
             << e << " instead of 24883200" << endl;
        ++result;
    }
    
    e = tgamma(ex(numeric(5, 2)))*tgamma(ex(numeric(9, 2)))*64;
    if (e != 315*Pi) {
        clog << "64*tgamma(5/2)*tgamma(9/2) erroneously returned "
             << e << " instead of 315*Pi" << endl;
        ++result;
    }
    
    e = tgamma(ex(numeric(-13, 2)));
    for (int i=-13; i<7; i=i+2)
        e += tgamma(ex(numeric(i, 2)));
    e = (e*tgamma(ex(numeric(15, 2)))*numeric(512));
    if (e != numeric(633935)*Pi) {
        clog << "512*(tgamma(-13/2)+...+tgamma(5/2))*tgamma(15/2) erroneously returned "
             << e << " instead of 633935*Pi" << endl;
        ++result;
    }
    
    return result;
}

/* Simple tests on the Psi-function (aka polygamma-function).  We stuff in
   arguments where the result exists in closed form and check if it's ok. */
static unsigned inifcns_consist_psi(void)
{
    unsigned result = 0;
    symbol x;
    ex e, f;
    
    // We check psi(1) and psi(1/2) implicitly by calculating the curious
    // little identity tgamma(1)'/tgamma(1) - tgamma(1/2)'/tgamma(1/2) == 2*log(2).
    e += (tgamma(x).diff(x)/tgamma(x)).subs(x==numeric(1));
    e -= (tgamma(x).diff(x)/tgamma(x)).subs(x==numeric(1,2));
    if (e!=2*log(2)) {
        clog << "tgamma(1)'/tgamma(1) - tgamma(1/2)'/tgamma(1/2) erroneously returned "
             << e << " instead of 2*log(2)" << endl;
        ++result;
    }
    
    return result;
}

/* Simple tests on the Riemann Zeta function.  We stuff in arguments where the
 * result exists in closed form and check if it's ok.  Of course, this checks
 * the Bernoulli numbers as a side effect. */
static unsigned inifcns_consist_zeta(void)
{
    unsigned result = 0;
    ex e;
    
    for (int i=0; i<13; i+=2)
        e += zeta(i)/pow(Pi,i);
    if (e!=numeric(-204992279,638512875)) {
        clog << "zeta(0) + zeta(2) + ... + zeta(12) erroneously returned "
             << e << " instead of -204992279/638512875" << endl;
        ++result;
    }
    
    e = 0;
    for (int i=-1; i>-16; i--)
        e += zeta(i);
    if (e!=numeric(487871,1633632)) {
        clog << "zeta(-1) + zeta(-2) + ... + zeta(-15) erroneously returned "
             << e << " instead of 487871/1633632" << endl;
        ++result;
    }
    
    return result;
}

unsigned exam_inifcns(void)
{
    unsigned result = 0;
    
    cout << "examining consistency of symbolic functions" << flush;
    clog << "----------consistency of symbolic functions:" << endl;
    
    result += inifcns_consist_trans();  cout << '.' << flush;
    result += inifcns_consist_gamma();  cout << '.' << flush;
    result += inifcns_consist_psi();  cout << '.' << flush;
    result += inifcns_consist_zeta();  cout << '.' << flush;
    
    if (!result) {
        cout << " passed " << endl;
        clog << "(no output)" << endl;
    } else {
        cout << " failed " << endl;
    }
    
    return result;
}