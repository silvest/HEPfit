/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gslpp_complex.h>
#include "PVfunctions.h"
#include "LoopToolsWrapper.h"
using namespace gslpp;
using namespace std;

void output(const complex func1, const complex func2) {
    cout << left
         << setw(25) << func1.real()
         << setw(25) << func2.real()
         << setw(25) << func2.real() - func1.real()
         << setw(25) << (func2.real() - func1.real())/(func2.real() + func1.real())
         << endl;
}


int main(int argc, char** argv) {

    cout.precision(16);

    PVfunctions PV(true);
    const double mu2 = 90.0*90.0;
    const double p2  = 95.0*95.0;
    const double m02 = 88.0*88.0;
    const double m12 = 93.0*93.0;
    double Epsilon;

    for (int i=-2; i<=20; i++) {
        Epsilon = pow(10.0, -i);
        cout << left << setw(10) << Epsilon;

        /* logarithm */
        //output(0.0, Epsilon*log(Epsilon));
        //output(1.0/m02, 1.0/(m02*(1.0+Epsilon) - m02)*log(m02*(1.0+Epsilon)/m02));

        /* A0 */
        output(PV.A0(mu2, 0.0), PV.A0(mu2, mu2*Epsilon));

        /* B0 */
        //output(PV.B0(mu2, 0.0, m02, m12), PV.B0(mu2, p2*Epsilon, m02, m12));
        //output(PV.B0(mu2, p2, 0.0, m12), PV.B0(mu2, p2, p2*Epsilon, m12));
        //output(PV.B0(mu2, p2, m02, 0.0), PV.B0(mu2, p2, m02, p2*Epsilon));
        //output(PV.B0(mu2, p2, m02, m02), PV.B0(mu2, p2, m02, m02*(1.0+Epsilon)));

    }

    cout << PV.B0(m02, m02 + m12 - 0.1, m02, m12) << endl;
    cout << PV.B0(m02, m02 + m12, m02, m12) << endl;
    cout << PV.B0(m02, m02 + m12 + 0.1, m02, m12) << endl;

    return (EXIT_SUCCESS);
}

