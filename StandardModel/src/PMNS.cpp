/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "PMNS.h"

PMNS::PMNS() : U(3, 3)
{
}

void PMNS::setPMNS(double s12_v, double s13_v, double s23_v, double delta_v, double alpha21_v, double alpha31_v)
{
    s12 = s12_v;
    s13 = s13_v;
    s23 = s23_v;
    c12 = sqrt(1. - s12 * s12);
    c13 = sqrt(1. - s13 * s13);
    c23 = sqrt(1. - s23 * s23);
    delta = delta_v;
    alpha21 = alpha21_v;
    alpha31 = alpha31_v;
    
    U.assign(0, 0, c12*c13);
    U.assign(0, 1, gslpp::complex(s12*c13, alpha21 / 2., true));
    U.assign(0, 2, gslpp::complex(s13, -delta + alpha31 / 2., true));

    U.assign(1, 0 , -s12 * c23 - gslpp::complex(c12 * s23*s13, delta, true));
    U.assign(1, 1, (c12 * c23 - gslpp::complex(s12 * s23*s13, delta, true)) *
                   gslpp::complex(1., alpha21 / 2., true));
    U.assign(1, 2, gslpp::complex(s23*c13, alpha31 / 2., true));

    U.assign(2, 0, s12 * s23 - gslpp::complex(c12 * c23*s13, delta, true));
    U.assign(2, 1, (-c12 * s23 - gslpp::complex(s12 * c23*s13, delta, true)) *
                    gslpp::complex(1., alpha21 / 2., true));
    U.assign(2, 2, gslpp::complex(c23*c13, alpha31 / 2., true));

}

// Gilman parameterization

double PMNS::gets12()
{
    return s12;
}

double PMNS::gets13()
{
    return s13;
}

double PMNS::gets23()
{
    return s23;
}

double PMNS::getc12()
{
    return c12;
}

double PMNS::getc23()
{
    return c23;
}

double PMNS::getc13()
{
    return c13;
}

double PMNS::getdelta()
{
    return delta;
}

double PMNS::getalpha21()
{
    return alpha21;
}

double PMNS::getalpha31()
{
    return alpha31;
}


//Absolute values of PMNS elements

double PMNS::getUe1()
{
    return U(0, 0).abs();
}

double PMNS::getUe2()
{
    return U(0, 1).abs();
}

double PMNS::getUe3()
{
    return U(0, 2).abs();
}

double PMNS::getUmu1()
{
    return U(1, 0).abs();
}

double PMNS::getUmu2()
{
    return U(1, 1).abs();
}

double PMNS::getUmu3()
{
    return U(1, 2).abs();
}

double PMNS::getUtau1()
{
    return U(2, 0).abs();
}

double PMNS::getUtau2()
{
    return U(2, 1).abs();
}

double PMNS::getUtau3()
{
    return U(2, 2).abs();
}

// Phases

double PMNS::getArgUe1()
{
    return U(0, 0).arg();
}

double PMNS::getArgUe2()
{
    return U(0, 1).arg();
}

double PMNS::getArgUe3()
{
    return U(0, 2).arg();
}

double PMNS::getArgUmu1()
{
    return U(1, 0).arg();
}

double PMNS::getArgUmu2()
{
    return U(1, 1).arg();
}

double PMNS::getArgUmu3()
{
    return U(1, 2).arg();
}

double PMNS::getArgUtau1()
{
    return U(2, 0).arg();
}

double PMNS::getArgUtau2()
{
    return U(2, 1).arg();
}

double PMNS::getArgUtau3()
{
    return U(2, 2).arg();
}


//Complex values of PMNS elements

gslpp::complex PMNS::U_e1()
{
    return U(0, 0);
}

gslpp::complex PMNS::U_e2()
{
    return U(0, 1);
}

gslpp::complex PMNS::U_e3()
{
    return U(0, 2);
}

gslpp::complex PMNS::U_mu1()
{
    return U(1, 0);
}

gslpp::complex PMNS::U_mu2()
{
    return U(1, 1);
}

gslpp::complex PMNS::U_mu3()
{
    return U(1, 2);
}

gslpp::complex PMNS::U_tau1()
{
    return U(2, 0);
}

gslpp::complex PMNS::U_tau2()
{
    return U(2, 1);
}

gslpp::complex PMNS::U_tau3()
{
    return U(2, 2);
}
