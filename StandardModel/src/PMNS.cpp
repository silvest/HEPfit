/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "PMNS.h"

PMNS::PMNS() : U(3, 3)
{}

void PMNS::computePMNS(double s12_v, double s13_v, double s23_v, double delta_v, double alpha21_v, double alpha31_v)
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