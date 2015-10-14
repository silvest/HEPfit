/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "PMNS.h"

PMNS::PMNS(std::string H_i)
{
    H = H_i;
}

PMNS::PMNS(const PMNS& orig)
{
    s12 = orig.s12;
    s13 = orig.s13;
    s23 = orig.s23;
    delta = orig.delta;
    c12 = orig.c12;
    c23 = orig.c23;
    c13 = orig.c13;

    U11 = orig.U11;
    U12 = orig.U12;
    U13 = orig.U13;
    U21 = orig.U21;
    U22 = orig.U22;
    U23 = orig.U23;
    U31 = orig.U31;
    U32 = orig.U32;
    U33 = orig.U33;

}

PMNS::~PMNS()
{
}

void PMNS::getPMNS(gslpp::matrix<gslpp::complex> & x) const
{
    x.assign(0, 0, U11);
    x.assign(0, 1, U12);
    x.assign(0, 2, U13);
    x.assign(1, 0, U21);
    x.assign(1, 1, U22);
    x.assign(1, 2, U23);
    x.assign(2, 0, U31);
    x.assign(2, 1, U32);
    x.assign(2, 2, U33);
}

void PMNS::setPMNS(double s12, double s13, double s23, double delta)
{
    c12 = sqrt(1. - pow(s12, 2.));
    c23 = sqrt(1. - pow(s23, 2.));
    c13 = sqrt(1. - pow(s13, 2.));

    U11 = gslpp::complex(c12*c13, 0.);
    U12 = gslpp::complex(s12*c13, 0.);
    U13 = gslpp::complex(s13, -delta, true);

    U21 = -s12 * c23 - gslpp::complex(c12 * s23*s13, delta, true);
    U22 = c12 * c23 - gslpp::complex(s12 * s23*s13, delta, true);
    U23 = gslpp::complex(s23*c13, 0.);

    U31 = s12 * s23 - gslpp::complex(c12 * c23*s13, delta, true);
    U32 = -c12 * s23 - gslpp::complex(s12 * c23*s13, delta, true);
    U33 = gslpp::complex(c23*c13, 0.);

    return;
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

std::string PMNS::getHierarchy()
{
    return H;
}

// J_CP

double PMNS::getJcp()
{
    throw std::runtime_error("PMNS()::Jcp Not implemented!");
    return 1;
}

// Sides
/*double CKM::getRt()
{
    return (Vtd*Vtb.conjugate()/(Vcd*Vcb.conjugate())).abs();
}

double CKM::getRb() 
{
    return (Vud*Vub.conjugate()/(Vcd*Vcb.conjugate())).abs();
}

double CKM::getRts() 
{
    return (Vts*Vtb.conjugate()/(Vcs*Vcb.conjugate())).abs();
}*/
