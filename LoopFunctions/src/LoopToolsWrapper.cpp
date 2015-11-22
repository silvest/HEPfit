/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <stdexcept>
#include <complex>
#include "PVfunctions.h" // needed for USE_LOOPTOOLS macro
#include "LoopToolsWrapper.h"
#ifdef USE_LOOPTOOLS
#include <clooptools.h>
static bool LoopToolsInit = false;
#endif


LoopToolsWrapper::LoopToolsWrapper()
{
#ifdef USE_LOOPTOOLS
    if (!LoopToolsInit) {
        //std::cout << std::endl;
        ltini();
        std::cout << std::endl;
        LoopToolsInit = true;
    }

    /* set the photon mass regulating IR divergences to 0, which means that
     * dimensional regularization is employed for IR divergences and that 
     * the finite piece is returned from an IR-divergent loop function. */
    setlambda(0.0);
#endif
}

LoopToolsWrapper::~LoopToolsWrapper()
{
    //#ifdef USE_LOOPTOOLS
    // for debug
    //std::cout << std::endl
    //          << "************* LoopTools ****************" << std::endl;
    //ltexi();
    //std::cout << "****************************************" << std::endl;
    //#endif
}

#ifdef USE_LOOPTOOLS

double LoopToolsWrapper::PV_A0(const double mu2, const double m2) const
{
    setmudim(mu2);
    std::complex<double> A0val = A0(m2);
    return ( A0val.real() );
}
    
gslpp::complex LoopToolsWrapper::PV_B0(const double mu2, const double p2,
                                const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B0val = B0(p2, m02, m12);
    return gslpp::complex( B0val.real(), B0val.imag(), false );
}

gslpp::complex LoopToolsWrapper::PV_B1(const double mu2, const double p2,
                                const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B1val = B1(p2, m02, m12);
    return gslpp::complex( B1val.real(), B1val.imag(), false );
}

gslpp::complex LoopToolsWrapper::PV_B11(const double mu2, const double p2,
                                 const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B11val = B11(p2, m02, m12);
    return gslpp::complex( B11val.real(), B11val.imag(), false );
}

gslpp::complex LoopToolsWrapper::PV_B00(const double mu2, const double p2,
                                 const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B00val = B00(p2, m02, m12);
    return gslpp::complex( B00val.real(), B00val.imag(), false );
}

gslpp::complex LoopToolsWrapper::PV_B0p(const double muIR2, const double p2,
                                 const double m02, const double m12) const
{
    setmudim(muIR2);
    std::complex<double> B0pval = DB0(p2, m02, m12);
    return gslpp::complex( B0pval.real(), B0pval.imag(), false );
}

gslpp::complex LoopToolsWrapper::PV_B1p(const double mu2, const double p2,
                                 const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B1pval = DB1(p2, m02, m12);
    return gslpp::complex( B1pval.real(), B1pval.imag(), false );
}

gslpp::complex LoopToolsWrapper::PV_B11p(const double mu2, const double p2,
                                  const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B11pval = DB11(p2, m02, m12);
    return gslpp::complex( B11pval.real(), B11pval.imag(), false );
}

gslpp::complex LoopToolsWrapper::PV_B00p(const double mu2, const double p2,
                                  const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B00pval = DB00(p2, m02, m12);
    return gslpp::complex( B00pval.real(), B00pval.imag(), false );
}

gslpp::complex LoopToolsWrapper::PV_C0(const double p2,
                                const double m02, const double m12, const double m22) const
{
    std::complex<double> C0val = C0(0.0, 0.0, p2, m02, m12, m22);
    return gslpp::complex( C0val.real(), C0val.imag(), false );
}
    
gslpp::complex LoopToolsWrapper::PV_D0(const double s, const double t, const double m02,
                                const double m12, const double m22, const double m32) const
{
    std::complex<double> D0val = D0(0.0, 0.0, 0.0, 0.0, s, t, m02, m12, m22, m32);
    return gslpp::complex( D0val.real(), D0val.imag(), false );
}

gslpp::complex LoopToolsWrapper::PV_D00(const double s, const double t, const double m02,
                                const double m12, const double m22, const double m32) const
{
// cannot be compiled?
//    std::complex<double> D00val = D0i(dd00, 0.0, 0.0, 0.0, 0.0, s, t, m02, m12, m22, m32);
//    return gslpp::complex( D00val.real(), D00val.imag(), false );

    throw std::runtime_error("LoopToolsWrapper::PV_D00: Not implemented!");
}

#endif // #ifdef USE_LOOPTOOLS


