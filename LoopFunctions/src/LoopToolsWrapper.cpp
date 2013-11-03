/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <complex>
#include <clooptools.h>
#include "LoopToolsWrapper.h"

static bool LoopToolsInit = false;

LoopToolsWrapper::LoopToolsWrapper()
{
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
}

LoopToolsWrapper::~LoopToolsWrapper()
{
    // for debug
    //std::cout << std::endl
    //          << "************* LoopTools ****************" << std::endl;
    //ltexi();
    //std::cout << "****************************************" << std::endl;
}

double LoopToolsWrapper::PV_A0(const double mu2, const double m2) const
{
    setmudim(mu2);
    std::complex<double> A0val = A0(m2);
    return ( - A0val.real() );
}
    
complex LoopToolsWrapper::PV_B0(const double mu2, const double p2,
                                const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B0val = B0(p2, m02, m12);
    return complex( B0val.real(), B0val.imag(), false );
}

complex LoopToolsWrapper::PV_B1(const double mu2, const double p2,
                                const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B1val = B1(p2, m02, m12);
    return complex( B1val.real(), B1val.imag(), false );
}

complex LoopToolsWrapper::PV_B21(const double mu2, const double p2,
                                 const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B21val = B11(p2, m02, m12);
    return complex( B21val.real(), B21val.imag(), false );
}

complex LoopToolsWrapper::PV_B22(const double mu2, const double p2,
                                 const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B22val = B00(p2, m02, m12);
    return complex( B22val.real(), B22val.imag(), false );
}

complex LoopToolsWrapper::PV_B0p(const double muIR2, const double p2,
                                 const double m02, const double m12) const
{
    setmudim(muIR2);
    std::complex<double> B0pval = DB0(p2, m02, m12);
    return complex( B0pval.real(), B0pval.imag(), false );
}

complex LoopToolsWrapper::PV_B1p(const double mu2, const double p2,
                                 const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B1pval = DB1(p2, m02, m12);
    return complex( B1pval.real(), B1pval.imag(), false );
}

complex LoopToolsWrapper::PV_B21p(const double mu2, const double p2,
                                  const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B21pval = DB11(p2, m02, m12);
    return complex( B21pval.real(), B21pval.imag(), false );
}

complex LoopToolsWrapper::PV_B22p(const double mu2, const double p2,
                                  const double m02, const double m12) const
{
    setmudim(mu2);
    std::complex<double> B22pval = DB00(p2, m02, m12);
    return complex( B22pval.real(), B22pval.imag(), false );
}

complex LoopToolsWrapper::PV_C0(const double p2,
                                const double m02, const double m12, const double m22) const
{
    std::complex<double> C0val = C0(0.0, 0.0, p2, m02, m12, m22);
    return complex( - C0val.real(), - C0val.imag(), false );
}
    
complex LoopToolsWrapper::PV_D0(const double s, const double t, const double m02,
                                const double m12, const double m22, const double m32) const
{
    std::complex<double> D0val = D0(0.0, 0.0, 0.0, 0.0, s, t, m02, m12, m22, m32);
    return complex( D0val.real(), D0val.imag(), false );
}


