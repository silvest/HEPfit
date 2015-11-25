/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EW_ABC_H
#define	EW_ABC_H

#include "gslpp.h"
#include "NPEpsilons.h"

/**
 * @class EW_ABC
 * @ingroup NewPhysics
 * @brief A test class for the electroweak precision observables. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A test class for the electroweak precision observables with the epsilon
 * parameters based on Eqs.(7)-(14) of IJMPA, 7, 1031-1058 (1998) by Altarelli,
 * Barbieri and Caravaglios. This class is used for test. If bAlternative=true
 * is passed to observables, the formulae in Eqs.(16)-(20) of IJMP, A7, 1031-1058
 * (1998) by the same authors are adopted instead. 
 */
class EW_ABC {
public:

    EW_ABC(const NPEpsilons& NPE_i);

    double Mw(const bool bAlternative) const;
    double Gamma_l(StandardModel::lepton l) const;
    double Gamma_q(QCD::quark q) const;
    double Gamma_b() const;

    double Gamma_had() const;
    double GammaZ(const bool bAlternative) const;

    double R_l(const bool bAlternative) const;
    double R_c() const;
    double R_b(const bool bAlternative) const;

    double sigma0_had(const bool bAlternative) const;

    double A_l(StandardModel::lepton l, const bool bAlternative) const;
    double A_q(QCD::quark q) const;
    double A_b() const;

    double AFB_l(StandardModel::lepton l, const bool bAlternative) const;
    double AFB_c() const;
    double AFB_b() const;

    double sin2thetaEff(const bool bAlternative) const;


    ////////////////////////////////////////////////////////////////////////

    double eps1() const;
    double eps2() const;
    double eps3() const;
    double epsb() const;


    ////////////////////////////////////////////////////////////////////////

    gslpp::complex gVl(StandardModel::lepton l) const;
    gslpp::complex gAl(StandardModel::lepton l) const;
    gslpp::complex gVl_over_gAl(StandardModel::lepton l) const;
    gslpp::complex gVq(QCD::quark q) const;
    gslpp::complex gAq(QCD::quark q) const;
    gslpp::complex gVq_over_gAq(QCD::quark q) const;
    gslpp::complex gVb() const;
    gslpp::complex gAb() const;
    gslpp::complex gVb_over_gAb() const;


    ////////////////////////////////////////////////////////////////////////
private:
    const NPEpsilons& NPE;

};

#endif	/* EW_ABC_H */

