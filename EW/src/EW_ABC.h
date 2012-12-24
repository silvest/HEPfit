/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EW_ABC_H
#define	EW_ABC_H

#include <gslpp.h>
#include <StandardModel.h>
#include <EWepsilons.h>
using namespace gslpp;


/**
 * @class EW_ABC
 * @brief a class for the electroweak precision observables with the epsilon 
 * parameters based on IJMPA, 7, 1031-1058 (1998) by G.Altarelli, R.Barbieri 
 * and F.Caravaglios. 
 */  
class EW_ABC {
public:

    EW_ABC(const StandardModel& SM_i) : SM(SM_i), myEWepsilons(SM_i) {};

    double Mw(const double eps1, const double eps2, const double eps3) const;
    double Gamma_l(StandardModel::lepton l, const double eps1, const double eps3) const;
    double Gamma_q(StandardModel::quark q, const double eps1, const double eps3) const;
    double Gamma_b(const double eps1, const double eps3, const double epsb) const;
    
    double Gamma_had(const double eps1, const double eps3, const double epsb) const;
    double GammaZ(const double eps1, const double eps3, const double epsb) const;    
    
    double R_l(const double eps1, const double eps3, const double epsb) const;
    double R_c(const double eps1, const double eps3, const double epsb) const;    
    double R_b(const double eps1, const double eps3, const double epsb) const;    

    double sigma0_had(const double eps1, const double eps3, const double epsb) const;
    
    double A_l(StandardModel::lepton l, const double eps1, const double eps3) const;
    double A_q(StandardModel::quark q, const double eps1, const double eps3) const;
    double A_b(const double eps1, const double eps3, const double epsb) const;
    
    double AFB_l(StandardModel::lepton l, const double eps1, const double eps3) const;
    double AFB_c(const double eps1, const double eps3) const;
    double AFB_b(const double eps1, const double eps3, const double epsb) const;
    
    double sin2thetaEff(const double eps1, const double eps3) const;

    //////////////////////////////////////////////////////////////////////// 
    
    complex gVl(StandardModel::lepton l, const double eps1, const double eps3) const;
    complex gAl(StandardModel::lepton l, const double eps1) const;
    complex gVl_over_gAl(StandardModel::lepton l, const double eps1, const double eps3) const;
    complex gVq(StandardModel::quark q, const double eps1, const double eps3) const;
    complex gAq(StandardModel::quark q, const double eps1) const;
    complex gVq_over_gAq(StandardModel::quark q, const double eps1, const double eps3) const;
    complex gVb(const double eps1, const double eps3, const double epsb) const;
    complex gAb(const double eps1, const double epsb) const;
    complex gVb_over_gAb(const double eps1, const double eps3, const double epsb) const;
    
private:
    const StandardModel& SM;
    const EWepsilons myEWepsilons;

};

#endif	/* EW_ABC_H */

