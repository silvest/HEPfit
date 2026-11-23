/* 
 * Copyright (C) 2012 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BSMUMU_H
#define	BSMUMU_H

class StandardModel;
class EvolBsmm;
#include <memory>
#include "ThObservable.h"
#include "QCD.h"
#include "OrderScheme.h"
#include "DGamma.h"

class Bsmumu : public ThObservable {
public:   
    /**
     * constructor
     * @param Flavour
     */
    Bsmumu(const StandardModel& SM_i, int obsFlag, QCD::lepton lep_i=QCD::MU);
    
    /**
     * 
     * @brief hep-ph/9512380v2
     * @return theoretical value of |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    double computeThValue();
    /**
     * @brief The mass-eigenstate rate asymmetry @f$A_{\Delta\Gamma}@f$ of @f$B_s \to \mu^+\mu^-@f$ (1 in the SM).
     */
    double computeAmumu(orders order);
    /**
     * @brief The CP asymmetry @f$S_{\ell\ell}@f$ of @f$B_s \to \mu^+\mu^-@f$ (observable Smumu_Bs_old).
     * @details LHCb convention (De Bruyn et al., arXiv:1204.1737, eqs. (16), (19), with the NP mixing phase as in
     * Buras et al., arXiv:1303.3820):
     * @f$[\Gamma(B_q(t)\to\ell\ell) - \Gamma(\bar B_q(t)\to\ell\ell)]/[\Gamma + \bar\Gamma]
     * = S_{\ell\ell}\sin(\Delta m_q t)/[\cosh(y_q\Gamma_q t) + A_{\Delta\Gamma}\sinh(y_q\Gamma_q t)]@f$ and
     * @f$S_{\ell\ell} = [|P|^2\sin(2\varphi_P - \phi_q^{NP}) - |S|^2\sin(2\varphi_S - \phi_q^{NP})]/(|P|^2+|S|^2)@f$,
     * where @f$\lambda = (q/p)\bar A/A = -1@f$ in the SM and @f$\phi_q^{NP}@f$ is the NP phase of @f$B_q@f$ mixing with
     * the sign of the HEPfit observables Phis_JPsiPhi (@f$\phi_s = \phi_s^{SM} + \phi_s^{NP}@f$, @f$\phi_s^{SM}
     * \simeq -2\beta_s@f$) and SJPsiK (@f$S_{\psi K} = \sin(2\beta + \phi_d^{NP})@f$). @f$S_{\ell\ell} = 0@f$ in the SM;
     * with an SM decay amplitude @f$S_{\ell\ell} = -\sin\phi_q^{NP}@f$, i.e. @f$+\sin(2\phi_{B_s})@f$ and
     * @f$-\sin(2\phi_{B_d})@f$ for the mixing phases getPhiBs(), getPhiBd() (e.g. the NPDF2 parameters PhiBs, PhiBd). (After upstream #24 and until this was
     * fixed the opposite sign, @f$2\,{\rm Im}\lambda/(1+|\lambda|^2)@f$, was returned.)
     */
    double computeSmumu(orders order);
    
    
protected:
    
    /**
     * 
     * @param order
     * @param order_qed
     * @return the short distance contribution to the 
     * |\f$ BR(B_s \rightarrow \mu \bar{\mu}) \f$|
     */
    void computeAmpSq(orders order, orders_qed order_qed, double mu);
    void computeObs(orders order, orders_qed order_qed, bool withMixing = true);
    
private:
    
    QCD::lepton lep; /**< Final leptons type. */
    double beta;
    double mBs;
    double mW;
    double mlep;
    double mb;
    double ms;
    double chiral;
    double absP;
    double argP;
    double absS;
    double argS;
    double ampSq;
    double Amumu;
    double Smumu;
//    double convertFromSingletoDoubleGF;
    double coupling;
    int obs;
    int leptonindex;
    gslpp::complex C_10;
    gslpp::complex C_10p;
    gslpp::complex C_S;
    gslpp::complex C_Sp;
    gslpp::complex C_P;
    gslpp::complex C_Pp;
    
    std::unique_ptr<EvolBsmm> evolbsmm;
    gslpp::vector<gslpp::complex> ** allcoeff;
 

};

#endif	/* BSMUMU_H */
