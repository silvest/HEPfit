/* 
 * File:   EWSMTwoFermionsLEP2.h
 * Author: mishima
 */

#ifndef EWSMTWOFERMIONSLEP2_H
#define	EWSMTWOFERMIONSLEP2_H

#include <gslpp.h>
#include "EWSMOneLoopEW_HV.h"
using namespace gslpp;


/**
 * @class EWSMTwoFermionsLEP2
 * @brief Cross sections and forward-backward asymmetries for e^+e^- -> f fbar at LEP-II
 * 
 * The formulae used here are referred to Hollik's pape, Fortschr. Phys 38 (1990) 3, 165. 
 */
class EWSMTwoFermionsLEP2 {
public:
    
    /**
     * @brief EWSMTwoFermionsLEP2 constructor
     * @param[in] SM_i reference to a StandardModel object
     */
    EWSMTwoFermionsLEP2(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////  
    
    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bDP with/without dressed gauge-boson propagators
     * @param[in] bQED with/without QED corrections
     * @return the total cross section for e^+ e^- -> l lbar in GeV^{-2}
     */
    double sigma_l(const StandardModel::lepton l, const double s,
                   const double Mw, const double GammaZ,
                   const bool bDP=true, const bool bQED=true) const;

    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] bDP with/without dressed gauge-boson propagators
     * @param[in] bQED with/without QED corrections
     * @return the total cross section for e^+ e^- -> q qbar in GeV^{-2}
     */
    double sigma_q(const StandardModel::quark q, const double s,
                   const double Mw, const double GammaZ,
                   const bool bDP=true, const bool bQED=true) const;

    
    ////////////////////////////////////////////////////////////////////////  
private:
    const StandardModel& SM;
    const EWSMOneLoopEW_HV myOneLoopEW_HV;

    
    ////////////////////////////////////////////////////////////////////////  
    
    /**
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] mf mass of f
     * @param[in] Qf electric charge of f
     * @param[in] I3f isospin of f
     * @param[in] Ncf Ncf=3 for f=q; 1 for f=l
     * @param[in] bQED with/without QED corrections
     * @param[in] bQED with/without QED corrections
     * @return the total cross section for e^+ e^- -> f fbar in GeV^{-2}
     */
    double sigma_f(const double s, const double Mw, const double GammaZ,
                   const double mf, const double Qf, const double I3f, const double Ncf, 
                   const bool bDP=true, const bool bQED=true) const;
    
    // Renormalized self-energies
    complex Sigma_hat_ZZ(const double mu, const double s, const double Mw) const;
    complex Sigma_hat_gZ(const double mu, const double s, const double Mw) const;
    complex Sigma_hat_gg(const double mu, const double s, const double Mw) const;
    
    // Dressed gauge-boson propagators
    complex chi_Z(const double mu, const double s, const double Mw) const;
    complex chi_gamma(const double mu, const double s, const double Mw) const;
    complex chi_gammaZ(const double mu, const double s, const double Mw) const;
    
    // QED corrections to e^+e^- -> f fbar cross sections
    double delta() const;
    double Bf(const double s, const double mf) const;
    double gamma_delta(const double s, const double mf, const double Qf) const;
    complex gamma_delta_int(const double s, const double GammaZ, const double mf, const double Qf) const;
    double gamma_delta_res(const double s, const double GammaZ, const double mf, const double Qf) const;
    double gamma_tail(const double s, const double GammaZ) const;
    double gamma_fin(const double s, const double mf, const double Qf) const;    
    double C11V(const double s, const double mf, const double Qf) const;    
    double C11A(const double s, const double mf, const double Qf) const;    
    complex C12V(const double s, const double GammaZ, const double mf, const double Qf) const;    
    complex C12A(const double s, const double mf, const double Qf) const;    
    double C22V(const double s, const double GammaZ, const double mf, const double Qf) const;    
    double C22A(const double s, const double mf, const double Qf) const;     

    
};

#endif	/* EWSMTWOFERMIONSLEP2_H */

