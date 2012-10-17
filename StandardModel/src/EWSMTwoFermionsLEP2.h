/* 
 * File:   EWSMTwoFermionsLEP2.h
 * Author: mishima
 */

#ifndef EWSMTWOFERMIONSLEP2_H
#define	EWSMTWOFERMIONSLEP2_H

#include <string>
#include <gslpp.h>
#include <Polylogarithms.h>
#include <PVfunctions.h>
#include "EWSMOneLoopEW.h"
using namespace gslpp;


/**
 * @class EWSMTwoFermionsLEP2
 * @brief Cross sections and forward-backward asymmetries for e^+e^- -> f fbar at LEP-II
 * 
 * Formulae used in the current class are calculated in the unitary gauge. 
 */
class EWSMTwoFermionsLEP2 {
public:

    /**
     * @brief EWSMTwoFermionsLEP2 constructor
     * @param[in] SM_i reference to a StandardModel object
     * @param[in] bDebug_i boolean value for debugging (true for debugging)
     */
    EWSMTwoFermionsLEP2(const StandardModel& SM_i, const bool bDebug_i=false);

    
    ////////////////////////////////////////////////////////////////////////  

    /**
     * @brief set a flag to control radiative corrections
     * @param[in] str "Weak", "WeakBox", "ISR", "QEDFSR" or "QCDFSR"
     * @param[in] flag boolean variable
     */
    void setFlag(const std::string str, const bool flag);
    

    ////////////////////////////////////////////////////////////////////////  
    
    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @return the total cross section for e^+ e^- -> l lbar in GeV^{-2}
     */
    double sigma_l(const StandardModel::lepton l, const double s, 
                   const double Mw, const double GammaZ) const;
    
    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @return the total cross section for e^+ e^- -> q qbar in GeV^{-2}
     */
    double sigma_q(const StandardModel::quark q, const double s, 
                   const double Mw, const double GammaZ) const;

    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @return the forward-backward asymmetry for e^+ e^- -> l lbar
     */
    double AFB_l(const StandardModel::lepton l, const double s, 
                 const double Mw, const double GammaZ) const;
    
    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @return the forward-backward asymmetry for e^+ e^- -> q qbar
     */
    double AFB_q(const StandardModel::quark q, const double s, 
                 const double Mw, const double GammaZ) const;

    
    ////////////////////////////////////////////////////////////////////////    

    /**
     * @param[in] l name of lepton
     * @return mass of lepton
     */
    double ml(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of quark
     * @param[in] mu renormalization scale
     * @param[in] order (=LO, NLO, NNLO, FULLNLO[defalut], FULLNNLO)
     * @return the MSbar mass of u, d, s, c, b or the pole mass of t
     */
    double mq(const StandardModel::quark q, const double mu, 
              const orders order=FULLNLO) const;

        
    ////////////////////////////////////////////////////////////////////////  

    complex Vpol_inv(const double s) const;
    
    complex chi_Z(const double s, const double Mw, const double GammaZ) const;

    complex rho_ef(const double s, const double Mw, const double I3f, 
                   const double Qf, const double mfp) const;
    complex kappa_e(const double s, const double Mw, const double I3f, 
                    const double Qf) const;
    complex kappa_f(const double s, const double Mw, const double I3f, 
                    const double Qf, const double mfp) const;
    complex kappa_ef(const double s, const double Mw, const double I3f, 
                     const double Qf, const double mfp) const;

    complex I2e(const double s, const double Mw, const double I3f, 
                const double Qf) const;
    complex I2f(const double s, const double Mw, const double I3f, 
                const double Qf, const double mfp) const;
    
    complex G_e(const double s, const double Mw, const double I3f, 
                const double Qf) const;
    complex G_f(const double s, const double Mw, const double I3f, 
                const double Qf, const double mfp) const;
    complex G_ef(const double s, const double Mw, const double I3f, 
                 const double Qf, const double mfp) const;
    
    double G_1(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mfp) const;
    double G_2(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mfp) const;
    double G_3(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mfp) const; 

    double sigma(const double s, const double Mw, const double GammaZ, 
                 const double I3f, const double Qf, const double mfp,
                 const double mf, const double Ncf) const;
    
    double AFB(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mfp,
               const double mf) const;

    
    ////////////////////////////////////////////////////////////////////////  
private:
    const bool bDebug; // true for debugging    
    bool bWeak, bWeakBox, bISR, bQEDFSR, bQCDFSR; // flags for radiative corrections
    
    const StandardModel& SM;
    const EWSMOneLoopEW myOneLoopEW;
    const Polylogarithms Polylog;
    const PVfunctions PV;        
    
};

#endif	/* EWSMTWOFERMIONSLEP2_H */

