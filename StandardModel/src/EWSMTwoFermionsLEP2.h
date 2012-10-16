/* 
 * File:   EWSMTwoFermionsLEP2.h
 * Author: mishima
 */

#ifndef EWSMTWOFERMIONSLEP2_H
#define	EWSMTWOFERMIONSLEP2_H

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
     */
    EWSMTwoFermionsLEP2(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////  

    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] bWEAK with/without weak corrections
     * @param[in] bWEAKBOX with/without weak box corrections
     * @param[in] bQED with/without QED corrections
     * @return the total cross section for e^+ e^- -> l lbar in GeV^{-2}
     */
    double sigma_l(const StandardModel::lepton l, const double s, 
                   const double Mw, const double GammaZ, 
                   const bool bWEAK, const bool bWEAKBOX, const bool bQED) const;
    
    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] bWEAK with/without weak corrections
     * @param[in] bWEAKBOX with/without weak box corrections
     * @param[in] bQED with/without QED corrections
     * @return the total cross section for e^+ e^- -> q qbar in GeV^{-2}
     */
    double sigma_q(const StandardModel::quark q, const double s, 
                   const double Mw, const double GammaZ, 
                   const bool bWEAK, const bool bWEAKBOX, const bool bQED) const;

    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] bWEAK with/without weak corrections
     * @param[in] bWEAKBOX with/without weak box corrections
     * @param[in] bQED with/without QED corrections
     * @return the forward-backward asymmetry for e^+ e^- -> l lbar
     */
    double AFB_l(const StandardModel::lepton l, const double s, 
                 const double Mw, const double GammaZ, 
                 const bool bWEAK, const bool bWEAKBOX, const bool bQED) const;
    
    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] bWEAK with/without weak corrections
     * @param[in] bWEAKBOX with/without weak box corrections
     * @param[in] bQED with/without QED corrections
     * @return the forward-backward asymmetry for e^+ e^- -> q qbar
     */
    double AFB_q(const StandardModel::quark q, const double s, 
                 const double Mw, const double GammaZ, 
                 const bool bWEAK, const bool bWEAKBOX, const bool bQED) const;

    
    ////////////////////////////////////////////////////////////////////////  
private:
    const StandardModel& SM;
    const EWSMOneLoopEW myOneLoopEW;
    const Polylogarithms Polylog;
    const PVfunctions PV;        
    
};

#endif	/* EWSMTWOFERMIONSLEP2_H */

