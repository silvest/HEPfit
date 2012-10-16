/* 
 * File:   EWSMTwoFermionsLEP2.h
 * Author: mishima
 */

#ifndef EWSMTWOFERMIONSLEP2_H
#define	EWSMTWOFERMIONSLEP2_H

#include <map>
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
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] flags a set of flags to control the inclusions of higher-order corrections
     * @return the total cross section for e^+ e^- -> l lbar in GeV^{-2}
     */
    double sigma_l(const StandardModel::lepton l, const double s, 
                   const double Mw, const double GammaZ, 
                   const std::map<std::string, bool>& flags) const;
    
    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] flags a set of flags to control the inclusions of higher-order corrections
     * @return the total cross section for e^+ e^- -> q qbar in GeV^{-2}
     */
    double sigma_q(const StandardModel::quark q, const double s, 
                   const double Mw, const double GammaZ, 
                   const std::map<std::string, bool>& flags) const;

    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] flags a set of flags to control the inclusions of higher-order corrections
     * @return the forward-backward asymmetry for e^+ e^- -> l lbar
     */
    double AFB_l(const StandardModel::lepton l, const double s, 
                 const double Mw, const double GammaZ, 
                 const std::map<std::string, bool>& flags) const;
    
    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width (used in the Born approximation/in the QED corrections)
     * @param[in] flags a set of flags to control the inclusions of higher-order corrections
     * @return the forward-backward asymmetry for e^+ e^- -> q qbar
     */
    double AFB_q(const StandardModel::quark q, const double s, 
                 const double Mw, const double GammaZ, 
                 const std::map<std::string, bool>& flags) const;


    ////////////////////////////////////////////////////////////////////////  

    complex chi_Z(const double s, const double Mw, const double GammaZ) const;


    
    
    ////////////////////////////////////////////////////////////////////////  
private:
    bool bDebug; // true for debugging    
    
    const StandardModel& SM;
    const EWSMOneLoopEW myOneLoopEW;
    const Polylogarithms Polylog;
    const PVfunctions PV;        
    
};

#endif	/* EWSMTWOFERMIONSLEP2_H */

