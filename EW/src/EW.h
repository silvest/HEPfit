/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EW_H
#define	EW_H

#include <stdexcept>
#include <ThObsType.h>
#include <StandardModel.h>
#include "EW_NPZff.h"

using namespace gslpp;

/**
 * @class EW
 * @ingroup EW 
 * @brief The base class for the electroweak precision observables. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class includes basic functions for the computation of the SM
 * predictions for electroweak precision pseudo observables, such as partial
 * decay widths of the \f$Z\f$ boson. 
 */
class EW : public ThObsType {
public:
    
    /**
     * @brief A constructor.
     * @param[in] SM_i a reference to an object of StandardModel class
     */
    EW(const StandardModel& SM_i);

    
    //////////////////////////////////////////////////////////////////////// 

    bool checkNPZff() const;
    
    /**
     * @return a reference to the StandardModel object in the current class
     */
    const StandardModel& getSM() const 
    {
        return SM;
    } 

    const EW_NPZff getMyEW_NPZff() const
    {
        return myEW_NPZff;
    }


    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @param[in] l name of a lepton
     * @return the effective weak mixing angle for the lepton "l"
     */
    double sin2thetaEff(const StandardModel::lepton l) const;
    
     /**
     * @param[in] q name of a quark
     * @return the effective weak mixing angle for the quark "q"
     */
    double sin2thetaEff(const StandardModel::quark q) const;   
    
    /**
     * @param[in] l name of a lepton
     * @return the \f$Z\rightarrow \mbox{l}\bar{\mbox{l}}\f$ partial decay width 
     */
    double Gamma_l(const StandardModel::lepton l) const;
    
    /**
     * @param[in] q name of a quark
     * @return the \f$Z\rightarrow \mbox{q}\bar{\mbox{q}}\f$ partial decay width 
     */
    double Gamma_q(const StandardModel::quark q) const;
    
    /**
     * @return the invisible decay width of the \f$Z\f$ boson
     */
    double Gamma_inv() const;

    /**
     * @return the hadronic decay width of the \f$Z\f$ boson
     */
    double Gamma_had() const;

    /**
     * @return the total decay width of the \f$Z\f$ boson
     */
    double Gamma_Z() const;
    
    /**
     * @param[in] l name of a lepton
     * @return the cross section for \f$e^+e^- \rightarrow Z \rightarrow l\bar{l}\f$
     */
    double sigma0_l(const StandardModel::lepton l) const;

    /**
     * @return the cross section for \f$e^+e^- \rightarrow Z \rightarrow \mathrm{hadrons}f$
     */
    double sigma0_had() const; 
 
    /**
     * @param[in] l name of a lepton
     * @return the asymmetry parameter for \f$Z\rightarrow l\bar{l}\f$, \f$A_l\f$
     */
    double A_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the asymmetry parameter for \f$Z\rightarrow q\bar{q}\f$, \f$A_q\f$
     */
    double A_q(const StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////
protected:
    const StandardModel& SM;
    const EW_NPZff myEW_NPZff;

};

#endif	/* EW_H */

