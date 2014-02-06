/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
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
 * decay widths of the \f$Z\f$ boson, left-right asymmetries and cross sections
 * at the \f$Z\f$ pole. 
 */
class EW : public ThObsType {
public:
    
    /**
     * @brief Constructor.
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
    // Final-state corrections to Z-decay widths

    /*
     * @param[in] q name of a quark.
     * @return non-factorizable EW-QCD corrections in GeV.
     */
    double Delta_EWQCD(const StandardModel::quark q) const;

    /**
     * @param[in] q name of a quark.
     * @return Radiator functions to the vector current due to the
     * final-state QED and QCD corrections.
     */
    double RVq(const StandardModel::quark q) const;

    /**
     * @param[in] q name of a quark.
     * @return Radiator functions to the axial-vector current due to the
     * final-state QED and QCD corrections.
     */
    double RAq(const StandardModel::quark q) const;

    /**
     * @return Singlet vector corrections to the width of Z to hadrons.
     */
    double RVh() const;


    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief The effective leptonic weak mixing angle, \f$\sin^2{\theta_{Eff}^\ell}\f$.
     * @param[in] l name of a lepton
     * @return the effective weak mixing angle for the lepton "l"
     */
    double sin2thetaEff(const StandardModel::lepton l) const;
    
     /**
     * @brief The effective quark weak mixing angle, \f$\sin^2{\theta_{Eff}^q}\f$.
     * @param[in] q name of a quark
     * @return the effective weak mixing angle for the quark "q"
     */
    double sin2thetaEff(const StandardModel::quark q) const;   
    
    /**
     * @brief The \f$Z\to\ell^+\ell^-\f$ partial decay width, \f$\Gamma_\ell\f$.
     * @param[in] l name of a lepton
     * @return the \f$Z\rightarrow \ell^+\ell^-\f$ partial decay width in GeV 
     */
    double Gamma_l(const StandardModel::lepton l) const;
    
    /**
     * @brief The \f$Z\to q\bar{q}\f$ partial decay width, \f$\Gamma_q\f$.
     * @param[in] q name of a quark
     * @return the \f$Z\rightarrow q\bar{q}\f$ partial decay width in GeV
     */
    double Gamma_q(const StandardModel::quark q) const;
    
    /**
     * @brief The \f$Z\f$-boson invisible partial decay width, \f$\Gamma_{inv}\f$.
     * @return the invisible decay width of the \f$Z\f$ boson in GeV
     */
    double Gamma_inv() const;

    /**
     * @brief The \f$Z\to\mbox{hadrons}\f$ partial decay width, \f$\Gamma_{had}\f$.
     * @return the hadronic decay width of the \f$Z\f$ boson in GeV
     */
    double Gamma_had() const;

    /**
     * @brief The total decay width of the \f$Z\f$ boson, \f$\Gamma_Z\f$.
     * @return the total decay width of the \f$Z\f$ boson in GeV
     */
    double Gamma_Z() const;
    
    /**
     * @brief The \f$Z\f$-pole leptonic cross section, \f$\sigma_{lept}^0\f$.
     * @param[in] l name of a lepton
     * @return the cross section for \f$e^+e^- \rightarrow Z \rightarrow l\bar{l}\f$ 
     * at the \f$Z\f$ pole in GeV\f$^{-2}\f$
     */
    double sigma0_l(const StandardModel::lepton l) const;

    /**
     * @brief The \f$Z\f$-pole hadronic cross section, \f$\sigma_h^0\f$.
     * @return the cross section for \f$e^+e^- \rightarrow Z \rightarrow \mathrm{hadrons}\f$
     * at the \f$Z\f$ pole in GeV\f$^{-2}\f$
     */
    double sigma0_had() const; 
 
    /**
     * @brief The \f$Z\f$-pole leptonic left-right asymmetry, \f$A_l\f$.
     * @param[in] l name of a lepton
     * @return the asymmetry parameter for \f$Z\rightarrow l\bar{l}\f$, \f$A_l\f$
     */
    double A_l(const StandardModel::lepton l) const;

    /**
     * @brief The \f$Z\f$-pole quark left-right asymmetry, \f$A_q\f$.
     * @param[in] q name of a quark
     * @return the asymmetry parameter for \f$Z\rightarrow q\bar{q}\f$, \f$A_q\f$
     */
    double A_q(const StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////
protected:
    const StandardModel& SM;///> A reference to an object of the StandardModel class.
    const EW_NPZff myEW_NPZff;

};

#endif	/* EW_H */

