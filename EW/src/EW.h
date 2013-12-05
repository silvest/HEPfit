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
#include "EW_CHMN.h"
#include "EW_ABC.h"
#include "EW_BURGESS.h"

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
 * 
 * 
 */
class EW : public ThObsType {
public:
    
    /**
     * EWDEFAULT : use our formulae
     * EWCHMN : use the formulae by Cho, Hagiwara, Matsumoto and Nomura in EW_CHMN class.
     * EWBURGESS : use the formulae by Burgess et al. 
     * EWABC : use the formulae in Eqs.(7)-(14) of IJMP, A7, 1031-1058 (1998) by Altarelli, Barbieri and Caravaglios. 
     * EWABC2 : use the approximate formulae in Eqs.(16)-(20) of IJMP, A7, 1031-1058 (1998) by Altarelli, Barbieri and Caravaglios. 
     */
    enum EWTYPE {EWDEFAULT, EWCHMN, EWBURGESS, EWABC, EWABC2};
    
    /**
     * @brief A constructor.
     * @param[in] SM_i A reference to an object of StandardModel class
     */
    EW(const StandardModel& SM_i);

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @attention This function has to be called after initializing model flags 
     * by InputParser::ReadParameters(). 
     * @return 
     */
    EWTYPE getEWTYPE() const;

    bool checkLEP1NP() const;
    
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

    const EW_ABC getMyEW_ABC() const 
    {
        return myEW_ABC;
    }

    const EW_CHMN getMyEW_CHMN() const 
    {
        return myEW_CHMN;
    }

    const EW_BURGESS getMyEW_BURGESS() const
    {
        return myEW_BURGESS;
    }

    /**
     * @param[in] l lepton type
     * @return the electric charge of the lepton "l"
     */
    double Ql(const StandardModel::lepton l) const;    

    /**
     * @param[in] q quark type
     * @return the electric charge of the quark "q"
     */
    double Qq(const StandardModel::quark q) const;    

    /**
     * @return \f$\alpha(M_Z^2)\f$ 
     */
    double alpha() const;

    /**
     * @return the \f$W\f$-boson mass in the SM 
     */
    double Mw_SM() const;    

    /**
     * @return the value of \f$\sin^2{\theta_W}\f$ in the SM
     */
    double sW2_SM() const;

    /**
     * @return the value of \f$\cos^2{\theta_W}\f$ in the SM
     */
    double cW2_SM() const;


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
    const EW_CHMN myEW_CHMN;
    const EW_ABC myEW_ABC;
    const EW_BURGESS myEW_BURGESS;
    
    ////////////////////////////////////////////////////////////////////////   
    
private:

};

#endif	/* EW_H */

