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
#include "EW_CHMN.h"
#include "EW_ABC.h"

using namespace gslpp;

/**
 * @class EW
 * @ingroup EW 
 * @brief The base class for the electroweak precision observables at the Z pole. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
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
    
    /**
     * @return boolean: true for the case where the oblique parameters are employed. 
     */
    bool checkModelForSTU() const;
    
    /**
     * @return a reference to the StandardModel object in the current class
     */
    const StandardModel& getSM() const 
    {
        return SM;
    } 

    const EW_ABC getMyEW_ABC() const 
    {
        return myEW_ABC;
    }

    const EW_CHMN getMyEW_CHMN() const 
    {
        return myEW_CHMN;
    }
    
    /**
     * @param[in] l lepton
     * @return electric charge of a lepton "l"
     */
    double Ql(const StandardModel::lepton l) const;    

    /**
     * @param[in] q quark
     * @return electric charge of a quark "q"
     */
    double Qq(const StandardModel::quark q) const;    

    /**
     * @return alpha(Mz2) 
     */
    double alpha() const;

    /**
     * @return the W boson mass in the SM. 
     */
    double Mw_SM() const;    

    /**
     * @return sin^2 theta_W in the SM.
     */
    double sW2_SM() const;

    /**
     * @return cos^2 theta_W in the SM.
     */
    double cW2_SM() const;

    /**
     * @return the oblique parameters S
     */
    double S() const 
    {
        return ( SM.obliqueS() );
    }
    
    /**
     * @return the oblique parameters T
     */    
    double T() const 
    {
        return ( SM.obliqueT() );
    }
    
    /**
     * @return the oblique parameters U
     */    
    double U() const 
    {
        return ( SM.obliqueU() );
    }

    /**
     * @return Oblique parameter \hat{S}
     */
    double Shat() const 
    {
        return ( SM.obliqueShat() );
    }

    /**
     * @return Oblique parameter \hat{T}
     */
    double That() const 
    {
        return ( SM.obliqueThat() );
    }

    /**
     * @return Oblique parameter \hat{U}
     */
    double Uhat() const 
    {
        return ( SM.obliqueUhat() );
    }

    /**
     * @return the oblique parameters V
     */    
    double V() const 
    {
        return ( SM.obliqueV() );
    }
    
    /**
     * @return the oblique parameters W
     */    
    double W() const 
    {
        return ( SM.obliqueW() );
    }
    
    /**
     * @return the oblique parameters X
     */    
    double X() const 
    {
        return ( SM.obliqueX() );
    }
    
    /**
     * @return the oblique parameters Y
     */    
    double Y() const 
    {
        return ( SM.obliqueY() );
    }    

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @param[in] l name of a lepton
     * @return the effective mixing angle for lepton "l"
     */
    double sin2thetaEff(const StandardModel::lepton l) const;
    
     /**
     * @param[in] q name of a quark
     * @return the effective mixing angle for quark "q"
     */
    double sin2thetaEff(const StandardModel::quark q) const;   
    
    /**
     * @param[in] l name of a lepton
     * @return the partial width of Z decay into an l\bar{l} pair 
     */
    double Gamma_l(const StandardModel::lepton l) const;
    
    /**
     * @param[in] q name of a quark
     * @return the partial width of Z decay into a q\bar{q} pair 
     */
    double Gamma_q(const StandardModel::quark q) const;
    
    /**
     * @return the partial width of Z decay into neutrinos
     */
    double Gamma_inv() const;

    /**
     * @return the hadronic width of the Z boson
     */
    double Gamma_had() const;

    /**
     * @return the total width of the Z boson
     */
    double Gamma_Z() const;
    
    /**
     * @param[in] l name of a lepton
     * @return the cross section for e^+e^- -> Z -> l\bar{l}
     */
    double sigma0_l(const StandardModel::lepton l) const;

    /**
     * @return the cross section e^+e^- -> Z -> hadrons
     */
    double sigma0_had() const; 
 
    /**
     * @param[in] l name of a lepton
     * @return asymmetry parameter for Z->l\bar{l}
     */
    double A_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return asymmetry parameter for Z->q\bar{q}
     */
    double A_q(const StandardModel::quark q) const;

    
    ////////////////////////////////////////////////////////////////////////   
    
protected:
    bool bDebug; // true for debugging    
    
    const StandardModel& SM;
    const EW_CHMN myEW_CHMN;
    const EW_ABC myEW_ABC;
    
    ////////////////////////////////////////////////////////////////////////   
    
private:

};

#endif	/* EW_H */

