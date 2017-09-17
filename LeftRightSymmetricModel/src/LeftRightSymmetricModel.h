/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEFTRIGHTSYMMETRICMODEL_H
#define	LEFTRIGHTSYMMETRICMODEL_H

#include "gslpp.h"
#include "StandardModel.h"
#include "LeftRightSymmetricModelMatching.h"

class LRSMquantities; // forward reference to the class that generates the neutral scalar spectrum

/**
 * @class LeftRightSymmetricModel
 * @ingroup LeftRightSymmetricModel
 * @brief A base class for generic Left Right Symmetric Models.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 *
 * 
 * @anchor LeftRightSymmetricModel
 * <h3>%Model parameters</h3>
 */
class LeftRightSymmetricModel: public StandardModel {
public:

    friend class LRSMquantities;

    static const int NLeftRightSymmetricModelvars = 25;
    static const std::string LeftRightSymmetricModelvars[NLeftRightSymmetricModelvars];
    
    /**
     * @brief LeftRightSymmetricModel constructor
     */
    LeftRightSymmetricModel();
    
    /**
     * @brief LeftRightSymmetricModel destructor
     */
    ~LeftRightSymmetricModel();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    void setModelLeftRightSymmetricModel(){
        flagLeftRightSymmetricModel = true;
    }

    /**
     * @brief A get method to access the member reference of type StandardModelMatching.
     * @return a reference to a StandardModelMatching object
     */
    virtual LeftRightSymmetricModelMatching& getMatching() const
    {
        return LRSMM.getObj();
    }

    virtual void setMatching(LeftRightSymmetricModelMatching& LRSMMr)
    {
        LRSMM.setObj(LRSMMr);
    }


    /*------------------------------------------------------------------------*/
    // Flags

    virtual bool setFlagStr(const std::string name, const std::string value);
    virtual bool setFlag(const std::string, const bool);

    /*------------------------------------------------------------------------*/

    /**
     *
     * @return Flag to switch on CP violation
     */
    bool getCPVflag() const {
//        if(alpha_LRSM==0.0 && delta2_LRSM==0.0 && thetaL==0.0) flag_CPV=false;
        return flag_CPV;
    }

    //Parametrization according to 0712.4218

    /**
     * @brief Gets the first scalar mass square.
     * @return 
     */
    double getmH0sq1() const
    {
        return mH0sq[0];
    }

    /**
     * @brief Gets the second scalar mass square.
     * @return 
     */
    double getmH0sq2() const
    {
        return mH0sq[1];
    }

    /**
     * @brief Gets the third scalar mass square.
     * @return 
     */
    double getmH0sq3() const
    {
        return mH0sq[2];
    }

    /**
     * @brief Gets the fourth scalar mass square.
     * @return 
     */
    double getmH0sq4() const
    {
        return mH0sq[3];
    }

    /**
     * @brief Gets the fifth scalar mass square.
     * @return 
     */
    double getmH0sq5() const
    {
        return mH0sq[4];
    }

    
    /**
     * @brief Gets the approximate first scalar mass square.
     * @return 
     */
    double getmH0sq1_app() const
    {
        return mH0sq_app[0];
    }

    /**
     * @brief Gets the approximate second scalar mass square.
     * @return 
     */
    double getmH0sq2_app() const
    {
        return mH0sq_app[1];
    }

    /**
     * @brief Gets the approximate third scalar mass square.
     * @return 
     */
    double getmH0sq3_app() const
    {
        return mH0sq_app[2];
    }

    /**
     * @brief Gets the approximate fourth scalar mass square.
     * @return 
     */
    double getmH0sq4_app() const
    {
        return mH0sq_app[3];
    }

    
    
    
    
    
    /**
     * @return @f$m_{H_1^+}^2@f$
     */
    double getmH1p_2() const {
        return mH1p_2;
    }

    /**
     * @return @f$m_{H_2^+}^2@f$
     */
    double getmH2p_2() const {
        return mH2p_2;
    }

    /**
     * @return @f$m_{\delta^{++}_R}^2@f$
     */
    double getmdeltappR_2() const {
        return mdeltappR_2;
    }

    /**
     * @return @f$\xi_{\text{LRSM}}@f$
     */
    double getxi_LRSM() const {
        return xi_LRSM;
    }

    /**
     * @return @f$mWR@f$
     */
    double getmWR() const {
        return mWR;
    }

    /**
     * @return @f$\lambda_1@f$
     */
    double getlambda1_LRSM() const {
        return lambda1_LRSM;
    }

    /**
     * @return @f$\lambda_2@f$
     */
    double getlambda2_LRSM() const {
        return lambda2_LRSM;
    }

    /**
     * @return @f$\lambda_3@f$
     */
    double getlambda3_LRSM() const {
        return lambda3_LRSM;
    }

    /**
     * @return @f$\lambda_4@f$
     */
    double getlambda4_LRSM() const {
        return lambda4_LRSM;
    }

    /**
     * @return @f$\rho_1@f$
     */
    double getrho1_LRSM() const {
        return rho1_LRSM;
    }

    /**
     * @return @f$\rho_4@f$
     */
    double getrho4_LRSM() const {
        return rho4_LRSM;
    }

    /**
     * @return @f$\alpha_1@f$
     */
    double getalpha1_LRSM() const {
        return alpha1_LRSM;
    }

    /**
     * @return @f$\alpha_2@f$
     */
    double getalpha2_LRSM() const {
        return alpha2_LRSM;
    }

    /**
     * @return @f$\alpha_@f$
     */
    double getalpha_LRSM() const {
        return alpha_LRSM;
    }

    /**
     * @return @f$\delta_2@f$
     */
    double getdelta2_LRSM() const {
        return delta2_LRSM;
    }

    /**
     * @return @f$\theta^R_{12}@f$
     */
    double getthetaR12_LRSM() const {
        return thetaR12_LRSM;
    }

    /**
     * @return @f$\theta^R_{13}@f$
     */
    double getthetaR13_LRSM() const {
        return thetaR13_LRSM;
    }

    /**
     * @return @f$\theta^R_{23}@f$
     */
    double getthetaR23_LRSM() const {
        return thetaR23_LRSM;
    }

    /**
     * @return @f$\phi^R_1 = \phi^R@f$
     */
    double getphiR1_LRSM() const {
        return phiR1_LRSM;
    }

    /**
     * @return @f$\phi^R_2 = \phi^u_2@f$
     */
    double getphiR2_LRSM() const {
        return phiR2_LRSM;
    }

    /**
     * @return @f$\phi^R_3 = \phi^u_3@f$
     */
    double getphiR3_LRSM() const {
        return phiR3_LRSM;
    }

    /**
     * @return @f$\phi^R_4 = \phi^d_1@f$
     */
    double getphiR4_LRSM() const {
        return phiR4_LRSM;
    }

    /**
     * @return @f$\phi^R_5 = \phi^d_2@f$
     */
    double getphiR5_LRSM() const {
        return phiR5_LRSM;
    }

    /**
     * @return @f$\phi^R_6 = \phi^d_3@f$
     */
    double getphiR6_LRSM() const {
        return phiR6_LRSM;
    }

    /**
     * @brief A get method to retrieve the %CKMR matrix. 
     * @return the %CKMR matrix
     */
    gslpp::matrix<gslpp::complex> getVCKMR() const
    {
        gslpp::matrix<gslpp::complex> VCKMR(3,3,0.);
        double c12 = cos(thetaR12_LRSM);
        double s12 = sin(thetaR12_LRSM);
        double c13 = cos(thetaR13_LRSM);
        double s13 = sin(thetaR13_LRSM);
        double c23 = cos(thetaR23_LRSM);
        double s23 = sin(thetaR23_LRSM);
        gslpp::complex im = gslpp::complex::i();
        VCKMR.assign(0,0, c12*c13*(cos(phiR4_LRSM)-im*sin(phiR4_LRSM)) );
        VCKMR.assign(0,1, s12*c13*(cos(phiR5_LRSM)-im*sin(phiR5_LRSM)) );
        VCKMR.assign(0,2, s13*(cos(phiR4_LRSM+phiR6_LRSM)-im*sin(phiR4_LRSM+phiR6_LRSM)) );
        VCKMR.assign(1,0, -(s12*c23+c12*s13*s23*(cos(phiR1_LRSM)+im*sin(phiR1_LRSM)))*(cos(phiR2_LRSM-phiR4_LRSM)+im*sin(phiR2_LRSM-phiR4_LRSM)) );
        VCKMR.assign(1,1, (c12*c23-s12*s13*s23*(cos(phiR1_LRSM)+im*sin(phiR1_LRSM)))*(cos(phiR2_LRSM-phiR5_LRSM)+im*sin(phiR2_LRSM-phiR5_LRSM)) );
        VCKMR.assign(1,2, c13*s23*(cos(phiR2_LRSM-phiR6_LRSM)+im*sin(phiR2_LRSM-phiR6_LRSM)) );
        VCKMR.assign(2,0, (s12*s23-c12*s13*c23*(cos(phiR1_LRSM)+im*sin(phiR1_LRSM)))*(cos(phiR3_LRSM-phiR4_LRSM)+im*sin(phiR3_LRSM-phiR4_LRSM)) );
        VCKMR.assign(2,1, -(c12*s23+s12*s13*c23*(cos(phiR1_LRSM)+im*sin(phiR1_LRSM)))*(cos(phiR3_LRSM-phiR5_LRSM)+im*sin(phiR3_LRSM-phiR5_LRSM)) );
        VCKMR.assign(2,2, c13*c23*(cos(phiR3_LRSM-phiR6_LRSM)+im*sin(phiR3_LRSM-phiR6_LRSM)) );
        return VCKMR;
    }

    /**
     * @return LRSM scale
     */
    double getQ_LRSM() const {
        return Q_LRSM;
    }

protected: 

    virtual void setParameter(const std::string, const double&);
    
    LRSMquantities* myLRSMquantities;

    double mH0sq[5];
    double mH0sq_app[4];
    gslpp::matrix<gslpp::complex> U;

    mutable Matching<LeftRightSymmetricModelMatching,LeftRightSymmetricModel> LRSMM; ///< An object of type Matching.

private:

    double mH1p_2, mH2p_2, mdeltappR_2, xi_LRSM, mWR;
    double lambda1_LRSM, lambda2_LRSM, lambda3_LRSM, lambda4_LRSM, rho1_LRSM, rho4_LRSM, alpha1_LRSM, alpha2_LRSM, alpha_LRSM, delta2_LRSM, Q_LRSM;
    double thetaR12_LRSM, thetaR13_LRSM, thetaR23_LRSM, phiR1_LRSM, phiR2_LRSM, phiR3_LRSM, phiR4_LRSM, phiR5_LRSM, phiR6_LRSM;
    bool flagLeftRightSymmetricModel;
    bool flag_CPV;

};

#endif	/* LEFTRIGHTSYMMETRICMODEL_H */
