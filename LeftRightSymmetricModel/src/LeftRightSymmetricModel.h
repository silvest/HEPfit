/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEFTRIGHTSYMMETRICMODEL_H
#define	LEFTRIGHTSYMMETRICMODEL_H

#include "StandardModel.h"
//#include "LeftRightSymmetricModelMatching.h"

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

    static const int NLeftRightSymmetricModelvars = 14;
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

//    /**
//     * @brief A get method to access the member reference of type StandardModelMatching.
//     * @return a reference to a StandardModelMatching object
//     */
//    virtual LeftRightSymmetricModelMatching& getMatching() const
//    {
//        return LRSMM.getObj();
//    }

//    virtual void setMatching(LeftRightSymmetricModelMatching& LRSMMr)
//    {
//        LRSMM.setObj(LRSMMr);
//    }


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
     * @brief Gets the first scalar mass square.
     * @return 
     */
    double getmH0sq1_app() const
    {
        return mH0sq_app[0];
    }

    /**
     * @brief Gets the second scalar mass square.
     * @return 
     */
    double getmH0sq2_app() const
    {
        return mH0sq_app[1];
    }

    /**
     * @brief Gets the third scalar mass square.
     * @return 
     */
    double getmH0sq3_app() const
    {
        return mH0sq_app[2];
    }

    /**
     * @brief Gets the fourth scalar mass square.
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

//    mutable Matching<LeftRightSymmetricModelMatching,LeftRightSymmetricModel> LRSMM; ///< An object of type Matching.

private:

    double mH1p_2, mH2p_2, mdeltappR_2, xi_LRSM, mWR, lambda1_LRSM, lambda2_LRSM, lambda3_LRSM, lambda4_LRSM, rho1_LRSM, rho4_LRSM, alpha1_LRSM, alpha2_LRSM, Q_LRSM;
    bool flagLeftRightSymmetricModel;
    bool flag_CPV;
};

#endif	/* LEFTRIGHTSYMMETRICMODEL_H */
