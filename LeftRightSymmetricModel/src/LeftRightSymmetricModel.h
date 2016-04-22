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

    static const int NLeftRightSymmetricModelvars = 18;
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
     * @return @f$\log_{10}(\tan \beta)@f$
     */
    double getlogtb() const {
        return logtb;
    }

    /**
     * @return @f$\tan \beta@f$
     */
    double gettanb() const {
        return tanb;
    }

    /**
     * @return @f$\sin \beta@f$
     */
    double getsinb() const {
        return sinb;
    }

    /**
     * @return @f$\cos \beta@f$
     */
    double getcosb() const {
        return cosb;
    }

    /**
     * @return mass squared charged Higgs
     */
    double getmHp2() const {
        return mHp2;
    }

    /**
     * @return mass of the charged Higgs
     */
    double getmHp() const {
    if(mHp2 < 0.) {
        return 0.;
    }
        return sqrt(mHp2);
    }

    /**
     * @return alpha1
     */
    double getalpha1() const {
        return alpha1;
    }

    /**
     * @return alpha2
     */
    double getalpha2() const {
        return alpha2;
    }

    /**
     * @return alpha3
     */
    double getalpha3() const {
        return alpha3;
    }

    /**
     * @return m2_2
     */
    double getm2_2() const {
        return m2_2;
    }

    /**
     * @return parameter of the Higgs potential @f$Re(m_{12}^2)@f$ 
     */
    double getRem12_2() const {
        return Rem12_2;
    }

    /**
     * @return Relambda6
     */
    double getRelambda6() const {
        return Relambda6;
    }

    /**
     * @return Imlambda6
     */
    double getImlambda6() const {
        return Imlambda6;
    }

    /**
     * @return Relambda7
     */
    double getRelambda7() const {
        return Relambda7;
    }

    /**
     * @return Imlambda7
     */
    double getImlambda7() const {
        return Imlambda7;
    }

    /**
     * @return abssigma_u
     */
    double getabssigma_u() const {
        return abssigma_u;
    }

    /**
     * @return phi_u
     */
    double getphi_u() const {
        return phi_u;
    }

    /**
     * @return abssigma_d
     */
    double getabssigma_d() const {
        return abssigma_d;
    }

    /**
     * @return phi_d
     */
    double getphi_d() const {
        return phi_d;
    }

    /**
     * @return abssigma_l
     */
    double getabssigma_l() const {
        return abssigma_l;
    }

    /**
     * @return phi_l
     */
    double getphi_l() const {
        return phi_l;
    }

    /**
     * @return General THDM scale
     */
    double getQ_GTHDM() const {
        return Q_GTHDM;
    }

protected: 

    virtual void setParameter(const std::string, const double&);
    

//    mutable Matching<LeftRightSymmetricModelMatching,LeftRightSymmetricModel> LRSMM; ///< An object of type Matching.

private:

    double logtb, tanb, sinb, cosb, mHp2, alpha1, alpha2, alpha3, m2_2, Rem12_2, Relambda6, Imlambda6, Relambda7, Imlambda7, abssigma_u, phi_u, abssigma_d, phi_d, abssigma_l, phi_l, Q_GTHDM;
    bool flagLeftRightSymmetricModel;
};

#endif	/* LEFTRIGHTSYMMETRICMODEL_H */
