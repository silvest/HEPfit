/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDM_H
#define	GENERALTHDM_H

#include "StandardModel.h"
//#include "GeneralTHDMMatching.h"

//class GeneralTHDMcache; //forward reference to GeneralTHDMcache class

/**
 * @class GeneralTHDM
 * @ingroup GeneralTHDM
 * @brief A base class for generic Two-Higgs-Doublet models.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 *
 * 
 * @anchor GeneralTHDMParameters
 * <h3>%Model parameters</h3>
 */
class GeneralTHDM: public StandardModel {
public:

    static const int NGeneralTHDMvars = 18;
    static const std::string GeneralTHDMvars[NGeneralTHDMvars];
    
    /**
     * @brief GeneralTHDM constructor
     */
    GeneralTHDM();
    
    /**
     * @brief GeneralTHDM destructor
     */
    ~GeneralTHDM();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

//    /**
//     * @brief A get method to access the member reference of type StandardModelMatching.
//     * @return a reference to a StandardModelMatching object
//     */
//    virtual THDMMatching& getMatching() const
//    {
//        return THDMM.getObj();
//    }

//    virtual void setMatching(THDMMatching& THDMMr)
//    {
//        THDMM.setObj(THDMMr);
//    }


    ///////////////////////////////////////////////////////////////////////////
    // Flags

//    virtual bool setFlagStr(const std::string name, const std::string value);
//    virtual bool setFlag(const std::string, const bool);
//
//    THDMcache* getMyTHDMCache() const
//    {
//        return myTHDMcache;
//    }

    ///////////////////////////////////////////////////////////////////////////

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

//    /**
//     * @brief A method to check if the model type name in string form is valid.
//     * @param[in] THDM model type name
//     * @return a boolean that is true if the model type name is valid
//     */
//    bool checkmodelType(const std::string modeltype) const
//    {
//        if (modeltype.compare("type1") == 0
//                || modeltype.compare("type2") == 0
//                || modeltype.compare("typeX") == 0
//                || modeltype.compare("typeY") == 0)
//            return true;
//        else
//            return false;
//    }

//    mutable Matching<THDMMatching,THDM> THDMM; ///< An object of type Matching.

private:

//    THDMcache* myTHDMcache;

    double logtb, tanb, sinb, cosb, mHp2, alpha1, alpha2, alpha3, m2_2, Rem12_2, Relambda6, Imlambda6, Relambda7, Imlambda7, abssigma_u, phi_u, abssigma_d, phi_d, abssigma_l, phi_l, Q_GTHDM;
};

#endif	/* GENERALTHDM_H */
