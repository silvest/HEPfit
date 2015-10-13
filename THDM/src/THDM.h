/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDM_H
#define	THDM_H

#include <StandardModel.h>
#include "THDMMatching.h"
#include "THDMcache.h"


using namespace gslpp;
class THDMcache; //forward reference to THDMcache class

/**
 * @class THDM
 * @brief Two Higgs Doublet model type II
 */
class THDM: public StandardModel {
public:

    static const int NTHDMvars = 9;
    //The parameters of the Higgs potential for THDM (there are several basis) according to ?
    //We choose the physical basis: logtb,bma,mHh,mA,mHp,m12_2,lambda6,lambda7
    static const std::string THDMvars[NTHDMvars];
    
    /**
     * @brief THDM constructor
     */
    THDM();
    
    /**
     * @brief THDM destructor
     */
    ~THDM();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual bool setFlag(const std::string name, const bool value);
    
    virtual THDMMatching* getMyMatching() const
    {
        return myTHDMMatching;
    }

    
    ///////////////////////////////////////////////////////////////////////////
    
    /**
     * 
     * @return the down-type VEV
     */
    double v1() const;

    /**
     *
     * @return the up-type VEV
     */
    double v2() const;
    
    ///////////////////////////////////////////////////////////////////////////

//    THDMcache GetMycache() const {
//        return mycache;
//    }

    /**
     *
     * @return \f$\log(\tan \beta)$\f
     */
    double getLogtb() const {
        return logtb;
    }

    /**
     *
     * @return \f$\tan \beta$\f
     */
    double getTanb() const {
        return tanb;
    }

    /**
     *
     * @return \f$\sin \beta$\f
     */
    double getSinb() const {
        return sinb;
    }

    /**
     *
     * @return \f$\cos \beta$\f
     */
    double getCosb() const {
        return cosb;
    }

    /**
     *
     * @return \f$\beta-\alpha$\f
     */
    double getbma() const {
        return bma;
    }

   /**
     *
     * @return \f$\sin(\beta-\alpha)$\f
     */
    double getsin_ba() const {
        return sin_ba;
    }

    double computeCosa() const;

    double computeSina() const;
 
    /**
     *
     * @return mass heavy neutral scalar Higgs
     */
    double getMHh() const {
        return mHh;
    }

    /**
     *
     * @return mass pseudoscalar Higgs A
     */
    double getMA() const {
        return mA;
    }

    /**
     *
     * @return charged Higgs mass
     */
    double getMHp() const {
        return mHp;
    }

    /**
     *
     * @return parameter of the Higgs potential \f$m_{12}^2$\f 
     */
    double getM12_2() const {
        return m12_2;
    }

    /**
     *
     * @return Higgs potential parameter \f$\lambda6$\f (reference arXiv:0902.0851v2)  
     */
    double getLambda6() const {
        return lambda6;
    }

    /**
     *
     * @return \f$lambda_7$\f 
     */
    double getLambda7() const {
        return lambda7;
    }
    
    /**
     *
     * @return modelType
     */
    int getModelType() const {
        return modelType;
    }

protected: 
    
    virtual void setParameter(const std::string, const double&);
    THDMcache * mycache;

private:

    THDMMatching* myTHDMMatching;
    
    double logtb, tanb, sinb, cosb, bma, sin_ba, mHh, mA, mHp, m12_2, lambda6, lambda7, mh, modelType;
    
      
};

#endif	/* THDM_H */

