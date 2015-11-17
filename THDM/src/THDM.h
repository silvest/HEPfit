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

class THDMcache; //forward reference to THDMcache class

/**
 * @class THDM
 * @brief Z2 symmetric Two Higgs Doublet models
 */
class THDM: public StandardModel {
public:

    static const int NTHDMvars = 18;
    //The parameters of the Higgs potential for THDM (there are several basis) according to ?
    //We choose a modified physical basis: log(tanbeta),beta-alpha,mHh^2,mHh^2-mA^2,mHh^2-mHp^2,m12^2,lambda6,lambda7
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

    
    ///////////////////////////////////////////////////////////////////////////
    // Flags

    virtual bool setFlagStr(const std::string name, const std::string value);

    virtual THDMMatching* getMyMatching() const
    {
        return myTHDMMatching;
    }

    
    ///////////////////////////////////////////////////////////////////////////
    
    /**
     * 
     * @return the VEV \f$v_1$\f
     */
    double getv1() const {
        return v() * cosb;
    }

    /**
     *
     * @return the VEV \f$v_2$\f
     */
    double getv2() const {
        return v() * sinb;
    }
    
    ///////////////////////////////////////////////////////////////////////////

//    THDMcache GetMycache() const {
//        return mycache;
//    }

    /**
     *
     * @return \f$\log(\tan \beta)$\f
     */
    std::string getModelTypeflag() const {
        return flag_model;
    }

    /**
     *
     * @return \f$\log(\tan \beta)$\f
     */
    double getlogtb() const {
        return logtb;
    }

    /**
     *
     * @return \f$\tan \beta$\f
     */
    double gettanb() const {
        return tanb;
    }

    /**
     *
     * @return \f$\sin \beta$\f
     */
    double getsinb() const {
        return sinb;
    }

    /**
     *
     * @return \f$\cos \beta$\f
     */
    double getcosb() const {
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

    /**
     *
     * @return \f$\cos \alpha$\f
     */
    double getcosa() const{
        return cos(atan(pow(10.,logtb))-bma);
    }

    /**
     *
     * @return \f$\sin \alpha$\f
     */
    double getsina() const{
        return sin(atan(pow(10.,logtb))-bma);
    }

    /**
     *
     * @return mass squared heavy neutral scalar Higgs
     */
    double getmHh2() const {
        return mHh2;
    }

    /**
     *
     * @return mass heavy neutral scalar Higgs
     */
    double getmHh() const {
        return sqrt(mHh2);
    }

    /**
     *
     * @return mass squared pseudoscalar Higgs A
     */
    double getmA2() const {
        return mHh2-mHh2mmA2;
    }

    /**
     *
     * @return mass pseudoscalar Higgs A
     */
    double getmA() const {
        return sqrt(mHh2-mHh2mmA2);
    }

    /**
     *
     * @return charged Higgs mass squared
     */
    double getmHp2() const {
        return mHh2-mHh2mmHp2;
    }

    /**
     *
     * @return charged Higgs mass
     */
    double getmHp() const {
        return sqrt(mHh2-mHh2mmHp2);
    }

    /**
     *
     * @return parameter of the Higgs potential \f$m_{12}^2$\f 
     */
    double getm12_2() const {
        return m12_2;
    }

    /**
     *
     * @return Higgs potential parameter \f$\lambda6$\f (reference arXiv:0902.0851v2)  
     */
    double getlambda6() const {
        return lambda6;
    }

    /**
     *
     * @return \f$lambda_7$\f 
     */
    double getlambda7() const {
        return lambda7;
    }

    /**
     *
     * @return BDtaunu SM expectation
     */
    double getBDtaunu_SM() const {
        return BDtaunu_SM;
    }

    /**
     *
     * @return BDtaunu coefficient A
     */
    double getBDtaunu_A() const {
        return BDtaunu_A;
    }

    /**
     *
     * @return BDtaunu coefficient B
     */
    double getBDtaunu_B() const {
        return BDtaunu_B;
    }

    /**
     *
     * @return BDstartaunu SM expectation
     */
    double getBDstartaunu_SM() const {
        return BDstartaunu_SM;
    }

    /**
     *
     * @return BDstartaunu coefficient A
     */
    double getBDstartaunu_A() const {
        return BDstartaunu_A;
    }

    /**
     *
     * @return BDtaunu coefficient B
     */
    double getBDstartaunu_B() const {
        return BDstartaunu_B;
    }

    /**
     *
     * @return BHatBsTHDM for DeltamBs
     */
    double getBHatBsTHDM() const {
        return BHatBsTHDM;
    }

    /**
     *
     * @return etaBsTHDM for DeltamBs
     */
    double getetaBsTHDM() const {
        return etaBsTHDM;
    }

    /**
     *
     * @return nuisance parameter for the theoretical error on bsgamma
     */
    double getbsgamma_theoryerror() const {
        return bsgamma_theoryerror;
    }

    /**
     *
     * @return THDM scale
     */
    double getQ_THDM() const {
        return Q_THDM;
    }

protected: 

    virtual void setParameter(const std::string, const double&);
    THDMcache * mycache;

private:

    THDMMatching* myTHDMMatching;

    double logtb, tanb, sinb, cosb, bma, sin_ba, mHh2, mHh2mmA2, mHh2mmHp2, m12_2, lambda6, lambda7/*, mh*/;
    double BDtaunu_SM, BDtaunu_A, BDtaunu_B, BDstartaunu_SM, BDstartaunu_A, BDstartaunu_B, BHatBsTHDM, etaBsTHDM, bsgamma_theoryerror, Q_THDM;
    std::string flag_model;
};

#endif	/* THDM_H */

