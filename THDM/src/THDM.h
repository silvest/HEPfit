/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDM_H
#define	THDM_H

#include <StandardModel.h>
#include <PVfunctions.h>
#include "THDMcache.h"


using namespace gslpp;
//class THDMcache; //forward reference to THDMcache class

/**
 * @class THDM
 * @brief  Two Higgs Doublet model type II
 */
class THDM: public StandardModel {
public:

    static const int NTHDMvars = 8;
    //The parameters of the Higgs potential for THDM (there are several basis)Is it right? 
    //I choose the physical basis: mHp,sin_ba,lambda6,lambda7,mA,m12_2,tanb,mH
    static const std::string THDMvars[NTHDMvars];
    
    /**
     * @brief SUSY constructor
     */
    THDM();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual std::string ModelName() const {
        return "THDM";
    }
    
    ///////////////////////////////////////////////////////////////////////////
    
    /**
     * 
     * @return the down-type VEV
     */
    double v1();

    /**
     *
     * @return the up-type VEV
     */
    double v2();
    
    ///////////////////////////////////////////////////////////////////////////

//    THDMcache GetMycache() const {
//        return mycache;
//    }

    /**
     *
     * @return cos beta
     */
    double GetCosb() const {
        return cosb;
    }
    
    /**
     * @brief set tan beta, sin beta and cos beta
     * @param cosb cos beta
     */
    void SetCosb(double cosb);

    /**
     *
     * @return charged Higgs mass
     */
    double GetMHp() const {
        return mHp;
    }

    
    /**
     * @brief set the charged Higgs mass
     * @param mHp charged Higgs mass
     */
    void SetMHp(double mHp) {
        this->mHp = mHp;
    }
    
    
    /**
     *
     * @return sin beta
     */
    double GetSinb() const {
        return sinb;
    }

    /**
     *
     * @return \f$\tan \beta$\f
     */
    double GetTanb() const {
        return tanb;
    }

    /**
     * @brief set tan beta, sin beta and cos beta
     * @param tanb tan beta
     */
    void SetTanb(double tanb) {
        this->tanb = tanb;
    }

    
    /**
     *
     * @return Higgs potential parameter \f$\lambda6$\f (reference arXiv:0902.0851v2)  
     */
    double GetLambda6() const {
        return lambda6;
    }

    /**
     * @brief set lambda6
     * @param lambda6 lambda6 
     */
    void SetLambda6(double lambda6) {
        this->lambda6 = lambda6;
    }

    /**
     *
     * @return \f$lambda_7$\f 
     */
    double GetLambda7() const {
        return lambda7;
    }

    /**
     * @brief set lambda7
     * @param lambda7 lambda7 
     */
    void SetLambda7(double lambda7) {
        this->lambda7 = lambda7;
    }

    /**
     *
     * @return parameter of the Higgs potential \f$m_{12}^2$\f 
     */
    double GetM12_2() const {
        return m12_2;
    }

    /**
     * @brief set \f$m_{12}^2$\f
     * @param m12_2 \f$m_{12}^2$\f
     */
    void SetM12_2(double m12_2) {
        this->m12_2 = m12_2;
    }

    /**
     *
     * @return mass pseudoscalar Higgs A 
     */
    double GetMA() const {
        return mA;
    }

    /**
     * @brief set mA
     * @param mA \f$m_A$\f 
     */
    void SetMA(double mA) {
        this->mA = mA;
    }

    /**
     *
     * @return mass heavy neutral scalar Higgs 
     */
    double GetMH() const {
        return mH;
    }

    /**
     * @brief set mH
     * @param mH \f$m_H$\f 
     */
    void SetMH(double mH) {
        this->mH = mH;
    }

    /**
     *
     * @return \f$\sin(\beta-\alpha)$\f 
     */
    double GetSin_ba() const {
        return sin_ba;
    }

    /**
     * @brief set sin_ba
     * @param sin_ba \f$\sin(\beta-\alpha)$\f 
     */
    void SetSin_ba(double sin_ba) {
        this->sin_ba = sin_ba;
    }
    
    
    ///////////////////////////////////////////////////////////////////////////
    /* Functions for EW precision observables */

    double obliqueS() const;
    double obliqueT() const;
    double obliqueU() const;
    
    
    ///////////////////////////////////////////////////////////////////////////

protected: 
    
    virtual void SetParameter(const std::string, const double&);
    THDMcache mycache;

    private:
        
    //PVfunctions PV;
    
    //double DeltaS, DeltaT, DeltaU;
    
    double mH, mA, mHp, tanb, sinb, cosb, sin_ba, m12_2, lambda6, lambda7, mh;
    
    //double Mz2;
    //double s_W2;//\sin^2(\theta_W)   
    //double Mw_i, Mw2;
    //double cos2_ba, sin2_ba;   
    
    bool computeCKM, computeYe, computeYn;
    
    ////////////////////////////////////////////////////////////////////////////
    /*One-loop functions*/
    
    /**
     * @brief function F(m0,m1) used for THDM. Remember that this function is
     * defined for THDM while for SUSY we have a multiplicative factor 2.
     * @param[in] m0 mass m_0
     * @param[in] m1 mass m_1
     * @return the function F for THDM 
     */
    double F(const double m0, const double m1) const;
    
//    complex B0_Mz_0_Mz_mH;
//    complex B0_Mz_0_Mz_mh;
//    complex B0_Mz_0_Mw_mH;
//    complex B0_Mz_0_Mw_mh;
//    
//    complex B22prime_Mz_Mz2_mH_mA;
//    complex B22prime_Mz_Mz2_mHp_mHp;
//    complex B22prime_Mz_Mz2_mh_mA;
//    complex B22prime_Mz_Mz2_Mz_mH;
//    complex B22prime_Mz_Mz2_Mz_mh;
//    complex B22prime_Mz_Mw2_mA_mHp;
//    complex B22prime_Mz_Mw2_mHp_mHp;
//    complex B22prime_Mz_Mw2_mH_mHp;
//    complex B22prime_Mz_Mw2_mh_mHp;
//    complex B22prime_Mz_Mw2_Mw_mH;
//    complex B22prime_Mz_Mw2_Mw_mh;
//    
//    complex B0prime_Mz_Mz2_Mz_mH;
//    complex B0prime_Mz_Mz2_Mz_mh;
//    complex B0prime_Mz_Mw2_Mw_mH;
//    complex B0prime_Mz_Mw2_Mw_mh;
    
    
};

#endif	/* THDM_H */

