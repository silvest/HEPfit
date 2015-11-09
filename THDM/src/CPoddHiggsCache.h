/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CPODDHIGGSCACHE_H
#define	CPODDHIGGSCACHE_H

#include <stdexcept>
//#include <gslpp.h>
#include <ThObservable.h>
#include "THDM.h"
#include "THDMfunctions.h"
#include "THDMcache.h"
#include "lightHiggs.h"

/**
 * @class CPoddHiggsCache
 * @brief .
 */
class CPoddHiggsCache : public ThObservable {
public:
    CPoddHiggsCache(const StandardModel& SM_i);
    virtual ~CPoddHiggsCache();
    void computeParameters();
    
    double computeThValue();

    double cos_ab;
    
protected:
    
    THDMcache * mycache;
    lightHiggs * mylightHiggs;

    double ggF_A_tautau_TH;
    double bbF_A_tautau_TH;
    double A_gaga_TH;
    double A_hZ_TH; 
    double A_hZ_tautauZ_TH;  
    double A_hZ_bbZ_TH;
    double A_tt_TH;
    double A_bb_TH;
    
    double ggF_A_tautau_EX;
    double bbF_A_tautau_EX;
    double A_gaga_EX;
    double A_hZ_EX;
    double A_hZ_tautauZ_EX;
    double A_hZ_bbZ_EX;
    double A_tt_EX;
    double A_bb_EX;  

private:
    const THDM * myTHDM;
    const StandardModel& mySM;
    
    double SigmaggF;
    double Sigmaggh_tt;
    double Sigmaggh_bb;
    
    double rA_QuQu; 
    double rA_QdQd;
    double rA_ll;
    double rA_gg; 
   
    double Gamma_Agaga; 
    double Gamma_AZga;
    double Gamma_Agg;
    
    double TAUu;
    double TAUc;
    double TAUt;
    double TAUd;
    double TAUs;
    double TAUb;
    double TAUe;
    double TAUmu;
    double TAUtau;
    double TAUw;
    double TAUhp;
    
    gslpp::complex I_A_f;
    gslpp::complex I_A_fU;
    gslpp::complex I_A_fD;
    gslpp::complex I_A_fL;
    gslpp::complex I_A_W;
    double g_A_HpHm;
    gslpp::complex I_A_Hp;
    
    double LAMu;
    double LAMc;
    double LAMt;
    double LAMd;
    double LAMs;
    double LAMb;
    double LAMe;
    double LAMmu;
    double LAMtau;
    double LAMw;
    double LAMhp;
    
    gslpp::complex A_A_F;
    gslpp::complex A_A_U;
    gslpp::complex A_A_D;
    gslpp::complex A_A_L;
    gslpp::complex A_A_W;
    gslpp::complex A_A_Hp;
    
    int modelType;
    
    double SigmaggF_A;
    double SigmabbF_A;
    
    double SigmaSum;
    
    double BrSM_Atocc;
    double BrSM_Atobb;
    double BrSM_Atott;
    double BrSM_Atomumu;
    double BrSM_Atotautau;
    
    double GammaAtotSM;
    double GammaAHZ;
    double GammaAhZ;
    double GammaAHpW;
    double GammaAtot;
    
    double Br_Atott;
    double Br_Atobb;
    double Br_Atotautau;
    double Br_Atogaga;
    double Br_AtohZ;
    double Br_htotautau;
    double Br_htobb;
    
    gslpp::complex I_HH_f;
    gslpp::complex I_HH_fU;
    gslpp::complex I_HH_fD;
    gslpp::complex I_HH_fL;
    gslpp::complex I_HH_W;

    int HSTheta (const double x) const;
    double KaellenFunction (const double a, const double b, const double c) const;


};



class  Hobs_ggF_A_tautau: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_A_tautau(const StandardModel& SM_i);
    
    /**
     * @return Hobs_ggF_A_tautau
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_bbF_A_tautau: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_bbF_A_tautau(const StandardModel& SM_i);
    
    /**
     * @return Hobs_bbF_A_tautau
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_A_gaga: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_A_gaga(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_gaga
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_A_hZ: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_A_hZ(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_hZ
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_A_hZ_tautauZ: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_A_hZ_tautauZ(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_hZ_tautauZ
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_A_hZ_bbZ: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_A_hZ_bbZ(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_hZ_bbZ
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_A_tt: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_A_tt(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_tt
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_A_bb: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_A_bb(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_bb
     */
    double computeThValue ();
    
private:
    
};



class  log10_ggF_A_tautau_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_ggF_A_tautau_TH(const StandardModel& SM_i);
    
    /**
     * @return ggF_A_tautau_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_bbF_A_tautau_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_bbF_A_tautau_TH(const StandardModel& SM_i);
    
    /**
     * @return bbF_A_tautau_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_A_gaga_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_A_gaga_TH(const StandardModel& SM_i);
    
    /**
     * @return A_gaga_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_A_hZ_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_A_hZ_TH(const StandardModel& SM_i);
    
    /**
     * @return A_hZ_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_A_hZ_tautauZ_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_A_hZ_tautauZ_TH(const StandardModel& SM_i);
    
    /**
     * @return A_hZ_tautauZ_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_A_hZ_bbZ_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_A_hZ_bbZ_TH(const StandardModel& SM_i);
    
    /**
     * @return A_hZ_bbZ_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_A_tt_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_A_tt_TH(const StandardModel& SM_i);
    
    /**
     * @return A_tt_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_A_bb_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_A_bb_TH(const StandardModel& SM_i);
    
    /**
     * @return A_bb_TH
     */
    double computeThValue ();
    
private:

};



#endif	/* CPODDHIGGSCACHE_H */