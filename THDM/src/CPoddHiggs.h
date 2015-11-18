/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CPODDHIGGSCACHE_H
#define	CPODDHIGGSCACHE_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"
#include "THDMfunctions.h"
#include "THDMcache.h"
#include "lightHiggs.h"

/**
 * @addtogroup THDM
 * @brief A module for the @f$Z_2@f$ symmetric Two-Higgs-Doublet models.
 * @{
 */

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
    
protected:
    
    THDMcache * mycache;
    lightHiggs * mylightHiggs;

    double ggF_A_tautau_TH;
    double bbF_A_tautau_TH;
    double ggF_A_gaga_TH;
    double ggF_A_hZ_bbll_TH;
    double ggF_A_hZ_bbZ_TH;
    double ggF_A_hZ_tautaull_TH;
    double ggF_A_hZ_tautauZ_TH;   
    double pp_A_tt_TH;
    double bbF_A_bb_TH;
     
    double ggF_A_tautau_EX_ATLAS;
    double ggF_A_tautau_EX_CMS;
    double bbF_A_tautau_EX_ATLAS;
    double bbF_A_tautau_EX_CMS;
    double ggF_A_gaga_EX_ATLAS;
    double ggF_A_gaga_EX_CMS;
    double ggF_A_hZ_bbll_EX_CMS;
    double ggF_A_hZ_tautauZ_EX_ATLAS;
    double ggF_A_hZ_tautaull_EX_CMS;
    double ggF_A_hZ_bbZ_EX_ATLAS;
    double pp_A_tt_EX_ATLAS;
    double bbF_A_bb_EX_CMS;     

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
    
    double TAUc;
    double TAUt;
    double TAUs;
    double TAUb;
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

    double LAMc;
    double LAMt;
    double LAMs;
    double LAMb;
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
    double Br_Ztoee; 
    double Br_Ztomumu;
    double Br_Ztotautau;
    
    gslpp::complex I_HH_f;
    gslpp::complex I_HH_fU;
    gslpp::complex I_HH_fD;
    gslpp::complex I_HH_fL;
    gslpp::complex I_HH_W;

    int HSTheta (const double x) const;
    double KaellenFunction (const double a, const double b, const double c) const;


};



class  Hobs_ggF_A_tautau_ATLAS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_A_tautau_ATLAS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_ggF_A_tautau_ATLAS
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_ggF_A_tautau_CMS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_A_tautau_CMS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_ggF_A_tautau
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_bbF_A_tautau_ATLAS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_bbF_A_tautau_ATLAS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_bbF_A_tautau_ATLAS
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_bbF_A_tautau_CMS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_bbF_A_tautau_CMS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_bbF_A_tautau
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_ggF_A_gaga_ATLAS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_A_gaga_ATLAS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_ggF_A_gaga_ATLAS
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_ggF_A_gaga_CMS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_A_gaga_CMS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_ggF_A_gaga_CMS
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_ggF_A_hZ_bbll_CMS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_A_hZ_bbll_CMS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_hZ
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_ggF_A_hZ_bbZ_ATLAS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_A_hZ_bbZ_ATLAS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_hZ_bbZ
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_ggF_A_hZ_tautaull_CMS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_A_hZ_tautaull_CMS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_hZ
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_ggF_A_hZ_tautauZ_ATLAS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_A_hZ_tautauZ_ATLAS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_A_hZ_tautauZ
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_pp_A_tt_ATLAS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_pp_A_tt_ATLAS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_pp_A_tt
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_bbF_A_bb_CMS: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_bbF_A_bb_CMS(const StandardModel& SM_i);
    
    /**
     * @return Hobs_bbF_A_bb
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



class  log10_ggF_A_gaga_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_ggF_A_gaga_TH(const StandardModel& SM_i);
    
    /**
     * @return ggF_A_gaga_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_ggF_A_hZ_bbll_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_ggF_A_hZ_bbll_TH(const StandardModel& SM_i);
    
    /**
     * @return A_hZ_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_ggF_A_hZ_bbZ_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_ggF_A_hZ_bbZ_TH(const StandardModel& SM_i);
    
    /**
     * @return A_hZ_bbZ_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_ggF_A_hZ_tautaull_TH: public CPoddHiggsCache {
 public:

  /**                                                                                                                                         
   * @brief Constructor.                                                                                                                      
   */
  log10_ggF_A_hZ_tautaull_TH(const StandardModel& SM_i);

  /**                                                                                                                                         
   * @return A_bb_TH                                                                                                                          
   */
  double computeThValue ();

 private:

};



class  log10_ggF_A_hZ_tautauZ_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_ggF_A_hZ_tautauZ_TH(const StandardModel& SM_i);
    
    /**
     * @return A_hZ_tautauZ_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_pp_A_tt_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_pp_A_tt_TH(const StandardModel& SM_i);
    
    /**
     * @return pp_A_tt_TH
     */
    double computeThValue ();
    
private:
    
};



class  log10_bbF_A_bb_TH: public CPoddHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_bbF_A_bb_TH(const StandardModel& SM_i);
    
    /**
     * @return bbF_A_bb_TH
     */
    double computeThValue ();
    
private:

};

/**
 * @}
 */

#endif	/* CPODDHIGGSCACHE_H */
