/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEAVYHIGGSCACHE_H
#define	HEAVYHIGGSCACHE_H

#include <stdexcept>
//#include <gslpp.h>
#include <ThObservable.h>
#include "THDM.h"
#include "THDMfunctions.h"
#include "THDMcache.h"
#include "lightHiggs.h"

/**
 * @class heavyHiggsCache
 * @brief .
 */
class heavyHiggsCache : public ThObservable {
public:
    heavyHiggsCache(const StandardModel& SM_i);
    virtual ~heavyHiggsCache();
    void computeParameters();
    
    double computeThValue();
    
    double cos_2b;
    double cos_ab;
  
protected:

    THDMcache * mycache;
    lightHiggs * mylightHiggs;
    
    double ggF_H_tautau_TH;
    double bbF_H_tautau_TH;
    double H_gaga_TH;  
    double H_ZZ_TH; 
    double ggF_H_WW_TH;
    double H_WW_TH;  
    double VBF_H_WW_TH;
    double H_hh_TH;  
    double H_hh_bbbb_TH;
    double H_tt_TH;
    double H_bb_TH;
    double H_hh_gagabb_TH;
    
    double ggF_H_tautau_EX;
    double bbF_H_tautau_EX;
    double H_gaga_EX;
    double H_ZZ_EX;
    double ggF_H_WW_EX;
    double H_WW_EX;
    double VBF_H_WW_EX;
    double H_hh_EX;
    double H_hh_bbbb_EX;
    double H_tt_EX;
    double H_bb_EX;
    double H_hh_gagabb_EX;
    
private:
    const THDM * myTHDM;
    const StandardModel& mySM;
    
    double SigmaggF;
    double Sigmaggh_tt;
    double Sigmaggh_bb;
    
    double rHH_QuQu; 
    double rHH_QdQd;
    double rHH_ll;
    double rHH_gg; 
    double rHH_VV;
    
    double Gamma_Hgaga;
    double Gamma_HZga;
    double Gamma_Hgg;
    
//    double TAUu;
    double TAUc;
    double TAUt;
//    double TAUd;
    double TAUs;
    double TAUb;
//    double TAUe;
    double TAUmu;
    double TAUtau;
    double TAUw;
    double TAUhp;
    
    gslpp::complex I_HH_f;
    gslpp::complex I_HH_fU;
    gslpp::complex I_HH_fD;
    gslpp::complex I_HH_fL;
    gslpp::complex I_HH_W;
    double g_HH_HpHm;
    gslpp::complex I_HH_Hp;
    
//    double LAMu;
    double LAMc;
    double LAMt;
//    double LAMd;
    double LAMs;
    double LAMb;
//    double LAMe;
    double LAMmu;
    double LAMtau;
    double LAMw;
    double LAMhp;
    
    gslpp::complex A_HH_F;
    gslpp::complex A_HH_U;
    gslpp::complex A_HH_D;
    gslpp::complex A_HH_L;
    gslpp::complex A_HH_W;
    gslpp::complex A_HH_Hp;
    
    double SigmaTotSM_H;
    double SigmaggF_H;  
    double SigmabbF_H;
    double SigmaVBF_H;
    double SigmattF_H;
    double SigmaVH_H;
    
    double SigmaSum;
     
    double BrSM_Htocc;
    double BrSM_Htobb;
    double BrSM_Htott;
    double BrSM_Htomumu;
    double BrSM_Htotautau;
    double BrSM_HtoWW;
    double BrSM_HtoZZ;

    double GammaHtotSM;

    double GammaHhh;
    double GammaHHpHm;
    double GammaHAA;
    double GammaHAZ;
    double GammaHHpW;

    double GammaHtot;

    double Br_Htott;
    double Br_Htobb;
    double Br_Htotautau;
    double Br_HtoWW;
    double Br_HtoZZ;
    double Br_Htogaga;
    double Br_Htohh;
    double Br_htobb;
    double Br_htogaga;
};



class  Hobs_ggF_H_tautau: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_H_tautau(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_bbF_H_tautau: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_bbF_H_tautau(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_H_gaga: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_H_gaga(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_H_ZZ: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_H_ZZ(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_ggF_H_WW: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_ggF_H_WW(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_H_WW: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_H_WW(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_VBF_H_WW: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_VBF_H_WW(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_H_hh: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_H_hh(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_H_hh_bbbb: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_H_hh_bbbb(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_H_tt: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_H_tt(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_H_bb: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_H_bb(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  Hobs_H_hh_gagabb: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    Hobs_H_hh_gagabb(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_ggF_H_tautau_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_ggF_H_tautau_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_bbF_H_tautau_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_bbF_H_tautau_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_H_gaga_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_H_gaga_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_H_ZZ_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_H_ZZ_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_ggF_H_WW_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_ggF_H_WW_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_H_WW_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_H_WW_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_VBF_H_WW_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_VBF_H_WW_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_H_hh_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_H_hh_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_H_hh_bbbb_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_H_hh_bbbb_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_H_tt_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_H_tt_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_H_bb_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_H_bb_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



class  log10_H_hh_gagabb_TH: public heavyHiggsCache {
public:
    
    /**
     * @brief Constructor.
     */
    log10_H_hh_gagabb_TH(const StandardModel& SM_i);
    
    /**
     * @return 
     */
    double computeThValue ();
    
private:
    
};



#endif	/* HEAVYHIGGSCACHE_H */
