/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LIGHTHIGGS_H
#define	LIGHTHIGGS_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"
#include "THDMfunctions.h"


/**
 * @class lightHiggs
 * @brief .
 */
class lightHiggs : public ThObservable {
public:
    lightHiggs(const StandardModel& SM_i);
    virtual ~lightHiggs();
    void computeParameters();
    
    double computeThValue();
    double THDM_BR_h_bb();
    double THDM_BR_h_gaga();
    double THDM_BR_h_tautau();
    
    double cos_bpa;

protected:
    
    double ggF_tth;
    double VBF_Vh;
    double sum;
    double rh_QdQd;
    double rh_VV;
    double rh_ll;
    double rh_gaga;
    double Br_htobb;
    double Br_htogaga;
    double Br_htotautau;

private:
    const THDM * myTHDM;
    const StandardModel& mySM;
    
    double Br_htoWW;
    double Br_htoZZ;
    double Br_htogg;
    double Br_htoZga;
    double Br_htocc;
    
    double SigmaggF;
    double Sigmaggh_tt;
    double Sigmaggh_bb;
    double Sigmatth;
    double SigmaVBF;
    double SigmaWh;
    double SigmaZh;
    double SigmaVh;
    
    double pc_ggF;
    double pc_tth;
    double pc_VBF;
    double pc_Vh;
    
    double rh_gg;
    double rh_Zga;
    double rh_QuQu;
    
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
    
    gslpp::complex I_h_ferm;
    gslpp::complex fermU;
    gslpp::complex fermD;
    gslpp::complex fermL;
    gslpp::complex I_h_W;
    double ghHpHm;
    gslpp::complex I_h_Hp;

    double ABS1;
    double ABS2;

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
    
    gslpp::complex A_h_F;
    gslpp::complex A_h_U;
    gslpp::complex A_h_D;
    gslpp::complex A_h_L;
    gslpp::complex A_h_W;
    gslpp::complex A_h_Hp;

    double ABS3;
    double ABS4;
};



class ggF_tth_htobb : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htobb(const StandardModel& SM_i);
    
    /**
     * @return ggF_tth_htobb
     */
    double computeThValue ();
    
private:
    
};



class ggF_tth_htoWW : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoWW(const StandardModel& SM_i);
    
    /**
     * @return ggF_tth_htoWW
     */
    double computeThValue ();
    
private:
    
};



class ggF_tth_htotautau : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htotautau(const StandardModel& SM_i);
    
    /**
     * @return ggF_tth_htotautau
     */
    double computeThValue ();
    
private:
    
};



class ggF_tth_htoZZ : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return ggF_tth_htoZZ
     */
    double computeThValue ();
    
private:
    
};



class ggF_tth_htogaga : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htogaga(const StandardModel& SM_i);
    
    /**
     * @return ggF_tth_htogaga
     */
    double computeThValue ();
    
private:
    
};



class VBF_Vh_htobb : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htobb(const StandardModel& SM_i);
    
    /**
     * @return VBF_Vh_htobb
     */
    double computeThValue ();
    
private:
    
};



class VBF_Vh_htoWW : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htoWW(const StandardModel& SM_i);
    
    /**
     * @return VBF_Vh_htoWW
     */
    double computeThValue ();
    
private:
    
};



class VBF_Vh_htotautau : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htotautau(const StandardModel& SM_i);
    
    /**
     * @return VBF_Vh_htotautau
     */
    double computeThValue ();
    
private:
    
};



class VBF_Vh_htoZZ : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return VBF_Vh_htoZZ
     */
    double computeThValue ();
    
private:
    
};



class VBF_Vh_htogaga : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htogaga(const StandardModel& SM_i);
    
    /**
     * @return VBF_Vh_htogaga
     */
    double computeThValue ();
    
private:
    
};


#endif	/* HIGGSSIGSTR_H */