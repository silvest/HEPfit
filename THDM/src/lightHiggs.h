/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LIGHTHIGGS_H
#define	LIGHTHIGGS_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"
//#include "THDMfunctions.h"
#include "THDMcache.h"

/**
 * @class lightHiggs
 * @ingroup THDM
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The squared relative @f$\gamma \gamma@f$ and @f$Z\gamma@f$ couplings are calculated at one-loop
 * following @cite Gunion:1989we.
 */
//class lightHiggs : public ThObservable {
//public:
//    lightHiggs(const StandardModel& SM_i);
//    virtual ~lightHiggs();
////    void computeSignalStrengthQuantities();
//    //void updateSignalStrengthParameters();
////    void updateSignalStrengthQuantities();
////    void updateSignalStrengthQuantities(
////            const double bma, const double tanb, const double m12_2,
////            const double mHp2, const double MW, const double cW2,
////            const double mHl, const double vev, const double Mt,
////            const double Mb, const double Mtau, const double Mc,
////            const double Ms, const double Mmu, const double Mu,
////            const double Md, const double Me, const double MZ);
//
//    /**
//     * @brief Empty function
//     */
//    double computeThValue();
//
////    mutable double SigStrcache[18];
//
//protected:
////    double bma;
////    double tanb;
////    double m12_2;
////    double mHp2;
////    double MW;
////    double cW2;
////    double mHl;
////    double vev;
////    double Mt;
////    double Mb;
////    double Mtau;
////    double Mc;
////    double Ms;
////    double Mmu;
////    double Mu;
////    double Md;
////    double Me;
////    double MZ;
////
////    /**
////     * @brief SM branching ratio of @f$h\to b \bar b@f$.
////     * @return @f$BR{\text SM}(h\to b \bar b)@f$
////     */
////    double BrSM_htobb;
////
////    /**
////     * @brief SM branching ratio of @f$h\to \gamma \gamma@f$.
////     * @return @f$BR{\text SM}(h\to \gamma \gamma)@f$
////     */
////    double BrSM_htogaga;
////
////    /**
////     * @brief SM branching ratio of @f$h\to \tau \tau@f$.
////     * @return @f$BR{\text SM}(h\to \tau \tau)@f$
////     */
////    double BrSM_htotautau;
////
////    /**
////     * @brief Squared relative coupling of @f$h@f$ to two down type quarks.
////     * @return @f$r^{(h)}_{Q_dQ_d}@f$
////     * @details Depends on the type of @f$Z_2@f$ symmetry.
////     */
////    double rh_QdQd;
////
////    /**
////     * @brief Squared relative coupling of @f$h@f$ to two massive vector bosons.
////     * @return @f$r^{(h)}_{WW}=r^{(h)}_{ZZ}@f$
////     */
////    double rh_VV;
////
////    /**
////     * @brief Squared relative coupling of @f$h@f$ to two charged leptons.
////     * @return @f$r^{(h)}_{\ell \ell}@f$
////     * @details Depends on the type of @f$Z_2@f$ symmetry.
////     */
////    double rh_ll;
////
////    /**
////     * @brief Squared relative coupling of @f$h@f$ to two photons.
////     * @return @f$r^{(h)}_{\gamma \gamma}@f$
////     * @details Depends on the type of @f$Z_2@f$ symmetry.
////     */
////    double rh_gaga;
////
////    /**
////     * @brief Squared relative coupling of @f$h@f$ to two gluons.
////     * @return @f$r^{(h)}_{gg}@f$
////     * @details Depends on the type of @f$Z_2@f$ symmetry.
////     */
////    double rh_gg;
////
////    /**
////     * @brief Ratio of THDM and SM cross sections for ggF and tth production of h.
////     * @return @f$\sigma^{\text THDM}_{\text ggF+tth}/\sigma^{\text SM}_{\text ggF+tth}@f$
////     */
////    double ggF_tth;
////
////    /**
////     * @brief Ratio of THDM and SM cross sections for VBF and Vh production of h.
////     * @return @f$\sigma^{\text THDM}_{\text VBF+Vh}/\sigma^{\text SM}_{\text VBF+Vh}@f$
////     */
////    double VBF_Vh;
////
////    /**
////     * @brief Sum of the modified branching ratios.
////     * @return @f$\sum _i r^{(h)}_{i} BR^{\text SM}(h\to i)@f$
////     */
////    double sumModBRs;
////
////    /**
////     * @brief Total h decay rate in the THDM.
////     * @return @f$\Gamma_h@f$
////     */
////    double Gamma_h;
////
////private:
////    const THDM& myTHDM;
////
////    /**
////     * @brief Squared relative coupling of @f$h@f$ to a @f$Z@f$ boson and a photon.
////     * @return @f$r^{(h)}_{Z\gamma}@f$
////     * @details Depends on the type of @f$Z_2@f$ symmetry.
////     */
////    double rh_Zga;
////
////    /**
////     * @brief Squared relative coupling of @f$h@f$ to two up type quarks.
////     * @return @f$r^{(h)}_{Q_uQ_u}@f$
////     */
////    double rh_QuQu;
//};

/**
 * @class THDM_BR_h_bb
 * @ingroup THDM
 * @brief THDM branching ratio of @f$h\to b \bar b@f$.
 */
class THDM_BR_h_bb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    THDM_BR_h_bb(const StandardModel& SM_i);
    
    /**
     * @return @return @f$BR{\text THDM}(h\to b \bar b)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class THDM_BR_h_gaga
 * @ingroup THDM
 * @brief THDM branching ratio of @f$h\to \gamma \gamma@f$.
 */
class THDM_BR_h_gaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    THDM_BR_h_gaga(const StandardModel& SM_i);
    
    /**
     * @return @f$BR{\text THDM}(h\to \gamma \gamma)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class THDM_BR_h_tautau
 * @ingroup THDM
 * @brief THDM branching ratio of @f$h\to \tau \tau@f$.
 */
class THDM_BR_h_tautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    THDM_BR_h_tautau(const StandardModel& SM_i);
    
    /**
     * @return @f$BR{\text THDM}(h\to \tau \tau)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htobb
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two b quarks.
 */
class ggF_tth_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to b\bar b)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htoWW
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$W@f$ bosons.
 */
class ggF_tth_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to WW)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htotautau
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$\tau@f$ leptons.
 */
class ggF_tth_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \tau\tau)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htoZZ
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$Z@f$ bosons.
 */
class ggF_tth_htoZZ : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to ZZ)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htogaga
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two photons.
 */
class ggF_tth_htogaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \gamma\gamma)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htobb
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two b quarks.
 */
class VBF_Vh_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to b\bar b)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htoWW
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$W@f$ bosons.
 */
class VBF_Vh_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to WW)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htotautau
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$\tau@f$ leptons.
 */
class VBF_Vh_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to \tau\tau)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htoZZ
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$Z@f$ bosons.
 */
class VBF_Vh_htoZZ : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to ZZ)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htogaga
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two photons.
 */
class VBF_Vh_htogaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to \gamma\gamma)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class Gamma_h_THDM
 * @ingroup THDM
 * @brief Total h decay rate in the THDM.
 */
class Gamma_h_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Gamma_h_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$\Gamma_h@f$ in units of GeV
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class rh_gaga_THDM
 * @ingroup THDM
 * @brief Squared relative coupling of @f$h@f$ to two photons.
 */
class rh_gaga_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    rh_gaga_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(h)}_{\gamma \gamma}@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class rh_gg_THDM
 * @ingroup THDM
 * @brief Squared relative coupling of @f$h@f$ to two gluons.
 */
class rh_gg_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    rh_gg_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(h)}_{gg}@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};


#endif	/* HIGGSSIGSTR_H */
