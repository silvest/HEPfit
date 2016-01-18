/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEAVYHIGGS_H
#define	HEAVYHIGGS_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"
#include "THDMcache.h"

/**
 * @class heavyHiggs
 * @ingroup THDM 
 * @brief Base class for direct heavy Higgs search observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The @f$\gamma \gamma@f$, @f$Z\gamma@f$ and @f$gg@f$ decay widths are calculated at one-loop
 * following @cite Gunion:1989we and @cite Aoki:2009ha.
 */

/**
 * @class Hobs_ggF_H_tautau_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$.
 */
class Hobs_ggF_H_tautau_ATLAS: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_tautau_ATLAS constructor.
     */
    Hobs_ggF_H_tautau_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_tautau_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$.
 */
class Hobs_ggF_H_tautau_CMS: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_tautau_CMS constructor.
     */
    Hobs_ggF_H_tautau_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_H_tautau_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to H\to \tau\tau@f$.
 */
class Hobs_bbF_H_tautau_ATLAS: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_H_tautau_ATLAS constructor.
     */
    Hobs_bbF_H_tautau_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_H_tautau_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to H\to \tau\tau@f$.
 */
class Hobs_bbF_H_tautau_CMS: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_H_tautau_CMS constructor.
     */
    Hobs_bbF_H_tautau_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_gaga_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to \gamma\gamma@f$.
 */
class Hobs_ggF_H_gaga_ATLAS: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_gaga_ATLAS constructor.
     */
    Hobs_ggF_H_gaga_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \gamma\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_gaga_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to H\to \gamma\gamma@f$.
 */
class Hobs_ggF_H_gaga_CMS: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_gaga_CMS constructor.
     */
    Hobs_ggF_H_gaga_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \gamma\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_ZZ_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the signal strength of the process @f$pp \to H\to ZZ@f$.
 */
class Hobs_pp_H_ZZ_CMS: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_ZZ_CMS constructor.
     */
    Hobs_pp_H_ZZ_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\mu_H^{\text{THDM}}(H\to ZZ)]_{\text{theo}} / [\mu_H(H\to ZZ)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_WW_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to WW@f$.
 */
class Hobs_ggF_H_WW_ATLAS: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_WW_ATLAS constructor.
     */
    Hobs_ggF_H_WW_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_VBF_H_WW_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV \to H\to WW@f$.
 */
class Hobs_VBF_H_WW_ATLAS: public ThObservable {
public:

    /**
     * @brief Hobs_VBF_H_WW_ATLAS constructor.
     */
    Hobs_VBF_H_WW_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} / [\sigma_{VV\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_hh_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to hh@f$.
 */
class Hobs_ggF_H_hh_ATLAS: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_hh_ATLAS constructor.
     */
    Hobs_ggF_H_hh_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to hh)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_hh_bbtautau_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to H\to hh\to b\bar b \tau\tau@f$.
 */
class Hobs_ggF_H_hh_bbtautau_CMS: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_hh_bbtautau_CMS constructor.
     */
    Hobs_ggF_H_hh_bbtautau_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to hh\to b\bar b \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_bbbb_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to H\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_pp_H_hh_bbbb_CMS: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_bbbb_CMS constructor.
     */
    Hobs_pp_H_hh_bbbb_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_gagabb_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to H\to hh\to \gamma\gamma b\bar b@f$.
 */
class Hobs_pp_H_hh_gagabb_CMS: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_gagabb_CMS constructor.
     */
    Hobs_pp_H_hh_gagabb_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to \gamma\gamma b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_tt_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp \to H\to t\bar t@f$.
 */
class Hobs_pp_H_tt_ATLAS: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_tt_ATLAS constructor.
     */
    Hobs_pp_H_tt_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to t\bar t))_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_H_bb_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b \to H\to b\bar b@f$.
 */
class Hobs_bbF_H_bb_CMS: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_H_bb_CMS constructor.
     */
    Hobs_bbF_H_bb_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_tautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$.
 */
class log10_ggF_H_tautau_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_tautau_TH constructor.
     */
    log10_ggF_H_tautau_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_bbF_H_tautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to H\to \tau\tau@f$.
 */
class log10_bbF_H_tautau_TH: public ThObservable {
public:

    /**
     * @brief log10_bbF_H_tautau_TH constructor.
     */
    log10_bbF_H_tautau_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_gaga_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to \gamma\gamma@f$.
 */
class log10_ggF_H_gaga_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_gaga_TH constructor.
     */
    log10_ggF_H_gaga_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_ZZ_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to ZZ@f$.
 */
class log10_pp_H_ZZ_TH: public ThObservable {
public:

    /**
     * @brief log10_pp_H_ZZ_TH constructor.
     */
    log10_pp_H_ZZ_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_WW_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to WW@f$.
 */
class log10_ggF_H_WW_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_WW_TH constructor.
     */
    log10_ggF_H_WW_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_VBF_H_WW_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to H\to WW@f$.
 */
class log10_VBF_H_WW_TH: public ThObservable {
public:

    /**
     * @brief log10_VBF_H_WW_TH constructor.
     */
    log10_VBF_H_WW_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_hh_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to hh@f$.
 */
class log10_ggF_H_hh_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_hh_TH constructor.
     */
    log10_ggF_H_hh_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_hh_bbtautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to hh\to b\bar b \tau\tau@f$.
 */
class log10_ggF_H_hh_bbtautau_TH: public ThObservable {
 public:

  /**                                                                                                                                         
   * @brief log10_ggF_H_hh_bbtautau_TH constructor.                                                                                                                      
   */
  log10_ggF_H_hh_bbtautau_TH(const StandardModel& SM_i);

  /**                                                                                                                                         
   * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)]@f$
   */
  double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_bbbb_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b b\bar b@f$.
 */
class log10_pp_H_hh_bbbb_TH: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_bbbb_TH constructor.
     */
    log10_pp_H_hh_bbbb_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_gagabb_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to \gamma\gamma b\bar b@f$.
 */
class log10_pp_H_hh_gagabb_TH: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_gagabb_TH constructor.
     */
    log10_pp_H_hh_gagabb_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_tt_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to t\bar t@f$.
 */
class log10_pp_H_tt_TH: public ThObservable {
public:

    /**
     * @brief log10_pp_H_tt_TH constructor.
     */
    log10_pp_H_tt_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_bbF_H_bb_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to H\to b\bar b@f$.
 */
class log10_bbF_H_bb_TH: public ThObservable {
public:

    /**
     * @brief log10_bbF_H_bb_TH constructor.
     */
    log10_bbF_H_bb_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Gamma_HH_THDM
 * @ingroup THDM
 * @brief Total H decay rate in the %THDM.
 */
class Gamma_HH_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Gamma_HH_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$\Gamma_H@f$ in units of GeV
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

///**
// * @class rHH_gaga_THDM
// * @ingroup THDM
// * @brief Squared relative coupling of @f$H@f$ to two photons.
// */
//class rHH_gaga_THDM : public ThObservable {
//public:
//    
//    /**
//     * @brief Constructor.
//     */
//    rHH_gaga_THDM(const StandardModel& SM_i);
//    
//    /**
//     * @return @f$r^{(H)}_{\gamma \gamma}@f$
//     */
//    double computeThValue ();
//};

/**
 * @class rHH_gg_THDM
 * @ingroup THDM
 * @brief Squared relative coupling of @f$H@f$ to two gluons.
 */
class rHH_gg_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor for the squared relative coupling of @f$H@f$ to two gluons.
     */
    rHH_gg_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(H)}_{gg}@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class BR_HH_hh_THDM
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$H@f$ to two @f$h@f$.
 */
class BR_HH_hh_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_HH_hh_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H\to hh)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class BR_HH_AA_THDM
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$H@f$ to two @f$A@f$.
 */
class BR_HH_AA_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_HH_AA_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H\to AA)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class BR_HH_HpHm_THDM
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$H@f$ to charged Higgs bosons.
 */
class BR_HH_HpHm_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_HH_HpHm_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H\to H^\pm H^\mp)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class BR_HH_AZ_THDM
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$H@f$ to an @f$A@f$ and a @f$Z@f$ boson.
 */
class BR_HH_AZ_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_HH_AZ_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H\to AZ)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class BR_HH_HpW_THDM
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$H@f$ to a charged Higgs boson and a @f$W@f$ boson.
 */
class BR_HH_HpW_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_HH_HpW_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H\to H^\pm W^\mp)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

#endif	/* HEAVYHIGGS_H */
