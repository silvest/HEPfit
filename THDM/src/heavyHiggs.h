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
#include "THDMfunctions.h"
#include "THDMcache.h"
#include "lightHiggs.h"

/**
 * @class heavyHiggs
 * @ingroup THDM 
 * @brief Base class for direct heavy Higgs search observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The @f$\gamma \gamma@f$, @f$Z\gamma@f$ and @f$gg@f$ decay widths are calculated at one-loop
 * following @cite Gunion:1989we and @cite Aoki:2009ha.
 */
class heavyHiggs : public ThObservable {
public:
    heavyHiggs(const StandardModel& SM_i);
    virtual ~heavyHiggs();
    void computeParameters();

    /**
     * @brief Empty function
     */
    double computeThValue();

protected:

    THDMcache * mycache;
    lightHiggs * mylightHiggs;

    /**
     * @brief The heavy CP-even Higgs mass. (Required for the experimental tables.)
     * @return @f$m_H@f$
     */
    double mHh;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)@f$
     */
    double ggF_H_tautau_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to H\to \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)@f$
     */
    double bbF_H_tautau_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to \gamma\gamma@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)@f$
     */
    double ggF_H_gaga_TH;  

    /**
     * @brief Signal strength for the process @f$pp\to H\to ZZ@f$ at the LHC with 8 TeV.
     * @return @f$\mu_H^{\text{THDM}}(H\to ZZ)=[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ)] / [\sigma^{\text{SM}}_{pp\to H}\cdot BR^{\text{SM}}(H\to ZZ)]@f$
     */
    double pp_H_ZZ_TH; 

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to WW@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)@f$
     */
    double ggF_H_WW_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$VV\to H\to WW@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)@f$
     */
    double VBF_H_WW_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to hh@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)@f$
     */
    double ggF_H_hh_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to H\to hh\to b\bar b \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)@f$
     */
    double ggF_H_hh_bbtautau_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh\to b\bar b b\bar b@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)@f$
     */
    double pp_H_hh_bbbb_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to hh\to \gamma\gamma b\bar b@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)@f$
     */
    double pp_H_hh_gagabb_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to H\to t\bar t@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)@f$
     */
    double pp_H_tt_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to H\to b\bar b@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)@f$
     */
    double bbF_H_bb_TH;

private:
    const THDM * myTHDM;
    const StandardModel& mySM;

    /**
     * @brief Total decay width of the heavy CP-even Higgs @f$H@f$.
     * @return @f$\Gamma^{\text tot}_H@f$
     */
    double GammaHtot;
};

/**
 * @class Hobs_ggF_H_tautau_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$.
 */
class Hobs_ggF_H_tautau_ATLAS: public heavyHiggs {
public:

    /**
     * @brief Hobs_ggF_H_tautau_ATLAS constructor.
     */
    Hobs_ggF_H_tautau_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_H_tautau_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$.
 */
class Hobs_ggF_H_tautau_CMS: public heavyHiggs {
public:

    /**
     * @brief Hobs_ggF_H_tautau_CMS constructor.
     */
    Hobs_ggF_H_tautau_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_bbF_H_tautau_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to H\to \tau\tau@f$.
 */
class Hobs_bbF_H_tautau_ATLAS: public heavyHiggs {
public:

    /**
     * @brief Hobs_bbF_H_tautau_ATLAS constructor.
     */
    Hobs_bbF_H_tautau_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_bbF_H_tautau_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to H\to \tau\tau@f$.
 */
class Hobs_bbF_H_tautau_CMS: public heavyHiggs {
public:

    /**
     * @brief Hobs_bbF_H_tautau_CMS constructor.
     */
    Hobs_bbF_H_tautau_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_H_gaga_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to \gamma\gamma@f$.
 */
class Hobs_ggF_H_gaga_ATLAS: public heavyHiggs {
public:

    /**
     * @brief Hobs_ggF_H_gaga_ATLAS constructor.
     */
    Hobs_ggF_H_gaga_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \gamma\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_H_gaga_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to H\to \gamma\gamma@f$.
 */
class Hobs_ggF_H_gaga_CMS: public heavyHiggs {
public:

    /**
     * @brief Hobs_ggF_H_gaga_CMS constructor.
     */
    Hobs_ggF_H_gaga_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \gamma\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_pp_H_ZZ_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the signal strength of the process @f$pp \to H\to ZZ@f$.
 */
class Hobs_pp_H_ZZ_CMS: public heavyHiggs {
public:

    /**
     * @brief Hobs_pp_H_ZZ_CMS constructor.
     */
    Hobs_pp_H_ZZ_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\mu_H^{\text{THDM}}(H\to ZZ)]_{\text{theo}} / [\mu_H(H\to ZZ)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_H_WW_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to WW@f$.
 */
class Hobs_ggF_H_WW_ATLAS: public heavyHiggs {
public:

    /**
     * @brief Hobs_ggF_H_WW_ATLAS constructor.
     */
    Hobs_ggF_H_WW_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_VBF_H_WW_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV \to H\to WW@f$.
 */
class Hobs_VBF_H_WW_ATLAS: public heavyHiggs {
public:

    /**
     * @brief Hobs_VBF_H_WW_ATLAS constructor.
     */
    Hobs_VBF_H_WW_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} / [\sigma_{VV\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_H_hh_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to hh@f$.
 */
class Hobs_ggF_H_hh_ATLAS: public heavyHiggs {
public:

    /**
     * @brief Hobs_ggF_H_hh_ATLAS constructor.
     */
    Hobs_ggF_H_hh_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to hh)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_H_hh_bbtautau_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to H\to hh\to b\bar b \tau\tau@f$.
 */
class Hobs_ggF_H_hh_bbtautau_CMS: public heavyHiggs {
public:

    /**
     * @brief Hobs_ggF_H_hh_bbtautau_CMS constructor.
     */
    Hobs_ggF_H_hh_bbtautau_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to hh\to b\bar b \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_pp_H_hh_bbbb_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to H\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_pp_H_hh_bbbb_CMS: public heavyHiggs {
public:

    /**
     * @brief Hobs_pp_H_hh_bbbb_CMS constructor.
     */
    Hobs_pp_H_hh_bbbb_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_pp_H_hh_gagabb_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to H\to hh\to \gamma\gamma b\bar b@f$.
 */
class Hobs_pp_H_hh_gagabb_CMS: public heavyHiggs {
public:

    /**
     * @brief Hobs_pp_H_hh_gagabb_CMS constructor.
     */
    Hobs_pp_H_hh_gagabb_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to \gamma\gamma b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_pp_H_tt_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp \to H\to t\bar t@f$.
 */
class Hobs_pp_H_tt_ATLAS: public heavyHiggs {
public:

    /**
     * @brief Hobs_pp_H_tt_ATLAS constructor.
     */
    Hobs_pp_H_tt_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to t\bar t))_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_bbF_H_bb_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b \to H\to b\bar b@f$.
 */
class Hobs_bbF_H_bb_CMS: public heavyHiggs {
public:

    /**
     * @brief Hobs_bbF_H_bb_CMS constructor.
     */
    Hobs_bbF_H_bb_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_H_tautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$.
 */
class log10_ggF_H_tautau_TH: public heavyHiggs {
public:

    /**
     * @brief log10_ggF_H_tautau_TH constructor.
     */
    log10_ggF_H_tautau_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_bbF_H_tautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to H\to \tau\tau@f$.
 */
class log10_bbF_H_tautau_TH: public heavyHiggs {
public:

    /**
     * @brief log10_bbF_H_tautau_TH constructor.
     */
    log10_bbF_H_tautau_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_H_gaga_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to \gamma\gamma@f$.
 */
class log10_ggF_H_gaga_TH: public heavyHiggs {
public:

    /**
     * @brief log10_ggF_H_gaga_TH constructor.
     */
    log10_ggF_H_gaga_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_pp_H_ZZ_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to ZZ@f$.
 */
class log10_pp_H_ZZ_TH: public heavyHiggs {
public:

    /**
     * @brief log10_pp_H_ZZ_TH constructor.
     */
    log10_pp_H_ZZ_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_H_WW_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to WW@f$.
 */
class log10_ggF_H_WW_TH: public heavyHiggs {
public:

    /**
     * @brief log10_ggF_H_WW_TH constructor.
     */
    log10_ggF_H_WW_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_VBF_H_WW_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to H\to WW@f$.
 */
class log10_VBF_H_WW_TH: public heavyHiggs {
public:

    /**
     * @brief log10_VBF_H_WW_TH constructor.
     */
    log10_VBF_H_WW_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_H_hh_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to hh@f$.
 */
class log10_ggF_H_hh_TH: public heavyHiggs {
public:

    /**
     * @brief log10_ggF_H_hh_TH constructor.
     */
    log10_ggF_H_hh_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_H_hh_bbtautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to hh\to b\bar b \tau\tau@f$.
 */
class log10_ggF_H_hh_bbtautau_TH: public heavyHiggs {
 public:

  /**                                                                                                                                         
   * @brief log10_ggF_H_hh_bbtautau_TH constructor.                                                                                                                      
   */
  log10_ggF_H_hh_bbtautau_TH(const StandardModel& SM_i);

  /**                                                                                                                                         
   * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)]@f$
   */
  double computeThValue();
};

/**
 * @class log10_pp_H_hh_bbbb_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b b\bar b@f$.
 */
class log10_pp_H_hh_bbbb_TH: public heavyHiggs {
public:

    /**
     * @brief log10_pp_H_hh_bbbb_TH constructor.
     */
    log10_pp_H_hh_bbbb_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_pp_H_hh_gagabb_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to \gamma\gamma b\bar b@f$.
 */
class log10_pp_H_hh_gagabb_TH: public heavyHiggs {
public:

    /**
     * @brief log10_pp_H_hh_gagabb_TH constructor.
     */
    log10_pp_H_hh_gagabb_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_pp_H_tt_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to t\bar t@f$.
 */
class log10_pp_H_tt_TH: public heavyHiggs {
public:

    /**
     * @brief log10_pp_H_tt_TH constructor.
     */
    log10_pp_H_tt_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_bbF_H_bb_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to H\to b\bar b@f$.
 */
class log10_bbF_H_bb_TH: public heavyHiggs {
public:

    /**
     * @brief log10_bbF_H_bb_TH constructor.
     */
    log10_bbF_H_bb_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]@f$
     */
    double computeThValue();
};

#endif	/* HEAVYHIGGS_H */
