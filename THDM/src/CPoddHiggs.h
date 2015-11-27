/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CPODDHIGGS_H
#define	CPODDHIGGS_H

#include <stdexcept>
#include "ThObservable.h"
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
 * @class CPoddHiggs
 * @ingroup THDM 
 * @brief Base class for direct CP-odd Higgs search observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The @f$\gamma \gamma@f$, @f$Z\gamma@f$ and @f$gg@f$ decay widths are calculated at one-loop
 * following @cite Gunion:1989we and @cite Aoki:2009ha.
 */
class CPoddHiggs : public ThObservable {
public:
    CPoddHiggs(const StandardModel& SM_i);
    virtual ~CPoddHiggs();
    void computeParameters();

    /**
     * @brief Empty function
     */
    double computeThValue();
    
protected:
    
    THDMcache * mycache;
    lightHiggs * mylightHiggs;

    /**
     * @brief The CP-odd Higgs mass. (Required for the experimental tables.)
     * @return @f$m_A@f$
     */
    double mA;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)@f$
     */
    double ggF_A_tautau_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to A\to \tau\tau@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)@f$
     */
    double bbF_A_tautau_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to \gamma\gamma@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)@f$
     */
    double ggF_A_gaga_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ \to b\bar b \ell \ell@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b \ell \ell)@f$
     */
    double ggF_A_hZ_bbll_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ \to b\bar b Z@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b Z)@f$
     */
    double ggF_A_hZ_bbZ_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ \to \tau \tau \ell \ell@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau \ell \ell)@f$
     */
    double ggF_A_hZ_tautaull_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$gg\to A\to hZ \to \tau \tau Z@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau Z)@f$
     */
    double ggF_A_hZ_tautauZ_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$pp\to A\to t\bar t@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to t\bar t)@f$
     */
    double pp_A_tt_TH;

    /**
     * @brief Cross section times branching ratio for the process @f$b\bar b\to A\to b\bar b@f$ at the LHC with 8 TeV.
     * @return @f$\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to b\bar b)@f$
     */
    double bbF_A_bb_TH;

private:
    const THDM * myTHDM;
    const StandardModel& mySM;

    /**
     * @brief Total decay width of the CP-odd Higgs @f$A@f$.
     * @return @f$\Gamma^{\text tot}_A@f$
     */
    double GammaAtot;
};

/**
 * @class Hobs_ggF_A_tautau_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to A\to \tau\tau@f$.
 */
class Hobs_ggF_A_tautau_ATLAS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_ggF_A_tautau_ATLAS constructor.
     */
    Hobs_ggF_A_tautau_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_A_tautau_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to A\to \tau\tau@f$.
 */
class Hobs_ggF_A_tautau_CMS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_ggF_A_tautau_CMS constructor.
     */
    Hobs_ggF_A_tautau_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_bbF_A_tautau_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to A\to \tau\tau@f$.
 */
class Hobs_bbF_A_tautau_ATLAS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_bbF_A_tautau_ATLAS constructor.
     */
    Hobs_bbF_A_tautau_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to A}\cdot BR(A\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_bbF_A_tautau_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to A\to \tau\tau@f$.
 */
class Hobs_bbF_A_tautau_CMS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_bbF_A_tautau_CMS constructor.
     */
    Hobs_bbF_A_tautau_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to A}\cdot BR(A\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_A_gaga_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to A\to \gamma\gamma@f$.
 */
class Hobs_ggF_A_gaga_ATLAS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_ggF_A_gaga_ATLAS constructor.
     */
    Hobs_ggF_A_gaga_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to \gamma\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_A_gaga_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to A\to \gamma\gamma@f$.
 */
class Hobs_ggF_A_gaga_CMS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_ggF_A_gaga_CMS constructor.
     */
    Hobs_ggF_A_gaga_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to \gamma\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_A_hZ_bbll_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to A\to hZ \to b\bar b \ell \ell@f$.
 */
class Hobs_ggF_A_hZ_bbll_CMS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_ggF_A_hZ_bbll_CMS constructor.
     */
    Hobs_ggF_A_hZ_bbll_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b \ell \ell)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to hZ \to b\bar b \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_A_hZ_bbZ_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to A\to hZ \to b\bar b Z@f$.
 */
class Hobs_ggF_A_hZ_bbZ_ATLAS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_ggF_A_hZ_bbZ_ATLAS constructor.
     */
    Hobs_ggF_A_hZ_bbZ_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b Z)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to hZ \to b\bar b Z)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_A_hZ_tautaull_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to A\to hZ \to \tau \tau \ell \ell@f$.
 */
class Hobs_ggF_A_hZ_tautaull_CMS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_ggF_A_hZ_tautaull_CMS constructor.
     */
    Hobs_ggF_A_hZ_tautaull_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau \ell \ell)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to hZ \to \tau \tau \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_ggF_A_hZ_tautauZ_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to A\to hZ \to \tau \tau Z@f$.
 */
class Hobs_ggF_A_hZ_tautauZ_ATLAS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_ggF_A_hZ_tautauZ_ATLAS constructor.
     */
    Hobs_ggF_A_hZ_tautauZ_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau b Z)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to hZ \to \tau \tau Z)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_pp_A_tt_ATLAS
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp \to A\to t\bar t@f$.
 */
class Hobs_pp_A_tt_ATLAS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_pp_A_tt_ATLAS constructor.
     */
    Hobs_pp_A_tt_ATLAS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to t\bar t)]_{\text{theo}} / [\sigma_{pp\to A}\cdot BR(A\to t\bar t))_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class Hobs_bbF_A_bb_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b \to A\to b\bar b@f$.
 */
class Hobs_bbF_A_bb_CMS: public CPoddHiggs {
public:

    /**
     * @brief Hobs_bbF_A_bb_CMS constructor.
     */
    Hobs_bbF_A_bb_CMS(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to A}\cdot BR(A\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_A_tautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to \tau\tau@f$.
 */
class log10_ggF_A_tautau_TH: public CPoddHiggs {
public:

    /**
     * @brief log10_ggF_A_tautau_TH constructor.
     */
    log10_ggF_A_tautau_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_bbF_A_tautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to A\to \tau\tau@f$.
 */
class log10_bbF_A_tautau_TH: public CPoddHiggs {
public:

    /**
     * @brief log10_bbF_A_tautau_TH constructor.
     */
    log10_bbF_A_tautau_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_A_gaga_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to \gamma\gamma@f$.
 */
class log10_ggF_A_gaga_TH: public CPoddHiggs {
public:

    /**
     * @brief log10_ggF_A_gaga_TH constructor.
     */
    log10_ggF_A_gaga_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_A_hZ_bbll_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to hZ \to b\bar b \ell \ell@f$.
 */
class log10_ggF_A_hZ_bbll_TH: public CPoddHiggs {
public:

    /**
     * @brief log10_ggF_A_hZ_bbll_TH constructor.
     */
    log10_ggF_A_hZ_bbll_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ\to b\bar b \ell \ell)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_A_hZ_bbZ_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to hZ \to b\bar b Z@f$.
 */
class log10_ggF_A_hZ_bbZ_TH: public CPoddHiggs {
public:

    /**
     * @brief log10_ggF_A_hZ_bbZ_TH constructor.
     */
    log10_ggF_A_hZ_bbZ_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ\to b\bar b Z)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_A_hZ_tautaull_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to hZ \to \tau \tau \ell \ell@f$.
 */
class log10_ggF_A_hZ_tautaull_TH: public CPoddHiggs {
public:

    /**
     * @brief log10_ggF_A_hZ_tautaull_TH constructor.
     */
    log10_ggF_A_hZ_tautaull_TH(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ\to \tau \tau \ell \ell)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_ggF_A_hZ_tautauZ_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to hZ \to \tau \tau Z@f$.
 */
class log10_ggF_A_hZ_tautauZ_TH: public CPoddHiggs {
public:

    /**
     * @brief log10_ggF_A_hZ_tautauZ_TH constructor.
     */
    log10_ggF_A_hZ_tautauZ_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ\to \tau \tau Z)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_pp_A_tt_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to A\to t\bar t@f$.
 */
class log10_pp_A_tt_TH: public CPoddHiggs {
public:

    /**
     * @brief log10_pp_A_tt_TH constructor.
     */
    log10_pp_A_tt_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to t\bar t)]@f$
     */
    double computeThValue();
};

/**
 * @class log10_bbF_A_bb_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to A\to b\bar b@f$.
 */
class log10_bbF_A_bb_TH: public CPoddHiggs {
public:

    /**
     * @brief log10_bbF_A_bb_TH constructor.
     */
    log10_bbF_A_bb_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to b\bar b)]@f$
     */
    double computeThValue();
};

/**
 * @}
 */

#endif	/* CPODDHIGGS_H */
