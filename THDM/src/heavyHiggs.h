/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
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
 * @class Hobs_ggF_H_tautau_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$.
 */
class Hobs_ggF_H_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_tautau_ATLAS8 constructor.
     */
    Hobs_ggF_H_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_tautau_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg\to H\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_tautau_ATLAS8 constructor.
     */
    Robs_ggF_H_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{gg\to H}\cdot BR(H\to \tau\tau)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_tautau_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$.
 */
class Hobs_ggF_H_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_tautau_CMS8 constructor.
     */
    Hobs_ggF_H_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_tautau_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$gg\to H\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_tautau_CMS8 constructor.
     */
    Robs_ggF_H_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{CMS,95\% observed}} \right) / [\sigma_{gg\to H}\cdot BR(H\to \tau\tau)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_H_tautau_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to H\to \tau\tau@f$.
 */
class Hobs_bbF_H_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_H_tautau_ATLAS8 constructor.
     */
    Hobs_bbF_H_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_bbF_H_tautau_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$b\bar b\to H\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_bbF_H_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_bbF_H_tautau_ATLAS8 constructor.
     */
    Robs_bbF_H_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} - [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_H_tautau_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to H\to \tau\tau@f$.
 */
class Hobs_bbF_H_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_H_tautau_CMS8 constructor.
     */
    Hobs_bbF_H_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_bbF_H_tautau_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$b\bar b\to H\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_bbF_H_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_bbF_H_tautau_CMS8 constructor.
     */
    Robs_bbF_H_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]_{\text{theo}} - [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{CMS,95\% observed}} \right) / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau\tau)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_gaga_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp \to H\to \gamma\gamma@f$.
 */
class Hobs_pp_H_gaga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_gaga_ATLAS8 constructor.
     */
    Hobs_pp_H_gaga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to \gamma\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_gaga_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$pp \to H\to \gamma\gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_gaga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_gaga_ATLAS8 constructor.
     */
    Robs_pp_H_gaga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_gaga_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to H\to \gamma\gamma@f$.
 */
class Hobs_ggF_H_gaga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_gaga_CMS8 constructor.
     */
    Hobs_ggF_H_gaga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \gamma\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_gaga_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$gg \to H\to \gamma\gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_gaga_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_gaga_CMS8 constructor.
     */
    Robs_ggF_H_gaga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_mu_pp_H_VV_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the signal strength of the process @f$pp \to H\to VV@f$.
 */
class Hobs_mu_pp_H_VV_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_mu_pp_H_VV_CMS8 constructor.
     */
    Hobs_mu_pp_H_VV_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\mu_H^{\text{THDM}}(H\to VV)]_{\text{theo}} / [\mu_H(H\to VV)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_mu_pp_H_VV_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp \to H\to VV@f$ assuming a Gaussian likelihood.
 */
class Robs_mu_pp_H_VV_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_mu_pp_H_VV_CMS8 constructor.
     */
    Robs_mu_pp_H_VV_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\mu_H^{\text{THDM}}(H\to VV)]_{\text{theo}} - [\mu_H^{\text{THDM}}(H\to VV)]_{\text{CMS,95\% observed}} \right) / [\mu_H^{\text{THDM}}(H\to VV)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_ZZ_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to ZZ@f$.
 */
class Hobs_ggF_H_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_ZZ_ATLAS8 constructor.
     */
    Hobs_ggF_H_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_ZZ_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg \to H\to ZZ@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_ZZ_ATLAS8 constructor.
     */
    Robs_ggF_H_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_VBF_H_ZZ_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV \to H\to ZZ@f$.
 */
class Hobs_VBF_H_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VBF_H_ZZ_ATLAS8 constructor.
     */
    Hobs_VBF_H_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} / [\sigma_{VV\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_VBF_H_ZZ_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$VV \to H\to ZZ@f$ assuming a Gaussian likelihood.
 */
class Robs_VBF_H_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_VBF_H_ZZ_ATLAS8 constructor.
     */
    Robs_VBF_H_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} - [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_WW_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to WW@f$.
 */
class Hobs_ggF_H_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_WW_ATLAS8 constructor.
     */
    Hobs_ggF_H_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_WW_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg \to H\to WW@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_WW_ATLAS8 constructor.
     */
    Robs_ggF_H_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_VBF_H_WW_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV \to H\to WW@f$.
 */
class Hobs_VBF_H_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VBF_H_WW_ATLAS8 constructor.
     */
    Hobs_VBF_H_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} / [\sigma_{VV\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_VBF_H_WW_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$VV \to H\to WW@f$ assuming a Gaussian likelihood.
 */
class Robs_VBF_H_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_VBF_H_WW_ATLAS8 constructor.
     */
    Robs_VBF_H_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} - [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_hh_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to hh@f$.
 */
class Hobs_ggF_H_hh_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_hh_ATLAS8 constructor.
     */
    Hobs_ggF_H_hh_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to hh)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_hh_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg \to H\to hh@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_hh_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_hh_ATLAS8 constructor.
     */
    Robs_ggF_H_hh_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to H\to hh@f$.
 */
class Hobs_pp_H_hh_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_CMS8 constructor.
     */
    Hobs_pp_H_hh_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_hh_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$pp \to H\to hh@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_hh_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_hh_CMS8 constructor.
     */
    Robs_pp_H_hh_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_hh_bbtautau_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to H\to hh\to b\bar b \tau\tau@f$.
 */
class Hobs_ggF_H_hh_bbtautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_hh_bbtautau_CMS8 constructor.
     */
    Hobs_ggF_H_hh_bbtautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to hh\to b\bar b \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_hh_bbtautau_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$gg \to H\to hh\to b\bar b \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_hh_bbtautau_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_hh_bbtautau_CMS8 constructor.
     */
    Robs_ggF_H_hh_bbtautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_bbbb_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to H\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_pp_H_hh_bbbb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_bbbb_CMS8 constructor.
     */
    Hobs_pp_H_hh_bbbb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_hh_bbbb_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp \to H\to hh\to b\bar b b\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_hh_bbbb_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_hh_bbbb_CMS8 constructor.
     */
    Robs_pp_H_hh_bbbb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_gagabb_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to H\to hh\to \gamma\gamma b\bar b@f$.
 */
class Hobs_pp_H_hh_gagabb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_gagabb_CMS8 constructor.
     */
    Hobs_pp_H_hh_gagabb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to \gamma\gamma b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_hh_gagabb_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp \to H\to hh\to \gamma\gamma b\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_hh_gagabb_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_hh_gagabb_CMS8 constructor.
     */
    Robs_pp_H_hh_gagabb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_tt_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to H\to t\bar t@f$.
 */
class Hobs_ggF_H_tt_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_tt_ATLAS8 constructor.
     */
    Hobs_ggF_H_tt_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to t\bar t))_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_tt_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg \to H\to t\bar t@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_tt_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_tt_ATLAS8 constructor.
     */
    Robs_ggF_H_tt_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_H_bb_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b \to H\to b\bar b@f$.
 */
class Hobs_bbF_H_bb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_H_bb_CMS8 constructor.
     */
    Hobs_bbF_H_bb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_bbF_H_bb_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$b\bar b \to H\to b\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_bbF_H_bb_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_bbF_H_bb_CMS8 constructor.
     */
    Robs_bbF_H_bb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ttF_H_tt_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$t\bar t \to H\to t\bar t@f$.
 */
class Hobs_ttF_H_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_ttF_H_tt_ATLAS13 constructor.
     */
    Hobs_ttF_H_tt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{t\bar t\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{theo}} / [\sigma_{t\bar t\to H}\cdot BR(H\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ttF_H_tt_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$t\bar t \to H\to t\bar t@f$ assuming a Gaussian likelihood.
 */
class Robs_ttF_H_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_ttF_H_tt_ATLAS13 constructor.
     */
    Robs_ttF_H_tt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{t\bar t\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{theo}} - [\sigma^{\text{THDM}}_{t\bar t\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{t\bar t\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_H_tt_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b \to H\to t\bar t@f$.
 */
class Hobs_bbF_H_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_H_tt_ATLAS13 constructor.
     */
    Hobs_bbF_H_tt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_bbF_H_tt_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$b\bar b \to H\to t\bar t@f$ assuming a Gaussian likelihood.
 */
class Robs_bbF_H_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_bbF_H_tt_ATLAS13 constructor.
     */
    Robs_bbF_H_tt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{theo}} - [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_tautau_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \tau \tau@f$.
 */
class Hobs_ggF_H_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_tautau_ATLAS13 constructor.
     */
    Hobs_ggF_H_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_tautau_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$gg\to H\to \tau \tau@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_tautau_ATLAS13 constructor.
     */
    Robs_ggF_H_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_H_tautau_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to H\to \tau \tau@f$.
 */
class Hobs_bbF_H_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_H_tautau_ATLAS13 constructor.
     */
    Hobs_bbF_H_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_bbF_H_tautau_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$b\bar b\to H\to \tau \tau@f$ assuming a Gaussian likelihood.
 */
class Robs_bbF_H_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_bbF_H_tautau_ATLAS13 constructor.
     */
    Robs_bbF_H_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{theo}} - [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_tautau_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \tau \tau@f$.
 */
class Hobs_ggF_H_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_tautau_CMS13 constructor.
     */
    Hobs_ggF_H_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_tautau_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$gg\to H\to \tau \tau@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_tautau_CMS13 constructor.
     */
    Robs_ggF_H_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_H_tautau_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to H\to \tau \tau@f$.
 */
class Hobs_bbF_H_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_H_tautau_CMS13 constructor.
     */
    Hobs_bbF_H_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{theo}} / [\sigma_{b\bar b\to H}\cdot BR(H\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_bbF_H_tautau_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$b\bar b\to H\to \tau \tau@f$ assuming a Gaussian likelihood.
 */
class Robs_bbF_H_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_bbF_H_tautau_CMS13 constructor.
     */
    Robs_bbF_H_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{theo}} - [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_gaga_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H\to \gamma \gamma@f$.
 */
class Hobs_pp_H_gaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_gaga_ATLAS13 constructor.
     */
    Hobs_pp_H_gaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_gaga_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$pp\to H\to \gamma \gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_gaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_gaga_ATLAS13 constructor.
     */
    Robs_pp_H_gaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_gaga_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \gamma \gamma@f$.
 */
class Hobs_ggF_H_gaga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_gaga_CMS13 constructor.
     */
    Hobs_ggF_H_gaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_gaga_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$gg\to H\to \gamma \gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_gaga_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_gaga_CMS13 constructor.
     */
    Robs_ggF_H_gaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_Zga_llga_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H\to Z\gamma@f$.
 */
class Hobs_pp_H_Zga_llga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_Zga_llga_ATLAS13 constructor.
     */
    Hobs_pp_H_Zga_llga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to Z\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_Zga_llga_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$pp\to H\to Z\gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_Zga_llga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_Zga_llga_ATLAS13 constructor.
     */
    Robs_pp_H_Zga_llga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_Zga_llga_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H\to Z\gamma@f$.
 */
class Hobs_pp_H_Zga_llga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_Zga_llga_CMS13 constructor.
     */
    Hobs_pp_H_Zga_llga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to Z\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_Zga_llga_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H\to Z\gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_Zga_llga_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_Zga_llga_CMS13 constructor.
     */
    Robs_pp_H_Zga_llga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_Zga_qqga_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H\to Z\gamma@f$.
 */
class Hobs_pp_H_Zga_qqga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_Zga_qqga_CMS13 constructor.
     */
    Hobs_pp_H_Zga_qqga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to Z\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_Zga_qqga_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H\to Z\gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_Zga_qqga_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_Zga_qqga_CMS13 constructor.
     */
    Robs_pp_H_Zga_qqga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_ZZ_llnunu_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to ZZ@f$.
 */
class Hobs_ggF_H_ZZ_llnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_ZZ_llnunu_ATLAS13 constructor.
     */
    Hobs_ggF_H_ZZ_llnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_ZZ_llnunu_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$gg\to H\to ZZ@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_ZZ_llnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_ZZ_llnunu_ATLAS13 constructor.
     */
    Robs_ggF_H_ZZ_llnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_ZZ_llll_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to ZZ \to \ell \ell \ell \ell@f$.
 */
class Hobs_ggF_H_ZZ_llll_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_ZZ_llll_ATLAS13 constructor.
     */
    Hobs_ggF_H_ZZ_llll_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_ZZ_llll_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$gg\to H\to ZZ \to \ell \ell \ell \ell@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_ZZ_llll_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_ZZ_llll_ATLAS13 constructor.
     */
    Robs_ggF_H_ZZ_llll_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_VBF_H_ZZ_llll_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV\to H\to ZZ \to \ell \ell \ell \ell@f$.
 */
class Hobs_VBF_H_ZZ_llll_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VBF_H_ZZ_llll_ATLAS13 constructor.
     */
    Hobs_VBF_H_ZZ_llll_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{theo}} / [\sigma_{VV\to H}\cdot BR(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_VBF_H_ZZ_llll_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$VV\to H\to ZZ \to \ell \ell \ell \ell@f$ assuming a Gaussian likelihood.
 */
class Robs_VBF_H_ZZ_llll_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_VBF_H_ZZ_llll_ATLAS13 constructor.
     */
    Robs_VBF_H_ZZ_llll_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{theo}} - [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_ZZ_llll_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H\to ZZ \to \ell \ell \ell \ell@f$.
 */
class Hobs_pp_H_ZZ_llll_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_ZZ_llll_CMS13 constructor.
     */
    Hobs_pp_H_ZZ_llll_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to ZZ \to \ell \ell \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_ZZ_llll_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H\to ZZ \to \ell \ell \ell \ell@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_ZZ_llll_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_ZZ_llll_CMS13 constructor.
     */
    Robs_pp_H_ZZ_llll_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_VBF_VH_H_ZZ_llll_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$(VBF+VH)\to H\to ZZ \to \ell \ell \ell \ell@f$.
 */
class Hobs_VBF_VH_H_ZZ_llll_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_VBF_VH_H_ZZ_llll_CMS13 constructor.
     */
    Hobs_VBF_VH_H_ZZ_llll_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{(VBF+VH)\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{theo}} / [\sigma_{(VBF+VH)\to H}\cdot BR(H\to ZZ \to \ell \ell \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_VBF_VH_H_ZZ_llll_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$(VBF+VH)\to H\to ZZ \to \ell \ell \ell \ell@f$ assuming a Gaussian likelihood.
 */
class Robs_VBF_VH_H_ZZ_llll_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_VBF_VH_H_ZZ_llll_CMS13 constructor.
     */
    Robs_VBF_VH_H_ZZ_llll_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{(VBF+VH)\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{theo}} - [\sigma^{\text{THDM}}_{(VBF+VH)\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{(VBF+VH)\to H}\cdot BR^{\text{THDM}}(H\to ZZ \to \ell \ell \ell \ell)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_ZZ_llqq_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to ZZ@f$.
 */
class Hobs_ggF_H_ZZ_llqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_ZZ_llqq_ATLAS13 constructor.
     */
    Hobs_ggF_H_ZZ_llqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_ZZ_llqq_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$gg\to H\to ZZ@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_ZZ_llqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_ZZ_llqq_ATLAS13 constructor.
     */
    Robs_ggF_H_ZZ_llqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_VBF_H_ZZ_llqq_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV\to H\to ZZ@f$.
 */
class Hobs_VBF_H_ZZ_llqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VBF_H_ZZ_llqq_ATLAS13 constructor.
     */
    Hobs_VBF_H_ZZ_llqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} / [\sigma_{VV\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_VBF_H_ZZ_llqq_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$VV\to H\to ZZ@f$ assuming a Gaussian likelihood.
 */
class Robs_VBF_H_ZZ_llqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_VBF_H_ZZ_llqq_ATLAS13 constructor.
     */
    Robs_VBF_H_ZZ_llqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} - [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_ZZ_nunuqq_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to ZZ@f$.
 */
class Hobs_ggF_H_ZZ_nunuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_ZZ_nunuqq_ATLAS13 constructor.
     */
    Hobs_ggF_H_ZZ_nunuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_ZZ_nunuqq_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$gg\to H\to ZZ@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_ZZ_nunuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_ZZ_nunuqq_ATLAS13 constructor.
     */
    Robs_ggF_H_ZZ_nunuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_ZZ_llqq_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H\to ZZ@f$.
 */
class Hobs_pp_H_ZZ_llqq_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_ZZ_llqq_CMS13 constructor.
     */
    Hobs_pp_H_ZZ_llqq_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to ZZ)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_ZZ_llqq_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H\to ZZ@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_ZZ_llqq_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_ZZ_llqq_CMS13 constructor.
     */
    Robs_pp_H_ZZ_llqq_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_WW_lnuqq_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to WW@f$.
 */
class Hobs_ggF_H_WW_lnuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_WW_lnuqq_ATLAS13 constructor.
     */
    Hobs_ggF_H_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_WW_lnuqq_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$gg\to H\to WW@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_WW_lnuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_WW_lnuqq_ATLAS13 constructor.
     */
    Robs_ggF_H_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_WW_enumunu_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to WW@f$.
 */
class Hobs_ggF_H_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_WW_enumunu_ATLAS13 constructor.
     */
    Hobs_ggF_H_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} / [\sigma_{gg\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_WW_enumunu_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$gg\to H\to WW@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_WW_enumunu_ATLAS13 constructor.
     */
    Robs_ggF_H_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} - [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_VBF_H_WW_enumunu_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV\to H\to WW@f$.
 */
class Hobs_VBF_H_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VBF_H_WW_enumunu_ATLAS13 constructor.
     */
    Hobs_VBF_H_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} / [\sigma_{VV\to H}\cdot BR(H\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_VBF_H_WW_enumunu_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$VV\to H\to WW@f$ assuming a Gaussian likelihood.
 */
class Robs_VBF_H_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_VBF_H_WW_enumunu_ATLAS13 constructor.
     */
    Robs_VBF_H_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{theo}} - [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_VBF_H_WW_lnulnu_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$(gg+VV)\to H\to WW\to \ell \nu \ell \nu@f$.
 */
class Hobs_ggF_VBF_H_WW_lnulnu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_VBF_H_WW_lnulnu_CMS13 constructor.
     */
    Hobs_ggF_VBF_H_WW_lnulnu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{(gg+VV)\to H}\cdot BR^{\text{THDM}}(H\to WW\to \ell \nu \ell \nu)]_{\text{theo}} / [\sigma_{(gg+VV)\to H}\cdot BR(H\to WW\to \ell \nu \ell \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_VBF_H_WW_lnulnu_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$(gg+VV)\to H\to WW\to \ell \nu \ell \nu@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_VBF_H_WW_lnulnu_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_VBF_H_WW_lnulnu_CMS13 constructor.
     */
    Robs_ggF_VBF_H_WW_lnulnu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{(gg+VV)\to H}\cdot BR^{\text{THDM}}(H\to WW\to \ell \nu \ell \nu)]_{\text{theo}} - [\sigma^{\text{THDM}}_{(gg+VV)\to H}\cdot BR^{\text{THDM}}(H\to WW\to \ell \nu \ell \nu)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{(gg+VV)\to H}\cdot BR^{\text{THDM}}(H\to WW\to \ell \nu \ell \nu)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_bbgaga_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H\to hh@f$.
 */
class Hobs_pp_H_hh_bbgaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_bbgaga_ATLAS13 constructor.
     */
    Hobs_pp_H_hh_bbgaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_hh_bbgaga_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$pp\to H\to hh@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_hh_bbgaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_hh_bbgaga_ATLAS13 constructor.
     */
    Robs_pp_H_hh_bbgaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_bbgaga_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H\to hh\to bb \gamma \gamma@f$.
 */
class Hobs_pp_H_hh_bbgaga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_bbgaga_CMS13 constructor.
     */
    Hobs_pp_H_hh_bbgaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to bb \gamma \gamma)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to bb \gamma \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_hh_bbgaga_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H\to hh\to bb \gamma \gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_hh_bbgaga_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_hh_bbgaga_CMS13 constructor.
     */
    Robs_pp_H_hh_bbgaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to bb \gamma \gamma)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to bb \gamma \gamma)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to bb \gamma \gamma)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_bbbb_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_pp_H_hh_bbbb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_bbbb_ATLAS13 constructor.
     */
    Hobs_pp_H_hh_bbbb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_hh_bbbb_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$pp\to H\to hh\to b\bar b b\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_hh_bbbb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_hh_bbbb_ATLAS13 constructor.
     */
    Robs_pp_H_hh_bbbb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_bbbb_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_pp_H_hh_bbbb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_bbbb_CMS13 constructor.
     */
    Hobs_pp_H_hh_bbbb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_hh_bbbb_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H\to hh\to b\bar b b\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_hh_bbbb_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_hh_bbbb_CMS13 constructor.
     */
    Robs_pp_H_hh_bbbb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_H_hh_gagaWW_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H\to hh@f$.
 */
class Hobs_ggF_H_hh_gagaWW_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H_hh_gagaWW_ATLAS13 constructor.
     */
    Hobs_ggF_H_hh_gagaWW_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_H_hh_gagaWW_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on signal strength of the process @f$pp\to H\to hh@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H_hh_gagaWW_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H_hh_gagaWW_ATLAS13 constructor.
     */
    Robs_ggF_H_hh_gagaWW_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_bbtautau_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b \tau \tau@f$.
 */
class Hobs_pp_H_hh_bbtautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_bbtautau_CMS13 constructor.
     */
    Hobs_pp_H_hh_bbtautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau \tau)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_hh_bbtautau_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H\to hh\to b\bar b \tau \tau@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_hh_bbtautau_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_hh_bbtautau_CMS13 constructor.
     */
    Robs_pp_H_hh_bbtautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau \tau)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau \tau)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_hh_bblnulnu_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu@f$.
 */
class Hobs_pp_H_hh_bblnulnu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_hh_bblnulnu_CMS13 constructor.
     */
    Hobs_pp_H_hh_bblnulnu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_hh_bblnulnu_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_hh_bblnulnu_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_hh_bblnulnu_CMS13 constructor.
     */
    Robs_pp_H_hh_bblnulnu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_H_bb_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H\to b\bar b@f$.
 */
class Hobs_pp_H_bb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H_bb_CMS13 constructor.
     */
    Hobs_pp_H_bb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{theo}} / [\sigma_{pp\to H}\cdot BR(H\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_H_bb_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H\to b\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H_bb_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H_bb_CMS13 constructor.
     */
    Robs_pp_H_bb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_tautau_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$ at 8 TeV.
 */
class log10_ggF_H_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_tautau_TH8 constructor.
     */
    log10_ggF_H_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_bbF_H_tautau_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to H\to \tau\tau@f$ at 8 TeV.
 */
class log10_bbF_H_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_bbF_H_tautau_TH8 constructor.
     */
    log10_bbF_H_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_gaga_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to \gamma\gamma@f$ at 8 TeV.
 */
class log10_pp_H_gaga_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H_gaga_TH8 constructor.
     */
    log10_pp_H_gaga_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_gaga_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to \gamma\gamma@f$ at 8 TeV.
 */
class log10_ggF_H_gaga_TH8: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_gaga_TH8 constructor.
     */
    log10_ggF_H_gaga_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma\gamma)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_mu_pp_H_VV_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the signal strength of the process @f$pp\to H\to VV@f$ at 8 TeV.
 */
class log10_mu_pp_H_VV_TH8: public ThObservable {
public:

    /**
     * @brief log10_mu_pp_H_VV_TH8 constructor.
     */
    log10_mu_pp_H_VV_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\mu_H^{\text{THDM}}(H\to VV)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_ZZ_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to ZZ@f$ at 8 TeV.
 */
class log10_ggF_H_ZZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_ZZ_TH8 constructor.
     */
    log10_ggF_H_ZZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_VBF_H_ZZ_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to H\to ZZ@f$ at 8 TeV.
 */
class log10_VBF_H_ZZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_VBF_H_ZZ_TH8 constructor.
     */
    log10_VBF_H_ZZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_WW_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to WW@f$ at 8 TeV.
 */
class log10_ggF_H_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_WW_TH8 constructor.
     */
    log10_ggF_H_WW_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_VBF_H_WW_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to H\to WW@f$ at 8 TeV.
 */
class log10_VBF_H_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_VBF_H_WW_TH8 constructor.
     */
    log10_VBF_H_WW_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_hh_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to hh@f$ at 8 TeV.
 */
class log10_ggF_H_hh_TH8: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_hh_TH8 constructor.
     */
    log10_ggF_H_hh_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh@f$ at 8 TeV.
 */
class log10_pp_H_hh_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_TH8 constructor.
     */
    log10_pp_H_hh_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_hh_bbtautau_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to hh\to b\bar b \tau\tau@f$ at 8 TeV.
 */
class log10_ggF_H_hh_bbtautau_TH8: public ThObservable {
 public:

  /**                                                                                                                                         
   * @brief log10_ggF_H_hh_bbtautau_TH8 constructor.                                                                                                                      
   */
  log10_ggF_H_hh_bbtautau_TH8(const StandardModel& SM_i);

  /**                                                                                                                                         
   * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau\tau)]@f$
   */
  double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_bbbb_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b b\bar b@f$ at 8 TeV.
 */
class log10_pp_H_hh_bbbb_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_bbbb_TH8 constructor.
     */
    log10_pp_H_hh_bbbb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_gagabb_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to \gamma\gamma b\bar b@f$ at 8 TeV.
 */
class log10_pp_H_hh_gagabb_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_gagabb_TH8 constructor.
     */
    log10_pp_H_hh_gagabb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma\gamma b\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_tt_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to t\bar t@f$ at 8 TeV.
 */
class log10_ggF_H_tt_TH8: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_tt_TH8 constructor.
     */
    log10_ggF_H_tt_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_bbF_H_bb_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to H\to b\bar b@f$ at 8 TeV.
 */
class log10_bbF_H_bb_TH8: public ThObservable {
public:

    /**
     * @brief log10_bbF_H_bb_TH8 constructor.
     */
    log10_bbF_H_bb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_tautau_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to \tau \tau@f$ at 13 TeV.
 */
class log10_ggF_H_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_tautau_TH13 constructor.
     */
    log10_ggF_H_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_bbF_H_tautau_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to H\to \tau \tau@f$ at 13 TeV.
 */
class log10_bbF_H_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_bbF_H_tautau_TH13 constructor.
     */
    log10_bbF_H_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to \tau \tau)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_gaga_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to \gamma \gamma@f$ at 13 TeV.
 */
class log10_pp_H_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_gaga_TH13 constructor.
     */
    log10_pp_H_gaga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_gaga_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to \gamma \gamma@f$ at 13 TeV.
 */
class log10_ggF_H_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_gaga_TH13 constructor.
     */
    log10_ggF_H_gaga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_Zga_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to Z\gamma@f$ at 13 TeV.
 */
class log10_pp_H_Zga_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_Zga_TH13 constructor.
     */
    log10_pp_H_Zga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to Z\gamma)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_ZZ_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to ZZ@f$ at 13 TeV.
 */
class log10_pp_H_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_ZZ_TH13 constructor.
     */
    log10_pp_H_ZZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_ZZ_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to ZZ@f$ at 13 TeV.
 */
class log10_ggF_H_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_ZZ_TH13 constructor.
     */
    log10_ggF_H_ZZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_VBF_H_ZZ_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to H\to ZZ@f$ at 13 TeV.
 */
class log10_VBF_H_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_VBF_H_ZZ_TH13 constructor.
     */
    log10_VBF_H_ZZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_ZZ_llll_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to ZZ\to \ell \ell \ell \ell@f$ at 13 TeV.
 */
class log10_ggF_H_ZZ_llll_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_ZZ_llll_TH13 constructor.
     */
    log10_ggF_H_ZZ_llll_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ\to \ell \ell \ell \ell)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_VBF_H_ZZ_llll_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to H\to ZZ\to \ell \ell \ell \ell@f$ at 13 TeV.
 */
class log10_VBF_H_ZZ_llll_TH13: public ThObservable {
public:

    /**
     * @brief log10_VBF_H_ZZ_llll_TH13 constructor.
     */
    log10_VBF_H_ZZ_llll_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to ZZ\to \ell \ell \ell \ell)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_ZZ_llll_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to ZZ\to \ell \ell \ell \ell@f$ at 13 TeV.
 */
class log10_pp_H_ZZ_llll_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_ZZ_llll_TH13 constructor.
     */
    log10_pp_H_ZZ_llll_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to ZZ\to \ell \ell \ell \ell)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_VBF_VH_H_ZZ_llll_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$(VV+VH)\to H\to ZZ\to \ell \ell \ell \ell@f$ at 13 TeV.
 */
class log10_VBF_VH_H_ZZ_llll_TH13: public ThObservable {
public:

    /**
     * @brief log10_VBF_VH_H_ZZ_llll_TH13 constructor.
     */
    log10_VBF_VH_H_ZZ_llll_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{(VV+VH)\to H}\cdot BR^{\text{THDM}}(H\to ZZ\to \ell \ell \ell \ell)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_WW_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to WW@f$ at 13 TeV.
 */
class log10_ggF_H_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_WW_TH13 constructor.
     */
    log10_ggF_H_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to WW)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_VBF_H_WW_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to H\to WW@f$ at 13 TeV.
 */
class log10_VBF_H_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_VBF_H_WW_TH13 constructor.
     */
    log10_VBF_H_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{VV\to H}\cdot BR^{\text{THDM}}(H\to WW)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_VBF_H_WW_lnulnu_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$(gg+VV)\to H\to WW\to \ell \nu \ell \nu@f$ at 13 TeV.
 */
class log10_ggF_VBF_H_WW_lnulnu_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggF_VBF_H_WW_lnulnu_TH13 constructor.
     */
    log10_ggF_VBF_H_WW_lnulnu_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{(gg+VV)\to H}\cdot BR^{\text{THDM}}(H\to WW\to \ell \nu \ell \nu)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_H_hh_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H\to hh@f$ at 13 TeV.
 */
class log10_ggF_H_hh_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggF_H_hh_TH13 constructor.
     */
    log10_ggF_H_hh_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to hh)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh@f$ at 13 TeV.
 */
class log10_pp_H_hh_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_TH13 constructor.
     */
    log10_pp_H_hh_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_bbbb_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b b\bar b@f$ at 13 TeV.
 */
class log10_pp_H_hh_bbbb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_bbbb_TH13 constructor.
     */
    log10_pp_H_hh_bbbb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_gagabb_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to \gamma \gamma b\bar b@f$ at 13 TeV.
 */
class log10_pp_H_hh_gagabb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_gagabb_TH13 constructor.
     */
    log10_pp_H_hh_gagabb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to \gamma \gamma b\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_bbtautau_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b \tau \tau@f$ at 13 TeV.
 */
class log10_pp_H_hh_bbtautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_bbtautau_TH13 constructor.
     */
    log10_pp_H_hh_bbtautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b \tau \tau)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_hh_bblnulnu_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu@f$ at 13 TeV.
 */
class log10_pp_H_hh_bblnulnu_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_hh_bblnulnu_TH13 constructor.
     */
    log10_pp_H_hh_bblnulnu_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to hh\to b\bar b WW\to b\bar b \ell \nu \ell \nu)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_tt_H_tt_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$t\bar t\to H\to t\bar t@f$ at 13 TeV.
 */
class log10_tt_H_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_tt_H_tt_TH13 constructor.
     */
    log10_tt_H_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{t\bar t\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_bb_H_tt_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to H\to t\bar t@f$ at 13 TeV.
 */
class log10_bb_H_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_H_tt_TH13 constructor.
     */
    log10_bb_H_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to H}\cdot BR^{\text{THDM}}(H\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_H_bb_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H\to b\bar b@f$ at 13 TeV.
 */
class log10_pp_H_bb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H_bb_TH13 constructor.
     */
    log10_pp_H_bb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H}\cdot BR^{\text{THDM}}(H\to b\bar b)]@f$
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







//
///**
// * @class LIMITmLIMEST
// * @ingroup THDM
// * @brief 
// */
//class LIMITmLIMEST : public ThObservable {
//public:
//    
//    /**
//     * @brief Constructor.
//     */
//    LIMITmLIMEST(const StandardModel& SM_i);
//    
//    /**
//     * @return 
//     */
//    double computeThValue ();
//private:
//    const THDM& myTHDM;
//};
//
///**
// * @class DEVIATIONoBANDSIZE
// * @ingroup THDM
// * @brief 
// */
//class DEVIATIONoBANDSIZE : public ThObservable {
//public:
//    
//    /**
//     * @brief Constructor.
//     */
//    DEVIATIONoBANDSIZE(const StandardModel& SM_i);
//    
//    /**
//     * @return 
//     */
//    double computeThValue ();
//private:
//    const THDM& myTHDM;
//};





#endif	/* HEAVYHIGGS_H */
