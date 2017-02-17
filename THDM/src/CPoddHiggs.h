/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CPODDHIGGS_H
#define	CPODDHIGGS_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"
#include "THDMcache.h"

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

/**
 * @class Hobs_ggF_A_tautau_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to A\to \tau\tau@f$.
 */
class Hobs_ggF_A_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_A_tautau_ATLAS8 constructor.
     */
    Hobs_ggF_A_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_A_tautau_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg\to A\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_A_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_A_tautau_ATLAS8 constructor.
     */
    Robs_ggF_A_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} - [\sigma_{gg\to A}\cdot BR(A\to \tau\tau)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{gg\to A}\cdot BR(A\to \tau\tau)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_A_tautau_CMS
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to A\to \tau\tau@f$.
 */
class Hobs_ggF_A_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_A_tautau_CMS8 constructor.
     */
    Hobs_ggF_A_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_A_tautau_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$gg\to A\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_A_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_A_tautau_CMS8 constructor.
     */
    Robs_ggF_A_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} - [\sigma_{gg\to A}\cdot BR(A\to \tau\tau)]_{\text{CMS,95\% observed}} \right) / [\sigma_{gg\to A}\cdot BR(A\to \tau\tau)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_A_tautau_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to A\to \tau\tau@f$.
 */
class Hobs_bbF_A_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_A_tautau_ATLAS8 constructor.
     */
    Hobs_bbF_A_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to A}\cdot BR(A\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_bbF_A_tautau_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$b\bar b\to A\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_bbF_A_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_bbF_A_tautau_ATLAS8 constructor.
     */
    Robs_bbF_A_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} - [\sigma_{b\bar b\to A}\cdot BR(A\to \tau\tau)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{b\bar b\to A}\cdot BR(A\to \tau\tau)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_A_tautau_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to A\to \tau\tau@f$.
 */
class Hobs_bbF_A_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_A_tautau_CMS8 constructor.
     */
    Hobs_bbF_A_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to A}\cdot BR(A\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_bbF_A_tautau_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$b\bar b\to A\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_bbF_A_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_bbF_A_tautau_CMS8 constructor.
     */
    Robs_bbF_A_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]_{\text{theo}} - [\sigma_{b\bar b\to A}\cdot BR(A\to \tau\tau)]_{\text{CMS,95\% observed}} \right) / [\sigma_{b\bar b\to A}\cdot BR(A\to \tau\tau)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_A_gaga_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp \to A\to \gamma\gamma@f$.
 */
class Hobs_pp_A_gaga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_A_gaga_ATLAS8 constructor.
     */
    Hobs_pp_A_gaga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)]_{\text{theo}} / [\sigma_{pp\to A}\cdot BR(A\to \gamma\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_A_gaga_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$pp \to A\to \gamma\gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_A_gaga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_A_gaga_ATLAS8 constructor.
     */
    Robs_pp_A_gaga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)]_{\text{theo}} - [\sigma_{pp\to A}\cdot BR(A\to \gamma\gamma)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{pp\to A}\cdot BR(A\to \gamma\gamma)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_A_gaga_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to A\to \gamma\gamma@f$.
 */
class Hobs_ggF_A_gaga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_A_gaga_CMS8 constructor.
     */
    Hobs_ggF_A_gaga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to \gamma\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_A_gaga_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$gg \to A\to \gamma\gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_A_gaga_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_A_gaga_CMS8 constructor.
     */
    Robs_ggF_A_gaga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)]_{\text{theo}} - [\sigma_{gg\to A}\cdot BR(A\to \gamma\gamma)]_{\text{CMS,95\% observed}} \right) / [\sigma_{gg\to A}\cdot BR(A\to \gamma\gamma)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_A_Zga_llga_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to Z\gamma \to \ell \ell \gamma@f$.
 */
class Hobs_pp_A_Zga_llga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_A_Zga_llga_CMS8 constructor.
     */
    Hobs_pp_A_Zga_llga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to Z\gamma \to \ell \ell \gamma )]_{\text{theo}} / [\sigma_{pp\to A}\cdot BR(A\to Z\gamma \to \ell \ell \gamma )]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_A_Zga_llga_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$gg \to A\to Z\gamma \to \ell \ell \gamma@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_A_Zga_llga_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_A_Zga_llga_CMS8 constructor.
     */
    Robs_pp_A_Zga_llga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to Z\gamma \to \ell \ell \gamma )]_{\text{theo}} - [\sigma_{pp\to A}\cdot BR(A\to Z\gamma \to \ell \ell \gamma )]_{\text{CMS,95\% observed}} \right) / [\sigma_{pp\to A}\cdot BR(A\to Z\gamma \to \ell \ell \gamma )]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_A_hZ_bbll_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to A\to hZ \to b\bar b \ell \ell@f$.
 */
class Hobs_ggF_A_hZ_bbll_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_A_hZ_bbll_CMS8 constructor.
     */
    Hobs_ggF_A_hZ_bbll_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b \ell \ell)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to hZ \to b\bar b \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_A_hZ_bbll_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$gg \to A\to hZ \to b\bar b \ell \ell@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_A_hZ_bbll_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_A_hZ_bbll_CMS8 constructor.
     */
    Robs_ggF_A_hZ_bbll_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b \ell \ell)]_{\text{theo}} - [\sigma_{gg\to A}\cdot BR(A\to hZ \to b\bar b \ell \ell)]_{\text{CMS,95\% observed}} \right) / [\sigma_{gg\to A}\cdot BR(A\to hZ \to b\bar b \ell \ell)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_A_hZ_bbZ_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to A\to hZ \to b\bar b Z@f$.
 */
class Hobs_ggF_A_hZ_bbZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_A_hZ_bbZ_ATLAS8 constructor.
     */
    Hobs_ggF_A_hZ_bbZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b Z)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to hZ \to b\bar b Z)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_A_hZ_bbZ_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg \to A\to hZ \to b\bar b Z@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_A_hZ_bbZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_A_hZ_bbZ_ATLAS8 constructor.
     */
    Robs_ggF_A_hZ_bbZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to b\bar b Z)]_{\text{theo}} - [\sigma_{gg\to A}\cdot BR(A\to hZ \to b\bar b Z)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{gg\to A}\cdot BR(A\to hZ \to b\bar b Z)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_A_hZ_tautaull_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg \to A\to hZ \to \tau \tau \ell \ell@f$.
 */
class Hobs_ggF_A_hZ_tautaull_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_A_hZ_tautaull_CMS8 constructor.
     */
    Hobs_ggF_A_hZ_tautaull_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau \ell \ell)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to hZ \to \tau \tau \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_A_hZ_tautaull_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$gg \to A\to hZ \to \tau \tau \ell \ell@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_A_hZ_tautaull_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_A_hZ_tautaull_CMS8 constructor.
     */
    Robs_ggF_A_hZ_tautaull_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau \ell \ell)]_{\text{theo}} - [\sigma_{gg\to A}\cdot BR(A\to hZ \to \tau \tau \ell \ell)]_{\text{CMS,95\% observed}} \right) / [\sigma_{gg\to A}\cdot BR(A\to hZ \to \tau \tau \ell \ell)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_A_hZ_tautauZ_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to A\to hZ \to \tau \tau Z@f$.
 */
class Hobs_ggF_A_hZ_tautauZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_A_hZ_tautauZ_ATLAS8 constructor.
     */
    Hobs_ggF_A_hZ_tautauZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau b Z)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to hZ \to \tau \tau Z)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_A_hZ_tautauZ_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg \to A\to hZ \to \tau \tau Z@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_A_hZ_tautauZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_A_hZ_tautauZ_ATLAS8 constructor.
     */
    Robs_ggF_A_hZ_tautauZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ \to \tau \tau b Z)]_{\text{theo}} - [\sigma_{gg\to A}\cdot BR(A\to hZ \to \tau \tau Z)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{gg\to A}\cdot BR(A\to hZ \to \tau \tau Z)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_ggF_A_tt_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to A\to t\bar t@f$.
 */
class Hobs_ggF_A_tt_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_A_tt_ATLAS8 constructor.
     */
    Hobs_ggF_A_tt_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to t\bar t)]_{\text{theo}} / [\sigma_{gg\to A}\cdot BR(A\to t\bar t))_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_ggF_A_tt_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg \to A\to t\bar t@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_A_tt_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_A_tt_ATLAS8 constructor.
     */
    Robs_ggF_A_tt_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to t\bar t)]_{\text{theo}} - [\sigma_{gg\to A}\cdot BR(A\to t\bar t)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{gg\to A}\cdot BR(A\to t\bar t)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_bbF_A_bb_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b \to A\to b\bar b@f$.
 */
class Hobs_bbF_A_bb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bbF_A_bb_CMS8 constructor.
     */
    Hobs_bbF_A_bb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to A}\cdot BR(A\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_bbF_A_bb_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$b\bar b \to A\to b\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_bbF_A_bb_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_bbF_A_bb_CMS8 constructor.
     */
    Robs_bbF_A_bb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to b\bar b)]_{\text{theo}} - [\sigma_{b\bar b\to A}\cdot BR(A\to b\bar b)]_{\text{CMS,95\% observed}} \right) / [\sigma_{b\bar b\to A}\cdot BR(A\to b\bar b)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_A_tautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to \tau\tau@f$.
 */
class log10_ggF_A_tautau_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_A_tautau_TH constructor.
     */
    log10_ggF_A_tautau_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_bbF_A_tautau_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to A\to \tau\tau@f$.
 */
class log10_bbF_A_tautau_TH: public ThObservable {
public:

    /**
     * @brief log10_bbF_A_tautau_TH constructor.
     */
    log10_bbF_A_tautau_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_A_gaga_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to A\to \gamma\gamma@f$.
 */
class log10_pp_A_gaga_TH: public ThObservable {
public:

    /**
     * @brief log10_pp_A_gaga_TH constructor.
     */
    log10_pp_A_gaga_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_A_gaga_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to \gamma\gamma@f$.
 */
class log10_ggF_A_gaga_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_A_gaga_TH constructor.
     */
    log10_ggF_A_gaga_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to \gamma\gamma)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_A_Zga_llga_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to Z\gamma \to \ell \ell \gamma@f$.
 */
class log10_pp_A_Zga_llga_TH: public ThObservable {
public:

    /**
     * @brief log10_pp_A_Zga_llga_TH constructor.
     */
    log10_pp_A_Zga_llga_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to A}\cdot BR^{\text{THDM}}(A\to Z\gamma \to \ell \ell \gamma )]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_A_hZ_bbll_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to hZ \to b\bar b \ell \ell@f$.
 */
class log10_ggF_A_hZ_bbll_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_A_hZ_bbll_TH constructor.
     */
    log10_ggF_A_hZ_bbll_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ\to b\bar b \ell \ell)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_A_hZ_bbZ_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to hZ \to b\bar b Z@f$.
 */
class log10_ggF_A_hZ_bbZ_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_A_hZ_bbZ_TH constructor.
     */
    log10_ggF_A_hZ_bbZ_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ\to b\bar b Z)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_A_hZ_tautaull_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to hZ \to \tau \tau \ell \ell@f$.
 */
class log10_ggF_A_hZ_tautaull_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_A_hZ_tautaull_TH constructor.
     */
    log10_ggF_A_hZ_tautaull_TH(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ\to \tau \tau \ell \ell)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_A_hZ_tautauZ_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to hZ \to \tau \tau Z@f$.
 */
class log10_ggF_A_hZ_tautauZ_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_A_hZ_tautauZ_TH constructor.
     */
    log10_ggF_A_hZ_tautauZ_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to hZ\to \tau \tau Z)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_ggF_A_tt_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to A\to t\bar t@f$.
 */
class log10_ggF_A_tt_TH: public ThObservable {
public:

    /**
     * @brief log10_ggF_A_tt_TH constructor.
     */
    log10_ggF_A_tt_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{gg\to A}\cdot BR^{\text{THDM}}(A\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_bbF_A_bb_TH
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to A\to b\bar b@f$.
 */
class log10_bbF_A_bb_TH: public ThObservable {
public:

    /**
     * @brief log10_bbF_A_bb_TH constructor.
     */
    log10_bbF_A_bb_TH(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{b\bar b\to A}\cdot BR^{\text{THDM}}(A\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Gamma_A_THDM
 * @ingroup THDM
 * @brief Total A decay rate in the %THDM.
 */
class Gamma_A_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Gamma_A_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$\Gamma_A@f$ in units of GeV
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

///**
// * @class rA_gaga_THDM
// * @ingroup THDM
// * @brief Squared relative coupling of @f$A@f$ to two photons.
// */
//class rA_gaga_THDM : public ThObservable {
//public:
//    
//    /**
//     * @brief Constructor.
//     */
//    rA_gaga_THDM(const StandardModel& SM_i);
//    
//    /**
//     * @return @f$r^{(A)}_{\gamma \gamma}@f$
//     */
//    double computeThValue ();
//};

/**
 * @class rA_gg_THDM
 * @ingroup THDM
 * @brief Squared relative coupling of @f$A@f$ to two gluons.
 */
class rA_gg_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor for the squared relative coupling of @f$A@f$ to two gluons.
     */
    rA_gg_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(A)}_{gg}@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class BR_A_HZ_THDM
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$A@f$ to an @f$H@f$ and a @f$Z@f$ boson.
 */
class BR_A_HZ_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_A_HZ_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(A\to HZ)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class BR_A_hZ_THDM
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$A@f$ to an @f$h@f$ and a @f$Z@f$ boson.
 */
class BR_A_hZ_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_A_hZ_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(A\to hZ)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class BR_A_HpW_THDM
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$A@f$ to a charged Higgs bosons and a @f$W@f$ boson.
 */
class BR_A_HpW_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_A_HpW_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(A\to H^\pm W^\mp)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @}
 */

#endif	/* CPODDHIGGS_H */
