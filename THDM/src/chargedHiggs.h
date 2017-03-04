/*
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHARGEDHIGGS_H
#define	CHARGEDHIGGS_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"
#include "THDMcache.h"

/**
 * @class chargedHiggs
 * @ingroup THDM 
 * @brief Base class for charged Higgs search observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */

/**
 * @class Hobs_pp_Hpm_taunu_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H^\pm \to \tau^\pm \nu@f$.
 */
class Hobs_pp_Hpm_taunu_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hpm_taunu_ATLAS8 constructor.
     */
    Hobs_pp_Hpm_taunu_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm \to \tau^\pm \nu)]_{\text{theo}} / [\sigma_{pp\to H^\pm}\cdot BR(H^\pm \to \tau^\pm \nu)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_Hpm_taunu_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$pp\to H^\pm \to \tau^\pm \nu@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_Hpm_taunu_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_Hpm_taunu_ATLAS8 constructor.
     */
    Robs_pp_Hpm_taunu_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm \to \tau^\pm \nu)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm \to \tau^\pm \nu)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{pp\to H^\pm}\cdot BR(H^\pm \to \tau^\pm \nu)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_Hp_taunu_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section of the process @f$pp\to H^+@f$.
 */
class Hobs_pp_Hp_taunu_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hp_taunu_CMS8 constructor.
     */
    Hobs_pp_Hp_taunu_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H^+}]_{\text{theo}} / [\sigma_{pp\to H^+}]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_Hp_taunu_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$pp\to H^+@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_Hp_taunu_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_Hp_taunu_CMS8 constructor.
     */
    Robs_pp_Hp_taunu_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H^+}]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H^+}]_{\text{CMS,95\% observed}} \right) / [\sigma_{pp\to H^+}]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_Hp_tb_ATLAS8
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H^+ \to t\bar b@f$.
 */
class Hobs_pp_Hp_tb_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hp_tb_ATLAS8 constructor.
     */
    Hobs_pp_Hp_tb_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{theo}} / [\sigma_{pp\to H^+}\cdot BR(H^+ \to t\bar b)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_Hp_tb_ATLAS8
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$pp\to H^+ \to t\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_Hp_tb_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_Hp_tb_ATLAS8 constructor.
     */
    Robs_pp_Hp_tb_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{pp\to H^+}\cdot BR(H^+ \to t\bar b)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_Hp_tb_CMS8
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section of the process @f$pp\to H^+@f$.
 */
class Hobs_pp_Hp_tb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hp_tb_CMS8 constructor.
     */
    Hobs_pp_Hp_tb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H^+}]_{\text{theo}} / [\sigma_{pp\to H^+}]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_Hp_tb_CMS8
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$pp\to H^+@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_Hp_tb_CMS8: public ThObservable {
public:

    /**
     * @brief Robs_pp_Hp_tb_CMS8 constructor.
     */
    Robs_pp_Hp_tb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H^+}]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H^+}]_{\text{CMS,95\% observed}} \right) / [\sigma_{pp\to H^+}]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_Hpm_taunu_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H^\pm \to \tau^\pm \nu@f$.
 */
class Hobs_pp_Hpm_taunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hpm_taunu_ATLAS13 constructor.
     */
    Hobs_pp_Hpm_taunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm \to \tau^\pm \nu)]_{\text{theo}} / [\sigma_{pp\to H^\pm}\cdot BR(H^\pm \to \tau^\pm \nu)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_Hpm_taunu_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$pp\to H^\pm \to \tau^\pm \nu@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_Hpm_taunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_Hpm_taunu_ATLAS13 constructor.
     */
    Robs_pp_Hpm_taunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm \to \tau^\pm \nu)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm \to \tau^\pm \nu)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{pp\to H^\pm}\cdot BR(H^\pm \to \tau^\pm \nu)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_Hpm_taunu_CMS13
 * @ingroup THDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H^\pm \to \tau^\pm \nu@f$.
 */
class Hobs_pp_Hpm_taunu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hpm_taunu_CMS13 constructor.
     */
    Hobs_pp_Hpm_taunu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm \to \tau^\pm \nu)]_{\text{theo}} / [\sigma_{pp\to H^\pm}\cdot BR(H^\pm \to \tau^\pm \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_Hpm_taunu_CMS13
 * @ingroup THDM
 * @brief Observable for the implementation of the CMS upper limit on the process @f$pp\to H^\pm \to \tau^\pm \nu@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_Hpm_taunu_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_Hpm_taunu_CMS13 constructor.
     */
    Robs_pp_Hpm_taunu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm \to \tau^\pm \nu)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm \to \tau^\pm \nu)]_{\text{CMS,95\% observed}} \right) / [\sigma_{pp\to H^\pm}\cdot BR(H^\pm \to \tau^\pm \nu)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_Hp_tb_ATLAS13
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H^+ \to t\bar b@f$.
 */
class Hobs_pp_Hp_tb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hp_tb_ATLAS13 constructor.
     */
    Hobs_pp_Hp_tb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{theo}} / [\sigma_{pp\to H^+}\cdot BR(H^+ \to t\bar b)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_Hp_tb_ATLAS13
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$pp\to H^+ \to t\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_Hp_tb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_Hp_tb_ATLAS13 constructor.
     */
    Robs_pp_Hp_tb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{pp\to H^+}\cdot BR(H^+ \to t\bar b)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_Hp_tb_ATLAS13_1
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H^+ \to t\bar b@f$.
 */
class Hobs_pp_Hp_tb_ATLAS13_1: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hp_tb_ATLAS13_1 constructor.
     */
    Hobs_pp_Hp_tb_ATLAS13_1(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{theo}} / [\sigma_{pp\to H^+}\cdot BR(H^+ \to t\bar b)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_Hp_tb_ATLAS13_1
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$pp\to H^+ \to t\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_Hp_tb_ATLAS13_1: public ThObservable {
public:

    /**
     * @brief Robs_pp_Hp_tb_ATLAS13_1 constructor.
     */
    Robs_pp_Hp_tb_ATLAS13_1(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{pp\to H^+}\cdot BR(H^+ \to t\bar b)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Hobs_pp_Hp_tb_ATLAS13_2
 * @ingroup THDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H^+ \to t\bar b@f$.
 */
class Hobs_pp_Hp_tb_ATLAS13_2: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hp_tb_ATLAS13_2 constructor.
     */
    Hobs_pp_Hp_tb_ATLAS13_2(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{theo}} / [\sigma_{pp\to H^+}\cdot BR(H^+ \to t\bar b)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Robs_pp_Hp_tb_ATLAS13_2
 * @ingroup THDM
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$pp\to H^+ \to t\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_Hp_tb_ATLAS13_2: public ThObservable {
public:

    /**
     * @brief Robs_pp_Hp_tb_ATLAS13_2 constructor.
     */
    Robs_pp_Hp_tb_ATLAS13_2(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{theo}} - [\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{pp\to H^+}\cdot BR(H^+ \to t\bar b)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_Hpm_taunu_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^\pm\to \tau^\pm \nu@f$ at 8 TeV.
 */
class log10_pp_Hpm_taunu_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_Hpm_taunu_TH8 constructor.
     */
    log10_pp_Hpm_taunu_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm\to \tau^\pm \nu)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_Hp_tb_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^+\to t\bar b@f$ at 8 TeV.
 */
class log10_pp_Hp_tb_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_Hp_tb_TH8 constructor.
     */
    log10_pp_Hp_tb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_Hp_TH8
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section of the process @f$pp\to H^+@f$ at 8 TeV.
 */
class log10_pp_Hp_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_Hp_TH8 constructor.
     */
    log10_pp_Hp_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H^+}]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_Hpm_taunu_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^\pm\to \tau^\pm \nu@f$ at 13 TeV.
 */
class log10_pp_Hpm_taunu_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_Hpm_taunu_TH13 constructor.
     */
    log10_pp_Hpm_taunu_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H^\pm}\cdot BR^{\text{THDM}}(H^\pm\to \tau^\pm \nu)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class log10_pp_Hp_tb_TH13
 * @ingroup THDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^+\to t\bar b@f$ at 13 TeV.
 */
class log10_pp_Hp_tb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_Hp_tb_TH13 constructor.
     */
    log10_pp_Hp_tb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDM}}_{pp\to H^+}\cdot BR^{\text{THDM}}(H^+ \to t\bar b)]@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class Gamma_Hp_THDM
 * @ingroup THDM
 * @brief Total Hp decay rate in the %THDM.
 */
class Gamma_Hp_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Gamma_Hp_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$\Gamma_H^+@f$ in units of GeV
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

#endif	/* CHARGEDHIGGS_H */
