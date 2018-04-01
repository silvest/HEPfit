/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GMDIRECTSEARCHES_H
#define	GMDIRECTSEARCHES_H

#include <stdexcept>
#include "ThObservable.h"
#include "GeorgiMachacek.h"
#include "GMcache.h"

/**
 * @class GMDirectSearches
 * @ingroup GeorgiMachacek
 * @brief Base class for direct Georgi-Machacek Higgs search observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */

/**
 * @class BR_H1_hh_GM
 * @ingroup GeorgiMachacek
 * @brief Branching ratio of @f$H_1@f$ to two @f$h@f$ in the GeorgiMachacek model.
 */
class BR_H1_hh_GM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_H1_hh_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1\to hh)@f$
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_ggF_H1_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H_1\to \tau\tau@f$.
 */
class Hobs_ggF_H1_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_ggF_H1_tautau_ATLAS8 constructor.
     */
    Hobs_ggF_H1_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H_1}\cdot BR(H_1\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Robs_ggF_H1_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg\to H_1\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_ggF_H1_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_ggF_H1_tautau_ATLAS8 constructor.
     */
    Robs_ggF_H1_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{theo}} - [\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{gg\to H_1}\cdot BR(H_1\to \tau\tau)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bbbb_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H_1\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_pp_H1_hh_bbbb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bbbb_CMS13 constructor.
     */
    Hobs_pp_H1_hh_bbbb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Robs_pp_H1_hh_bbbb_CMS13
 * @ingroup GeorgiMachacek
 * @brief Observable for the implementation of the CMS upper limit on signal strength of the process @f$pp\to H_1\to hh\to b\bar b b\bar b@f$ assuming a Gaussian likelihood.
 */
class Robs_pp_H1_hh_bbbb_CMS13: public ThObservable {
public:

    /**
     * @brief Robs_pp_H1_hh_bbbb_CMS13 constructor.
     */
    Robs_pp_H1_hh_bbbb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b b\bar b)]_{\text{theo}} - [\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b b\bar b)]_{\text{CMS,95\% observed}} \right) / [\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b b\bar b)]_{\text{CMS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_ggF_H1_tautau_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H_1\to \tau\tau@f$ at 8 TeV.
 */
class log10_ggF_H1_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_ggF_H1_tautau_TH8 constructor.
     */
    log10_ggF_H1_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_hh_bbbb_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H_1\to hh\to b\bar b b\bar b@f$ at 13 TeV.
 */
class log10_pp_H1_hh_bbbb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_hh_bbbb_TH13 constructor.
     */
    log10_pp_H1_hh_bbbb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

#endif	/* GMDIRECTSEARCHES_H */
