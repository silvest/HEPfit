/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWDIRECTSEARCHES_H
#define	THDMWDIRECTSEARCHES_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDMW.h"
#include "THDMWcache.h"

/**
 * @class THDMWdirectsearches
 * @ingroup THDMW 
 * @brief Base class for direct heavy Higgs search observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */

/**
 * @class Hobs_ggF_H_tautau_ATLAS8
 * @ingroup THDMW
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H\to \tau\tau@f$.
 */
class Hobs_pp_Sr_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Z'_ttbar_ATLAS13 constructor.
     */
    Hobs_pp_Sr_tt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return xsection times Br ratio for pp -> Sr -> t tbar
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};


/**
 * @class log10_pp_Sr_tt_TH13
 * @ingroup THDMW
 * @brief Decadic logarithm of the cross section times branching ratio of the process pp -> Sr -> t tbar at 13 TeV.
 */
class log10_pp_Sr_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_Sr_tt_TH13 constructor.
     */
    log10_pp_Sr_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDMW}}_{pp\to Sr}\cdot BR^{\text{THDMW}}(Sr\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};







#endif	/* THDMWDIRECTSEARCHES_H */