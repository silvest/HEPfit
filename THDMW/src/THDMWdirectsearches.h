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
 * @class Hobs_pp_Sr_tt_ATLAS13
 * @ingroup THDMW
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to Gkk\to t\bar t@f$.
 */
class Hobs_pp_Sr_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Sr_tt_ATLAS13 constructor.
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

/**
 * @class Hobs_pp_Srtt_tttt_ATLAS13
 * @ingroup THDMW
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to H t t\bar \to t\bar t  t\bar t@f$.
 */
class Hobs_pp_Srtt_tttt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Stt_tttt_ATLAS13 constructor.
     */
    Hobs_pp_Srtt_tttt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return xsection times Br ratio for pp -> Sr -> t tbar
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};


/**
 * @class log10_pp_Srtt_tttt_TH13
 * @ingroup THDMW
 * @brief Decadic logarithm of the cross section times branching ratio of the process pp -> Sr t tbar -> t tbar t tbar at 13 TeV.
 */
class log10_pp_Srtt_tttt_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_Srtt_tttt_TH13 constructor.
     */
    log10_pp_Srtt_tttt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDMW}}_{pp\to Sr t t\bar}\cdot BR^{\text{THDMW}}(Sr\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};







/**
 * @class Hobs_pp_Sr_jj_CMS13
 * @ingroup THDMW
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to Sr \to j j@f$.
 */
class Hobs_pp_Sr_jj_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_S_jj_ATLAS13 constructor.
     */
    Hobs_pp_Sr_jj_CMS13(const StandardModel& SM_i);

    /**
     * @return xsection times Br ratio for pp -> Sr -> j j
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};


/**
 * @class log10_pp_Sr_jj_TH13
 * @ingroup THDMW
 * @brief Decadic logarithm of the cross section times branching ratio of the process pp -> Sr  -> j j at 13 TeV.
 */
class log10_pp_Sr_jj_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_Sr_jj_TH13 constructor.
     */
    log10_pp_Sr_jj_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDMW}}_{pp\to Sr}\cdot BR^{\text{THDMW}}(Sr\to j j)]@f$
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};








/**
 * @class Hobs_pp_SrSr_jjjj_ATLAS13
 * @ingroup THDMW
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to Sr Sr \to j j j j@f$.
 */
class Hobs_pp_SrSr_jjjj_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_SS_jjjj_ATLAS13 constructor.
     */
    Hobs_pp_SrSr_jjjj_ATLAS13(const StandardModel& SM_i);

    /**
     * @return xsection times Br ratio for pp -> Sr Sr -> j j j j
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};


/**
 * @class log10_pp_SrSr_jjjj_TH13
 * @ingroup THDMW
 * @brief Decadic logarithm of the cross section times branching ratio of the process pp -> Sr Sr -> j j j j at 13 TeV.
 */
class log10_pp_SrSr_jjjj_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_SrSr_jjjj_TH13 constructor.
     */
    log10_pp_SrSr_jjjj_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDMW}}_{pp\to Sr Sr}\cdot BR^{\text{THDMW}}(Sr\to j j)]@f$
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class log10_pp_SrSr_jjjj_TH13
 * @ingroup THDMW
 * @brief Decadic logarithm of the cross section times branching ratio of the process pp -> Sr Sr -> j j j j at 13 TeV.
 */

//class logpp_SrSr_jjjj_TH13: public ThObservable {
//public:

    /**
     * @brief log10_pp_SrSr_jjjj_TH13 constructor.
     */
//    logpp_SrSr_jjjj_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{THDMW}}_{pp\to Sr Sr}\cdot BR^{\text{THDMW}}(Sr\to j j)]@f$
     */
//    double computeThValue();
//private:
//    const THDMW& myTHDMW;
//};





#endif	/* THDMWDIRECTSEARCHES_H */