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

class Hobs_tt_H1_tt_ATLAS13: public ThObservable {
public:

    Hobs_tt_H1_tt_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H1_tt_ATLAS13: public ThObservable {
public:

    Hobs_bb_H1_tt_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_tt_H3_tt_ATLAS13: public ThObservable {
public:

    Hobs_tt_H3_tt_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H3_tt_ATLAS13: public ThObservable {
public:

    Hobs_bb_H3_tt_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H1_bb_CMS8: public ThObservable {
public:

    Hobs_bb_H1_bb_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_bb_CMS8: public ThObservable {
public:

    Hobs_gg_H1_bb_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_bb_CMS13: public ThObservable {
public:

    Hobs_pp_H1_bb_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H1_bb_CMS13: public ThObservable {
public:

    Hobs_bb_H1_bb_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H3_bb_CMS8: public ThObservable {
public:

    Hobs_bb_H3_bb_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_bb_CMS8: public ThObservable {
public:

    Hobs_gg_H3_bb_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3_bb_CMS13: public ThObservable {
public:

    Hobs_pp_H3_bb_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H3_bb_CMS13: public ThObservable {
public:

    Hobs_bb_H3_bb_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_tautau_CMS8: public ThObservable {
public:

    Hobs_gg_H1_tautau_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H1_tautau_CMS8: public ThObservable {
public:

    Hobs_bb_H1_tautau_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_tautau_ATLAS13: public ThObservable {
public:

    Hobs_gg_H1_tautau_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_tautau_CMS13: public ThObservable {
public:

    Hobs_gg_H1_tautau_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H1_tautau_ATLAS13: public ThObservable {
public:

    Hobs_bb_H1_tautau_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H1_tautau_CMS13: public ThObservable {
public:

    Hobs_bb_H1_tautau_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to H_1\to \tau\tau@f$.
 */
class Hobs_gg_H1_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_tautau_ATLAS8 constructor.
     */
    Hobs_gg_H1_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H_1}\cdot BR(H_1\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Robs_gg_H1_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$gg\to H_1\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_gg_H1_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_gg_H1_tautau_ATLAS8 constructor.
     */
    Robs_gg_H1_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{theo}} - [\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{gg\to H_1}\cdot BR(H_1\to \tau\tau)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H1_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to H_1\to \tau\tau@f$.
 */
class Hobs_bb_H1_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H1_tautau_ATLAS8 constructor.
     */
    Hobs_bb_H1_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to H_1}\cdot BR(H_1\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Robs_bb_H1_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Observable for the implementation of the ATLAS upper limit on the process @f$b\bar b\to H_1\to \tau\tau@f$ assuming a Gaussian likelihood.
 */
class Robs_bb_H1_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Robs_bb_H1_tautau_ATLAS8 constructor.
     */
    Robs_bb_H1_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{theo}} - [\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{ATLAS,95\% observed}} \right) / [\sigma_{b\bar b\to H_1}\cdot BR(H_1\to \tau\tau)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_tautau_ATLAS8: public ThObservable {
public:

    Hobs_gg_H3_tautau_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_tautau_CMS8: public ThObservable {
public:

    Hobs_gg_H3_tautau_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H3_tautau_ATLAS8: public ThObservable {
public:

    Hobs_bb_H3_tautau_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H3_tautau_CMS8: public ThObservable {
public:

    Hobs_bb_H3_tautau_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_tautau_ATLAS13: public ThObservable {
public:

    Hobs_gg_H3_tautau_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_tautau_CMS13: public ThObservable {
public:

    Hobs_gg_H3_tautau_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H3_tautau_ATLAS13: public ThObservable {
public:

    Hobs_bb_H3_tautau_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H3_tautau_CMS13: public ThObservable {
public:

    Hobs_bb_H3_tautau_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_gaga_ATLAS8: public ThObservable {
public:

    Hobs_gg_H1_gaga_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_gaga_ATLAS13: public ThObservable {
public:

    Hobs_pp_H1_gaga_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_gaga_CMS13: public ThObservable {
public:

    Hobs_gg_H1_gaga_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_gaga_ATLAS8: public ThObservable {
public:

    Hobs_gg_H3_gaga_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3_gaga_ATLAS13: public ThObservable {
public:

    Hobs_pp_H3_gaga_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_gaga_CMS13: public ThObservable {
public:

    Hobs_gg_H3_gaga_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_gaga_ATLAS8: public ThObservable {
public:

    Hobs_gg_H5_gaga_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_gaga_ATLAS13: public ThObservable {
public:

    Hobs_pp_H5_gaga_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_gaga_CMS13: public ThObservable {
public:

    Hobs_gg_H5_gaga_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_Zga_llga_ATLAS8: public ThObservable {
public:

    Hobs_pp_H1_Zga_llga_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_Zga_llga_CMS8: public ThObservable {
public:

    Hobs_pp_H1_Zga_llga_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_Zga_llga_ATLAS13: public ThObservable {
public:

    Hobs_gg_H1_Zga_llga_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_Zga_CMS13: public ThObservable {
public:

    Hobs_gg_H1_Zga_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3_Zga_llga_ATLAS8: public ThObservable {
public:

    Hobs_pp_H3_Zga_llga_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3_Zga_llga_CMS8: public ThObservable {
public:

    Hobs_pp_H3_Zga_llga_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_Zga_llga_ATLAS13: public ThObservable {
public:

    Hobs_gg_H3_Zga_llga_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_Zga_CMS13: public ThObservable {
public:

    Hobs_gg_H3_Zga_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_Zga_llga_ATLAS8: public ThObservable {
public:

    Hobs_pp_H5_Zga_llga_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_Zga_llga_CMS8: public ThObservable {
public:

    Hobs_pp_H5_Zga_llga_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_Zga_llga_ATLAS13: public ThObservable {
public:

    Hobs_gg_H5_Zga_llga_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_Zga_CMS13: public ThObservable {
public:

    Hobs_gg_H5_Zga_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_ZZ_ATLAS8: public ThObservable {
public:

    Hobs_gg_H1_ZZ_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H1_ZZ_ATLAS8: public ThObservable {
public:

    Hobs_VV_H1_ZZ_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    Hobs_gg_H1_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H1_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    Hobs_VV_H1_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_ZZ_qqllnunu_ATLAS13: public ThObservable {
public:

    Hobs_gg_H1_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H1_ZZ_qqllnunu_ATLAS13: public ThObservable {
public:

    Hobs_VV_H1_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_ZZ_llqqnunull_CMS13: public ThObservable {
public:

    Hobs_pp_H1_ZZ_llqqnunull_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H1_ZZ_llqqnunull_CMS13: public ThObservable {
public:

    Hobs_VV_H1_ZZ_llqqnunull_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_ZZ_qqnunu_CMS13: public ThObservable {
public:

    Hobs_pp_H1_ZZ_qqnunu_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_ZZ_ATLAS8: public ThObservable {
public:

    Hobs_gg_H5_ZZ_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H5_ZZ_ATLAS8: public ThObservable {
public:

    Hobs_VV_H5_ZZ_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    Hobs_gg_H5_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H5_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    Hobs_VV_H5_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_ZZ_qqllnunu_ATLAS13: public ThObservable {
public:

    Hobs_gg_H5_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H5_ZZ_qqllnunu_ATLAS13: public ThObservable {
public:

    Hobs_VV_H5_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_ZZ_llqqnunull_CMS13: public ThObservable {
public:

    Hobs_pp_H5_ZZ_llqqnunull_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H5_ZZ_llqqnunull_CMS13: public ThObservable {
public:

    Hobs_VV_H5_ZZ_llqqnunull_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_ZZ_qqnunu_CMS13: public ThObservable {
public:

    Hobs_pp_H5_ZZ_qqnunu_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_WW_ATLAS8: public ThObservable {
public:

    Hobs_gg_H1_WW_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H1_WW_ATLAS8: public ThObservable {
public:

    Hobs_VV_H1_WW_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_WW_enumunu_ATLAS13: public ThObservable {
public:

    Hobs_gg_H1_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H1_WW_enumunu_ATLAS13: public ThObservable {
public:

    Hobs_VV_H1_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_WW_lnuqq_ATLAS13: public ThObservable {
public:

    Hobs_gg_H1_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H1_WW_lnuqq_ATLAS13: public ThObservable {
public:

    Hobs_VV_H1_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_ggVV_H1_WW_lnulnu_CMS13: public ThObservable {
public:

    Hobs_ggVV_H1_WW_lnulnu_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_WW_lnuqq_CMS13: public ThObservable {
public:

    Hobs_pp_H1_WW_lnuqq_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_WW_ATLAS8: public ThObservable {
public:

    Hobs_gg_H5_WW_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H5_WW_ATLAS8: public ThObservable {
public:

    Hobs_VV_H5_WW_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_WW_enumunu_ATLAS13: public ThObservable {
public:

    Hobs_gg_H5_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H5_WW_enumunu_ATLAS13: public ThObservable {
public:

    Hobs_VV_H5_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_WW_lnuqq_ATLAS13: public ThObservable {
public:

    Hobs_gg_H5_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_VV_H5_WW_lnuqq_ATLAS13: public ThObservable {
public:

    Hobs_VV_H5_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_ggVV_H5_WW_lnulnu_CMS13: public ThObservable {
public:

    Hobs_ggVV_H5_WW_lnulnu_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_WW_lnuqq_CMS13: public ThObservable {
public:

    Hobs_pp_H5_WW_lnuqq_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_mu_pp_H1_VV_CMS8: public ThObservable {
public:

    Hobs_mu_pp_H1_VV_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_VV_qqqq_ATLAS13: public ThObservable {
public:

    Hobs_pp_H1_VV_qqqq_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_mu_pp_H5_VV_CMS8: public ThObservable {
public:

    Hobs_mu_pp_H5_VV_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_VV_qqqq_ATLAS13: public ThObservable {
public:

    Hobs_pp_H5_VV_qqqq_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_hh_ATLAS8: public ThObservable {
public:

    Hobs_gg_H1_hh_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_hh_bbbb_CMS8: public ThObservable {
public:

    Hobs_pp_H1_hh_bbbb_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_hh_gagabb_CMS8: public ThObservable {
public:

    Hobs_pp_H1_hh_gagabb_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_hh_bbtautau_CMS8: public ThObservable {
public:

    Hobs_gg_H1_hh_bbtautau_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_hh_bbtautau_CMS8: public ThObservable {
public:

    Hobs_pp_H1_hh_bbtautau_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_hh_bbbb_ATLAS13: public ThObservable {
public:

    Hobs_pp_H1_hh_bbbb_ATLAS13(const StandardModel& SM_i);

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

class Hobs_gg_H1_hh_bbbb_CMS13: public ThObservable {
public:

    Hobs_gg_H1_hh_bbbb_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_hh_gagabb_ATLAS13: public ThObservable {
public:

    Hobs_pp_H1_hh_gagabb_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_hh_gagabb_CMS13: public ThObservable {
public:

    Hobs_pp_H1_hh_gagabb_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_hh_bbtautau_CMS13: public ThObservable {
public:

    Hobs_pp_H1_hh_bbtautau_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_hh_bblnulnu_CMS13: public ThObservable {
public:

    Hobs_pp_H1_hh_bblnulnu_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H1_hh_gagaWW_ATLAS13: public ThObservable {
public:

    Hobs_gg_H1_hh_gagaWW_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_hh_ATLAS8: public ThObservable {
public:

    Hobs_gg_H5_hh_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_hh_bbbb_CMS8: public ThObservable {
public:

    Hobs_pp_H5_hh_bbbb_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_hh_gagabb_CMS8: public ThObservable {
public:

    Hobs_pp_H5_hh_gagabb_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_hh_bbtautau_CMS8: public ThObservable {
public:

    Hobs_gg_H5_hh_bbtautau_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_hh_bbtautau_CMS8: public ThObservable {
public:

    Hobs_pp_H5_hh_bbtautau_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_hh_bbbb_ATLAS13: public ThObservable {
public:

    Hobs_pp_H5_hh_bbbb_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5_hh_bbbb_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to H_5\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_pp_H5_hh_bbbb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5_hh_bbbb_CMS13 constructor.
     */
    Hobs_pp_H5_hh_bbbb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H_5}\cdot BR^{\text{GM}}(H_5\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H_5}\cdot BR(H_5\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_hh_bbbb_CMS13: public ThObservable {
public:

    Hobs_gg_H5_hh_bbbb_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_hh_gagabb_ATLAS13: public ThObservable {
public:

    Hobs_pp_H5_hh_gagabb_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_hh_gagabb_CMS13: public ThObservable {
public:

    Hobs_pp_H5_hh_gagabb_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_hh_bbtautau_CMS13: public ThObservable {
public:

    Hobs_pp_H5_hh_bbtautau_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_hh_bblnulnu_CMS13: public ThObservable {
public:

    Hobs_pp_H5_hh_bblnulnu_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H5_hh_gagaWW_ATLAS13: public ThObservable {
public:

    Hobs_gg_H5_hh_gagaWW_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_hZ_bbZ_ATLAS8: public ThObservable {
public:

    Hobs_gg_H3_hZ_bbZ_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_hZ_bbll_CMS8: public ThObservable {
public:

    Hobs_gg_H3_hZ_bbll_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_hZ_tautauZ_ATLAS8: public ThObservable {
public:

    Hobs_gg_H3_hZ_tautauZ_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_hZ_tautaull_CMS8: public ThObservable {
public:

    Hobs_gg_H3_hZ_tautaull_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_gg_H3_hZ_bbZ_ATLAS13: public ThObservable {
public:

    Hobs_gg_H3_hZ_bbZ_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_bb_H3_hZ_bbZ_ATLAS13: public ThObservable {
public:

    Hobs_bb_H3_hZ_bbZ_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3_H1Z_bbll_CMS8: public ThObservable {
public:

    Hobs_pp_H3_H1Z_bbll_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3_H5Z_bbll_CMS8: public ThObservable {
public:

    Hobs_pp_H3_H5Z_bbll_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H1_H3Z_bbll_CMS8: public ThObservable {
public:

    Hobs_pp_H1_H3Z_bbll_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5_H3Z_bbll_CMS8: public ThObservable {
public:

    Hobs_pp_H5_H3Z_bbll_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3pm_taunu_ATLAS8: public ThObservable {
public:

    Hobs_pp_H3pm_taunu_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3p_taunu_CMS8: public ThObservable {
public:

    Hobs_pp_H3p_taunu_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3pm_taunu_ATLAS13: public ThObservable {
public:

    Hobs_pp_H3pm_taunu_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3pm_taunu_CMS13: public ThObservable {
public:

    Hobs_pp_H3pm_taunu_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3pm_tb_ATLAS8: public ThObservable {
public:

    Hobs_pp_H3pm_tb_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3p_tb_CMS8: public ThObservable {
public:

    Hobs_pp_H3p_tb_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3p_tb1_ATLAS13: public ThObservable {
public:

    Hobs_pp_H3p_tb1_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H3p_tb2_ATLAS13: public ThObservable {
public:

    Hobs_pp_H3p_tb2_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_WZ_H5p_WZ_qqll_ATLAS8: public ThObservable {
public:

    Hobs_WZ_H5p_WZ_qqll_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_WZ_H5p_WZ_lnull_CMS8: public ThObservable {
public:

    Hobs_WZ_H5p_WZ_lnull_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8: public ThObservable {
public:

    Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8: public ThObservable {
public:

    Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8: public ThObservable {
public:

    Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13: public ThObservable {
public:

    Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5ppmm_WW_jjll_CMS8: public ThObservable {
public:

    Hobs_pp_H5ppmm_WW_jjll_CMS8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class Hobs_pp_H5ppmm_WW_jjll_CMS13: public ThObservable {
public:

    Hobs_pp_H5ppmm_WW_jjll_CMS13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

    //log10(sigma*BR)

class log10_tt_H1_tt_TH13: public ThObservable {
public:

    log10_tt_H1_tt_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H1_tt_TH13: public ThObservable {
public:

    log10_bb_H1_tt_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_tt_H3_tt_TH13: public ThObservable {
public:

    log10_tt_H3_tt_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H3_tt_TH13: public ThObservable {
public:

    log10_bb_H3_tt_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H1_bb_TH8: public ThObservable {
public:

    log10_bb_H1_bb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_bb_TH8: public ThObservable {
public:

    log10_gg_H1_bb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_bb_TH13: public ThObservable {
public:

    log10_pp_H1_bb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H1_bb_TH13: public ThObservable {
public:

    log10_bb_H1_bb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H3_bb_TH8: public ThObservable {
public:

    log10_bb_H3_bb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_bb_TH8: public ThObservable {
public:

    log10_gg_H3_bb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3_bb_TH13: public ThObservable {
public:

    log10_pp_H3_bb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H3_bb_TH13: public ThObservable {
public:

    log10_bb_H3_bb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_tautau_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H_1\to \tau\tau@f$ at 8 TeV.
 */
class log10_gg_H1_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_tautau_TH8 constructor.
     */
    log10_gg_H1_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H1_tautau_TH8: public ThObservable {
public:

    log10_bb_H1_tautau_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_tautau_TH13: public ThObservable {
public:

    log10_gg_H1_tautau_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H1_tautau_TH13: public ThObservable {
public:

    log10_bb_H1_tautau_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_tautau_TH8: public ThObservable {
public:

    log10_gg_H3_tautau_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H3_tautau_TH8: public ThObservable {
public:

    log10_bb_H3_tautau_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_tautau_TH13: public ThObservable {
public:

    log10_gg_H3_tautau_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H3_tautau_TH13: public ThObservable {
public:

    log10_bb_H3_tautau_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_gaga_TH8: public ThObservable {
public:

    log10_gg_H1_gaga_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_gaga_TH13: public ThObservable {
public:

    log10_pp_H1_gaga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_gaga_TH13: public ThObservable {
public:

    log10_gg_H1_gaga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_gaga_TH8: public ThObservable {
public:

    log10_gg_H3_gaga_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3_gaga_TH13: public ThObservable {
public:

    log10_pp_H3_gaga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_gaga_TH13: public ThObservable {
public:

    log10_gg_H3_gaga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_gaga_TH8: public ThObservable {
public:

    log10_gg_H5_gaga_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_gaga_TH13: public ThObservable {
public:

    log10_pp_H5_gaga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_gaga_TH13: public ThObservable {
public:

    log10_gg_H5_gaga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_Zga_llga_TH8: public ThObservable {
public:

    log10_pp_H1_Zga_llga_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_Zga_llga_TH13: public ThObservable {
public:

    log10_gg_H1_Zga_llga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_Zga_TH13: public ThObservable {
public:

    log10_gg_H1_Zga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3_Zga_llga_TH8: public ThObservable {
public:

    log10_pp_H3_Zga_llga_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_Zga_llga_TH13: public ThObservable {
public:

    log10_gg_H3_Zga_llga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_Zga_TH13: public ThObservable {
public:

    log10_gg_H3_Zga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_Zga_llga_TH8: public ThObservable {
public:

    log10_pp_H5_Zga_llga_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_Zga_llga_TH13: public ThObservable {
public:

    log10_gg_H5_Zga_llga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_Zga_TH13: public ThObservable {
public:

    log10_gg_H5_Zga_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_ZZ_TH8: public ThObservable {
public:

    log10_gg_H1_ZZ_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H1_ZZ_TH8: public ThObservable {
public:

    log10_VV_H1_ZZ_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_ZZ_llllnunu_TH13: public ThObservable {
public:

    log10_gg_H1_ZZ_llllnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H1_ZZ_llllnunu_TH13: public ThObservable {
public:

    log10_VV_H1_ZZ_llllnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_ZZ_qqllnunu_TH13: public ThObservable {
public:

    log10_gg_H1_ZZ_qqllnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H1_ZZ_qqllnunu_TH13: public ThObservable {
public:

    log10_VV_H1_ZZ_qqllnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_ZZ_llqqnunull_TH13: public ThObservable {
public:

    log10_pp_H1_ZZ_llqqnunull_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H1_ZZ_llqqnunull_TH13: public ThObservable {
public:

    log10_VV_H1_ZZ_llqqnunull_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_ZZ_qqnunu_TH13: public ThObservable {
public:

    log10_pp_H1_ZZ_qqnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_ZZ_TH8: public ThObservable {
public:

    log10_gg_H5_ZZ_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H5_ZZ_TH8: public ThObservable {
public:

    log10_VV_H5_ZZ_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_ZZ_llllnunu_TH13: public ThObservable {
public:

    log10_gg_H5_ZZ_llllnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H5_ZZ_llllnunu_TH13: public ThObservable {
public:

    log10_VV_H5_ZZ_llllnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_ZZ_qqllnunu_TH13: public ThObservable {
public:

    log10_gg_H5_ZZ_qqllnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H5_ZZ_qqllnunu_TH13: public ThObservable {
public:

    log10_VV_H5_ZZ_qqllnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_ZZ_llqqnunull_TH13: public ThObservable {
public:

    log10_pp_H5_ZZ_llqqnunull_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H5_ZZ_llqqnunull_TH13: public ThObservable {
public:

    log10_VV_H5_ZZ_llqqnunull_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_ZZ_qqnunu_TH13: public ThObservable {
public:

    log10_pp_H5_ZZ_qqnunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_WW_TH8: public ThObservable {
public:

    log10_gg_H1_WW_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H1_WW_TH8: public ThObservable {
public:

    log10_VV_H1_WW_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_WW_enumunu_TH13: public ThObservable {
public:

    log10_gg_H1_WW_enumunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H1_WW_enumunu_TH13: public ThObservable {
public:

    log10_VV_H1_WW_enumunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_WW_lnuqq_TH13: public ThObservable {
public:

    log10_gg_H1_WW_lnuqq_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H1_WW_lnuqq_TH13: public ThObservable {
public:

    log10_VV_H1_WW_lnuqq_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_ggVV_H1_WW_lnulnu_TH13: public ThObservable {
public:

    log10_ggVV_H1_WW_lnulnu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_WW_lnuqq_TH13: public ThObservable {
public:

    log10_pp_H1_WW_lnuqq_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_WW_TH8: public ThObservable {
public:

    log10_gg_H5_WW_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H5_WW_TH8: public ThObservable {
public:

    log10_VV_H5_WW_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_WW_enumunu_TH13: public ThObservable {
public:

    log10_gg_H5_WW_enumunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H5_WW_enumunu_TH13: public ThObservable {
public:

    log10_VV_H5_WW_enumunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_WW_lnuqq_TH13: public ThObservable {
public:

    log10_gg_H5_WW_lnuqq_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_VV_H5_WW_lnuqq_TH13: public ThObservable {
public:

    log10_VV_H5_WW_lnuqq_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_ggVV_H5_WW_lnulnu_TH13: public ThObservable {
public:

    log10_ggVV_H5_WW_lnulnu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_WW_lnuqq_TH13: public ThObservable {
public:

    log10_pp_H5_WW_lnuqq_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_VV_TH8: public ThObservable {
public:

    log10_pp_H1_VV_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_mu_pp_H1_VV_TH8: public ThObservable {
public:

    log10_mu_pp_H1_VV_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_VV_TH13: public ThObservable {
public:

    log10_pp_H1_VV_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_VV_TH8: public ThObservable {
public:

    log10_pp_H5_VV_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_mu_pp_H5_VV_TH8: public ThObservable {
public:

    log10_mu_pp_H5_VV_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_VV_TH13: public ThObservable {
public:

    log10_pp_H5_VV_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_hh_TH8: public ThObservable {
public:

    log10_gg_H1_hh_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_hh_bbbb_TH8: public ThObservable {
public:

    log10_pp_H1_hh_bbbb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_hh_gagabb_TH8: public ThObservable {
public:

    log10_pp_H1_hh_gagabb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_hh_bbtautau_TH8: public ThObservable {
public:

    log10_gg_H1_hh_bbtautau_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_hh_bbtautau_TH8: public ThObservable {
public:

    log10_pp_H1_hh_bbtautau_TH8(const StandardModel& SM_i);

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

class log10_gg_H1_hh_bbbb_TH13: public ThObservable {
public:

    log10_gg_H1_hh_bbbb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_hh_gagabb_TH13: public ThObservable {
public:

    log10_pp_H1_hh_gagabb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_hh_bbtautau_TH13: public ThObservable {
public:

    log10_pp_H1_hh_bbtautau_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_hh_bblnulnu_TH13: public ThObservable {
public:

    log10_pp_H1_hh_bblnulnu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H1_hh_gagaWW_TH13: public ThObservable {
public:

    log10_gg_H1_hh_gagaWW_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_hh_TH8: public ThObservable {
public:

    log10_gg_H5_hh_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_hh_bbbb_TH8: public ThObservable {
public:

    log10_pp_H5_hh_bbbb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_hh_gagabb_TH8: public ThObservable {
public:

    log10_pp_H5_hh_gagabb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_hh_bbtautau_TH8: public ThObservable {
public:

    log10_gg_H5_hh_bbtautau_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_hh_bbtautau_TH8: public ThObservable {
public:

    log10_pp_H5_hh_bbtautau_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_hh_bbbb_TH13: public ThObservable {
public:

    log10_pp_H5_hh_bbbb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_hh_bbbb_TH13: public ThObservable {
public:

    log10_gg_H5_hh_bbbb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_hh_gagabb_TH13: public ThObservable {
public:

    log10_pp_H5_hh_gagabb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_hh_bbtautau_TH13: public ThObservable {
public:

    log10_pp_H5_hh_bbtautau_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_hh_bblnulnu_TH13: public ThObservable {
public:

    log10_pp_H5_hh_bblnulnu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H5_hh_gagaWW_TH13: public ThObservable {
public:

    log10_gg_H5_hh_gagaWW_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_hZ_bbZ_TH8: public ThObservable {
public:

    log10_gg_H3_hZ_bbZ_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_hZ_bbll_TH8: public ThObservable {
public:

    log10_gg_H3_hZ_bbll_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_hZ_tautauZ_TH8: public ThObservable {
public:

    log10_gg_H3_hZ_tautauZ_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_hZ_tautaull_TH8: public ThObservable {
public:

    log10_gg_H3_hZ_tautaull_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_gg_H3_hZ_bbZ_TH13: public ThObservable {
public:

    log10_gg_H3_hZ_bbZ_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_bb_H3_hZ_bbZ_TH13: public ThObservable {
public:

    log10_bb_H3_hZ_bbZ_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3_H1Z_bbll_TH8: public ThObservable {
public:

    log10_pp_H3_H1Z_bbll_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3_H5Z_bbll_TH8: public ThObservable {
public:

    log10_pp_H3_H5Z_bbll_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H1_H3Z_bbll_TH8: public ThObservable {
public:

    log10_pp_H1_H3Z_bbll_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5_H3Z_bbll_TH8: public ThObservable {
public:

    log10_pp_H5_H3Z_bbll_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3pm_taunu_TH8: public ThObservable {
public:

    log10_pp_H3pm_taunu_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3p_taunu_TH8: public ThObservable {
public:

    log10_pp_H3p_taunu_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3pm_taunu_TH13: public ThObservable {
public:

    log10_pp_H3pm_taunu_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3pm_tb_TH8: public ThObservable {
public:

    log10_pp_H3pm_tb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3p_tb_TH8: public ThObservable {
public:

    log10_pp_H3p_tb_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H3p_tb_TH13: public ThObservable {
public:

    log10_pp_H3p_tb_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_WZ_H5pm_WZ_TH8: public ThObservable {
public:

    log10_WZ_H5pm_WZ_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_WZ_H5pm_WZ_TH13: public ThObservable {
public:

    log10_WZ_H5pm_WZ_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5ppmmH5mmpp_TH8: public ThObservable {
public:

    log10_pp_H5ppmmH5mmpp_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5ppmmH5mmpp_TH13: public ThObservable {
public:

    log10_pp_H5ppmmH5mmpp_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5ppmm_WW_TH8: public ThObservable {
public:

    log10_pp_H5ppmm_WW_TH8(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

class log10_pp_H5ppmm_WW_TH13: public ThObservable {
public:

    log10_pp_H5ppmm_WW_TH13(const StandardModel& SM_i);

    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

#endif	/* GMDIRECTSEARCHES_H */
