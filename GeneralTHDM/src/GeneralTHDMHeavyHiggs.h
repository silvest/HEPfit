/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMHEAVYHIGGS_H
#define GENERALTHDMHEAVYHIGGS_H


#include "ThObservable.h"

class GeneralTHDM;
class GeneralTHDMcache;

/**
 * @class GeneralTHDMHeavyHiggs
 * @ingroup GeneralTHDM 
 * @brief Base class for direct heavy Higgs search observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The @f$\gamma \gamma@f$, @f$Z\gamma@f$ and @f$gg@f$ decay widths are calculated at one-loop
 * following @cite Gunion:1989we and @cite Aoki:2009ha.
 */


/**
 * @class Hobs_tt_phi2_tt_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$t\bar t \to phi2\to t\bar t@f$.
 */
class Hobs_tt_phi2_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_tt_phi2_tt_ATLAS13 constructor.
     */
    Hobs_tt_phi2_tt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{t\bar t\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to t\bar t)]_{\text{theo}} / [\sigma_{t\bar t\to phi2}\cdot BR(phi2\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi2_tt_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b \to phi2\to t\bar t@f$.
 */
class Hobs_bb_phi2_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi2_tt_ATLAS13 constructor.
     */
    Hobs_bb_phi2_tt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to t\bar t)]_{\text{theo}} / [\sigma_{b\bar b\to phi2}\cdot BR(phi2\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi2_bb_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b \to phi2\to b\bar b@f$.
 */
class Hobs_bb_phi2_bb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi2_bb_CMS8 constructor.
     */
    Hobs_bb_phi2_bb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to phi2}\cdot BR(phi2\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_bb_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b \to phi2\to b\bar b@f$.
 */
class Hobs_gg_phi2_bb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi2_bb_CMS8 constructor.
     */
    Hobs_gg_phi2_bb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to phi2}\cdot BR(phi2\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_bb_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to b\bar b@f$.
 */
class Hobs_pp_phi2_bb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_bb_CMS13 constructor.
     */
    Hobs_pp_phi2_bb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi2_bb_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to b\bar b@f$.
 */
class Hobs_bb_phi2_bb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_bb_CMS13 constructor.
     */
    Hobs_bb_phi2_bb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_tautau_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to \tau\tau@f$.
 */
class Hobs_gg_phi2_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_tautau_ATLAS8 constructor.
     */
    Hobs_gg_phi2_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi2_tautau_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to phi2\to \tau\tau@f$.
 */
class Hobs_bb_phi2_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi2_tautau_ATLAS8 constructor.
     */
    Hobs_bb_phi2_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to phi2}\cdot BR(phi2\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};




/**
 * @class Hobs_gg_phi2_tautau_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to \tau\tau@f$.
 */
class Hobs_gg_phi2_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_tautau_CMS8 constructor.
     */
    Hobs_gg_phi2_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
        const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi2_tautau_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to phi2\to \tau\tau@f$.
 */
class Hobs_bb_phi2_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi2_tautau_CMS8 constructor.
     */
    Hobs_bb_phi2_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to phi2}\cdot BR(phi2\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};




/**
 * @class Hobs_gg_phi2_tautau_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to \tau \tau@f$.
 */
class Hobs_gg_phi2_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_tautau_ATLAS13 constructor.
     */
    Hobs_gg_phi2_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau \tau)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi2_tautau_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to phi2\to \tau \tau@f$.
 */
class Hobs_bb_phi2_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi2_tautau_ATLAS13 constructor.
     */
    Hobs_bb_phi2_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau \tau)]_{\text{theo}} / [\sigma_{b\bar b\to phi2}\cdot BR(phi2\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_tautau_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to \tau \tau@f$.
 */
class Hobs_gg_phi2_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_tautau_CMS13 constructor.
     */
    Hobs_gg_phi2_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau \tau)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi2_tautau_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to phi2\to \tau \tau@f$.
 */
class Hobs_bb_phi2_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi2_tautau_CMS13 constructor.
     */
    Hobs_bb_phi2_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau \tau)]_{\text{theo}} / [\sigma_{b\bar b\to phi2}\cdot BR(phi2\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_gaga_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to \gamma \gamma@f$.
 */
class Hobs_gg_phi2_gaga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_gaga_ATLAS8 constructor.
     */
    Hobs_gg_phi2_gaga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \gamma \gamma)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_gaga_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to \gamma \gamma@f$.
 */
class Hobs_pp_phi2_gaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_gaga_ATLAS13 constructor.
     */
    Hobs_pp_phi2_gaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \gamma \gamma)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi2_gaga_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to \gamma \gamma@f$.
 */
class Hobs_gg_phi2_gaga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_gaga_CMS13 constructor.
     */
    Hobs_gg_phi2_gaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \gamma \gamma)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_Zga_llga_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to Z\gamma@f$.
 */
class Hobs_pp_phi2_Zga_llga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_Zga_llga_ATLAS13 constructor.
     */
    Hobs_pp_phi2_Zga_llga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to Z\gamma)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to Z\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_Zga_llga_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to Z\gamma@f$.
 */
class Hobs_pp_phi2_Zga_llga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_Zga_llga_CMS8 constructor.
     */
    Hobs_pp_phi2_Zga_llga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to Z\gamma)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to Z\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_Zga_llga_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to Z\gamma@f$.
 */
class Hobs_gg_phi2_Zga_llga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_Zga_llga_ATLAS13 constructor.
     */
    Hobs_gg_phi2_Zga_llga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to Z\gamma)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to Z\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi2_Zga_qqga_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to Z\gamma@f$.
 */
class Hobs_gg_phi2_Zga_qqga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_Zga_qqga_ATLAS13 constructor.
     */
    Hobs_gg_phi2_Zga_qqga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to Z\gamma)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to Z\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_Zga_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to Z\gamma@f$.
 */
class Hobs_gg_phi2_Zga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_Zga_CMS13 constructor.
     */
    Hobs_gg_phi2_Zga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to Z\gamma)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to Z\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class Hobs_gg_phi2_ZZ_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to phi2\to ZZ@f$.
 */
class Hobs_gg_phi2_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_ZZ_ATLAS8 constructor.
     */
    Hobs_gg_phi2_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_VV_phi2_ZZ_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV \to phi2\to ZZ@f$.
 */
class Hobs_VV_phi2_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi2_ZZ_ATLAS8 constructor.
     */
    Hobs_VV_phi2_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{VV\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]_{\text{theo}} / [\sigma_{VV\to phi2}\cdot BR(phi2\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_ZZ_llllnunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to ZZ@f$.
 */
class Hobs_gg_phi2_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_ZZ_llllnunu_ATLAS13 constructor.
     */
    Hobs_gg_phi2_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi2_ZZ_llllnunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV\to phi2\to ZZ@f$.
 */
class Hobs_VV_phi2_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi2_ZZ_llllnunu_ATLAS13 constructor.
     */
    Hobs_VV_phi2_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{VV\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]_{\text{theo}} / [\sigma_{VV\to phi2}\cdot BR(phi2\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi2_ZZ_llnunuqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to ZZ@f$.
 */
class Hobs_gg_phi2_ZZ_llnunuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_ZZ_llnunuqq_ATLAS13 constructor.
     */
    Hobs_gg_phi2_ZZ_llnunuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi2_ZZ_llnunuqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to ZZ@f$.
 */
class Hobs_VV_phi2_ZZ_llnunuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi2_ZZ_llnunuqq_ATLAS13 constructor.
     */
    Hobs_VV_phi2_ZZ_llnunuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi2_ZZ_llqqnunull_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to ZZ@f$.
 */
class Hobs_pp_phi2_ZZ_llqqnunull_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_ZZ_llqqnunull_CMS13 constructor.
     */
    Hobs_pp_phi2_ZZ_llqqnunull_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to ZZ)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_ZZ_qqnunu_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to ZZ@f$.
 */
class Hobs_pp_phi2_ZZ_qqnunu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_ZZ_qqnunu_CMS13 constructor.
     */
    Hobs_pp_phi2_ZZ_qqnunu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to ZZ)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_WW_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to phi2\to WW@f$.
 */
class Hobs_gg_phi2_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_WW_ATLAS8 constructor.
     */
    Hobs_gg_phi2_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to WW)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi2_WW_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV \to phi2\to WW@f$.
 */
class Hobs_VV_phi2_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi2_WW_ATLAS8 constructor.
     */
    Hobs_VV_phi2_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{VV\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to WW)]_{\text{theo}} / [\sigma_{VV\to phi2}\cdot BR(phi2\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi2_WW_enumunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_gg_phi2_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_WW_enumunu_ATLAS13 constructor.
     */
    Hobs_gg_phi2_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi2_WW_enumunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_VV_phi2_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi2_WW_enumunu_ATLAS13 constructor.
     */
    Hobs_VV_phi2_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_ggVV_phi2_WW_lnulnu_CMS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_ggVV_phi2_WW_lnulnu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggVV_phi2_WW_lnulnu_CMS13 constructor.
     */
    Hobs_ggVV_phi2_WW_lnulnu_CMS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_WW_lnuqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_gg_phi2_WW_lnuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_WW_lnuqq_ATLAS13 constructor.
     */
    Hobs_gg_phi2_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi2_WW_lnuqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_VV_phi2_WW_lnuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi2_WW_lnuqq_ATLAS13 constructor.
     */
    Hobs_VV_phi2_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_WW_lnuqq_CMS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_pp_phi2_WW_lnuqq_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_WW_lnuqq_CMS13 constructor.
     */
    Hobs_pp_phi2_WW_lnuqq_CMS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_mu_pp_phi2_VV_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the signal strength of the process @f$pp \to phi2\to VV@f$.
 */
class Hobs_mu_pp_phi2_VV_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_mu_pp_phi2_VV_CMS8 constructor.
     */
    Hobs_mu_pp_phi2_VV_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\mu_H^{\text{GTHDM}}(phi2\to VV)]_{\text{theo}} / [\mu_H(phi2\to VV)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi2_VV_qqqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to WW,ZZ@f$.
 */
class Hobs_pp_phi2_VV_qqqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_VV_qqqq_ATLAS13 constructor.
     */
    Hobs_pp_phi2_VV_qqqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot [BR^{\text{GTHDM}}(phi2\to WW)+BR^{\text{GTHDM}}(phi2\to ZZ)]]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot [BR^{\text{GTHDM}}(phi2\to WW)+BR^{\text{GTHDM}}(phi2\to ZZ)]]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi2_phi1phi1_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to phi2\to phi1phi1@f$.
 */
class Hobs_gg_phi2_phi1phi1_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_phi1phi1_ATLAS8 constructor.
     */
    Hobs_gg_phi2_phi1phi1_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to phi1phi1)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi1phi1_bbbb_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to phi2\to phi1phi1\to b\bar b b\bar b@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbbb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbbb_CMS8 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbbb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi1phi1\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi2_phi1phi1_bbgaga_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to phi2\to phi1phi1\to \gamma\gamma b\bar b@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbgaga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbgaga_CMS8 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbgaga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to \gamma\gamma b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi1phi1\to \gamma\gamma b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_phi1phi1_bbtautau_CMS8
 * @ingroup GeneralTHDM
 * @brief 
 * 
 * **/

class Hobs_gg_phi2_phi1phi1_bbtautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_phi1phi1_bbtautau_CMS8 constructor.
     */
    Hobs_gg_phi2_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi1phi1_bbtautau_CMS8
 * @ingroup GeneralTHDM
 * @brief 
 * 
 * **/

class Hobs_pp_phi2_phi1phi1_bbtautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbtautau_CMS8 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b b\bar b@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbbb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi1phi1\to b\bar b b\bar b)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi1phi1_bbbb_1_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b b\bar b@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbbb_1_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbbb_1_CMS13 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbbb_1_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi1phi1\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi2_phi1phi1_bbbb_2_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi2\to phi1phi1\to b\bar b b\bar b@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbbb_2_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbbb_2_CMS13 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbbb_2_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{gg\to phi2}\cdot BR(phi2\to phi1phi1\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbgaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi1phi1)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi1phi1_bbgaga_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to bb \gamma \gamma@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbgaga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbgaga_CMS13 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbgaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to bb \gamma \gamma)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi1phi1\to bb \gamma \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi1phi1_bbtautau_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbtautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbtautau_ATLAS13 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbtautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b \tau \tau)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi1phi1\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi2_phi1phi1_bbtautau_1_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbtautau_1_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_pp_phi2_phi1phi1_bbtautau_1_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi1phi1_bbtautau_2_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbtautau_2_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_pp_phi2_phi1phi1_bbtautau_2_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi2_phi1phi1_bbVV_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b VV\to b\bar b \ell \ell \nu \nu@f$.
 */
class Hobs_pp_phi2_phi1phi1_bbVV_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbVV_CMS13 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbVV_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b VV\to b\bar b \ell \ell \nu \nu)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi1phi1\to b\bar b VV\to b\bar b \ell \ell \nu \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi1phi1_bbWW_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1 [\to b\bar b WW]\f$.
 */
class Hobs_pp_phi2_phi1phi1_bbWW_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi1phi1_bbWW_ATLAS13 constructor.
     */
    Hobs_pp_phi2_phi1phi1_bbWW_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b WW)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi1phi1\to b\bar b WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_phi1phi1_gagaWW_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi2_phi1phi1_gagaWW_ATLAS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi2_phi1phi1_gagaWW_ATLAS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_phi1Z_bbZ_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi2_phi1Z_bbZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi2_phi1Z_bbZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_phi1Z_bbll_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi2_phi1Z_bbll_CMS8: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi2_phi1Z_bbll_CMS8(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_phi1Z_tautauZ_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi2_phi1Z_tautauZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi2_phi1Z_tautauZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};
/**
 * @class Hobs_gg_phi2_phi1Z_tautaull_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi2_phi1Z_tautaull_CMS8: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi2_phi1Z_tautaull_CMS8(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi2_phi1Z_bbZ_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi2_phi1Z_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi2_phi1Z_bbZ_ATLAS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi2_phi1Z_bbZ_1_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi2_phi1Z_bbZ_1_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi2_phi1Z_bbZ_1_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_phi1Z_bbZ_2_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi2_phi1Z_bbZ_2_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi2_phi1Z_bbZ_2_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi2_phi1Z_bbZ_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_bb_phi2_phi1Z_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_bb_phi2_phi1Z_bbZ_ATLAS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi2_phi1Z_bbZ_1_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_bb_phi2_phi1Z_bbZ_1_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_bb_phi2_phi1Z_bbZ_1_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi2_phi1Z_bbZ_2_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_bb_phi2_phi1Z_bbZ_2_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_bb_phi2_phi1Z_bbZ_2_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_tt_phi3_tt_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$t\bar t \to phi3\to t\bar t@f$.
 */
class Hobs_tt_phi3_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_tt_phi3_tt_ATLAS13 constructor.
     */
    Hobs_tt_phi3_tt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{t\bar t\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to t\bar t)]_{\text{theo}} / [\sigma_{t\bar t\to phi3}\cdot BR(phi3\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi3_tt_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b \to phi3\to t\bar t@f$.
 */
class Hobs_bb_phi3_tt_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi3_tt_ATLAS13 constructor.
     */
    Hobs_bb_phi3_tt_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to t\bar t)]_{\text{theo}} / [\sigma_{b\bar b\to phi3}\cdot BR(phi3\to t\bar t)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi3_bb_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b \to phi3\to b\bar b@f$.
 */
class Hobs_bb_phi3_bb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi3_bb_CMS8 constructor.
     */
    Hobs_bb_phi3_bb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to phi3}\cdot BR(phi3\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_bb_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b \to phi3\to b\bar b@f$.
 */
class Hobs_gg_phi3_bb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi3_bb_CMS8 constructor.
     */
    Hobs_gg_phi3_bb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to phi3}\cdot BR(phi3\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_bb_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to b\bar b@f$.
 */
class Hobs_pp_phi3_bb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_bb_CMS13 constructor.
     */
    Hobs_pp_phi3_bb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi3_bb_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to b\bar b@f$.
 */
class Hobs_bb_phi3_bb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_bb_CMS13 constructor.
     */
    Hobs_bb_phi3_bb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_tautau_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to \tau\tau@f$.
 */
class Hobs_gg_phi3_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_tautau_ATLAS8 constructor.
     */
    Hobs_gg_phi3_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi3_tautau_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to phi3\to \tau\tau@f$.
 */
class Hobs_bb_phi3_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi3_tautau_ATLAS8 constructor.
     */
    Hobs_bb_phi3_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to phi3}\cdot BR(phi3\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};




/**
 * @class Hobs_gg_phi3_tautau_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to \tau\tau@f$.
 */
class Hobs_gg_phi3_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_tautau_CMS8 constructor.
     */
    Hobs_gg_phi3_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
        const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi3_tautau_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to phi3\to \tau\tau@f$.
 */
class Hobs_bb_phi3_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi3_tautau_CMS8 constructor.
     */
    Hobs_bb_phi3_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to phi3}\cdot BR(phi3\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};




/**
 * @class Hobs_gg_phi3_tautau_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to \tau \tau@f$.
 */
class Hobs_gg_phi3_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_tautau_ATLAS13 constructor.
     */
    Hobs_gg_phi3_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau \tau)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_bb_phi3_tautau_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$b\bar b\to phi3\to \tau \tau@f$.
 */
class Hobs_bb_phi3_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi3_tautau_ATLAS13 constructor.
     */
    Hobs_bb_phi3_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau \tau)]_{\text{theo}} / [\sigma_{b\bar b\to phi3}\cdot BR(phi3\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_tautau_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to \tau \tau@f$.
 */
class Hobs_gg_phi3_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_tautau_CMS13 constructor.
     */
    Hobs_gg_phi3_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau \tau)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi3_tautau_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$b\bar b\to phi3\to \tau \tau@f$.
 */
class Hobs_bb_phi3_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi3_tautau_CMS13 constructor.
     */
    Hobs_bb_phi3_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau \tau)]_{\text{theo}} / [\sigma_{b\bar b\to phi3}\cdot BR(phi3\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_gaga_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to \gamma \gamma@f$.
 */
class Hobs_pp_phi3_gaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_gaga_ATLAS13 constructor.
     */
    Hobs_pp_phi3_gaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \gamma \gamma)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_gaga_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to \gamma \gamma@f$.
 */
class Hobs_gg_phi3_gaga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_gaga_ATLAS8 constructor.
     */
    Hobs_gg_phi3_gaga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \gamma \gamma)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_gaga_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to \gamma \gamma@f$.
 */
class Hobs_gg_phi3_gaga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_gaga_CMS13 constructor.
     */
    Hobs_gg_phi3_gaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \gamma \gamma)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_Zga_llga_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to Z\gamma@f$.
 */
class Hobs_pp_phi3_Zga_llga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_Zga_llga_ATLAS13 constructor.
     */
    Hobs_pp_phi3_Zga_llga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to Z\gamma)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to Z\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_Zga_llga_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to Z\gamma@f$.
 */
class Hobs_pp_phi3_Zga_llga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_Zga_llga_CMS8 constructor.
     */
    Hobs_pp_phi3_Zga_llga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to Z\gamma)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to Z\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_Zga_llga_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to Z\gamma@f$.
 */
class Hobs_gg_phi3_Zga_llga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_Zga_llga_ATLAS13 constructor.
     */
    Hobs_gg_phi3_Zga_llga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to Z\gamma)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to Z\gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi3_Zga_qqga_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to Z\gamma@f$.
 */
class Hobs_gg_phi3_Zga_qqga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_Zga_qqga_ATLAS13 constructor.
     */
    Hobs_gg_phi3_Zga_qqga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to Z\gamma)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to Z\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_Zga_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to Z\gamma@f$.
 */
class Hobs_gg_phi3_Zga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_Zga_CMS13 constructor.
     */
    Hobs_gg_phi3_Zga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to Z\gamma)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to Z\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_ZZ_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to phi3\to ZZ@f$.
 */
class Hobs_gg_phi3_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_ZZ_ATLAS8 constructor.
     */
    Hobs_gg_phi3_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_VV_phi3_ZZ_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV \to phi3\to ZZ@f$.
 */
class Hobs_VV_phi3_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi3_ZZ_ATLAS8 constructor.
     */
    Hobs_VV_phi3_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{VV\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]_{\text{theo}} / [\sigma_{VV\to phi3}\cdot BR(phi3\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_ZZ_llllnunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to ZZ@f$.
 */
class Hobs_gg_phi3_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_ZZ_llllnunu_ATLAS13 constructor.
     */
    Hobs_gg_phi3_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi3_ZZ_llllnunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV\to phi3\to ZZ@f$.
 */
class Hobs_VV_phi3_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi3_ZZ_llllnunu_ATLAS13 constructor.
     */
    Hobs_VV_phi3_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{VV\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]_{\text{theo}} / [\sigma_{VV\to phi3}\cdot BR(phi3\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi3_ZZ_llnunuqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to ZZ@f$.
 */
class Hobs_gg_phi3_ZZ_llnunuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_ZZ_llnunuqq_ATLAS13 constructor.
     */
    Hobs_gg_phi3_ZZ_llnunuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi3_ZZ_llnunuqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to ZZ@f$.
 */
class Hobs_VV_phi3_ZZ_llnunuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi3_ZZ_llnunuqq_ATLAS13 constructor.
     */
    Hobs_VV_phi3_ZZ_llnunuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to ZZ)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi3_ZZ_llqqnunull_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to ZZ@f$.
 */
class Hobs_pp_phi3_ZZ_llqqnunull_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_ZZ_llqqnunull_CMS13 constructor.
     */
    Hobs_pp_phi3_ZZ_llqqnunull_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to ZZ)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_ZZ_qqnunu_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to ZZ@f$.
 */
class Hobs_pp_phi3_ZZ_qqnunu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_ZZ_qqnunu_CMS13 constructor.
     */
    Hobs_pp_phi3_ZZ_qqnunu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to ZZ)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_WW_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to phi3\to WW@f$.
 */
class Hobs_gg_phi3_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_WW_ATLAS8 constructor.
     */
    Hobs_gg_phi3_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to WW)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi3_WW_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$VV \to phi3\to WW@f$.
 */
class Hobs_VV_phi3_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi3_WW_ATLAS8 constructor.
     */
    Hobs_VV_phi3_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{VV\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to WW)]_{\text{theo}} / [\sigma_{VV\to phi3}\cdot BR(phi3\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi3_WW_enumunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_gg_phi3_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_WW_enumunu_ATLAS13 constructor.
     */
    Hobs_gg_phi3_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi3_WW_enumunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_VV_phi3_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi3_WW_enumunu_ATLAS13 constructor.
     */
    Hobs_VV_phi3_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_ggVV_phi3_WW_lnulnu_CMS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_ggVV_phi3_WW_lnulnu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggVV_phi3_WW_lnulnu_CMS13 constructor.
     */
    Hobs_ggVV_phi3_WW_lnulnu_CMS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_WW_lnuqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_gg_phi3_WW_lnuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_WW_lnuqq_ATLAS13 constructor.
     */
    Hobs_gg_phi3_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_VV_phi3_WW_lnuqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_VV_phi3_WW_lnuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_phi3_WW_lnuqq_ATLAS13 constructor.
     */
    Hobs_VV_phi3_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_WW_lnuqq_CMS13
 * @ingroup GeneralTHDM
 * @brief 
 */
class Hobs_pp_phi3_WW_lnuqq_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_WW_lnuqq_CMS13 constructor.
     */
    Hobs_pp_phi3_WW_lnuqq_CMS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_mu_pp_phi3_VV_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the signal strength of the process @f$pp \to phi3\to VV@f$.
 */
class Hobs_mu_pp_phi3_VV_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_mu_pp_phi3_VV_CMS8 constructor.
     */
    Hobs_mu_pp_phi3_VV_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\mu_H^{\text{GTHDM}}(phi3\to VV)]_{\text{theo}} / [\mu_H(phi3\to VV)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi3_VV_qqqq_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to WW,ZZ@f$.
 */
class Hobs_pp_phi3_VV_qqqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_VV_qqqq_ATLAS13 constructor.
     */
    Hobs_pp_phi3_VV_qqqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot [BR^{\text{GTHDM}}(phi3\to WW)+BR^{\text{GTHDM}}(phi3\to ZZ)]]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot [BR^{\text{GTHDM}}(phi3\to WW)+BR^{\text{GTHDM}}(phi3\to ZZ)]]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi3_phi1phi1_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$gg \to phi3\to phi1phi1@f$.
 */
class Hobs_gg_phi3_phi1phi1_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_phi1phi1_ATLAS8 constructor.
     */
    Hobs_gg_phi3_phi1phi1_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to phi1phi1)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi1phi1_bbbb_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to phi3\to phi1phi1\to b\bar b b\bar b@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbbb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbbb_CMS8 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbbb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi1phi1_bbgaga_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to phi3\to phi1phi1\to \gamma\gamma b\bar b@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbgaga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbgaga_CMS8 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbgaga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to \gamma\gamma b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1\to \gamma\gamma b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_phi1phi1_bbtautau_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio
 */
class Hobs_gg_phi3_phi1phi1_bbtautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_phi1phi1_bbtautau_CMS8 constructor.
     */
    Hobs_gg_phi3_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi1phi1_bbtautau_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp \to phi3\to phi1phi1\to \tau\tau b\bar b@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbtautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbtautau_CMS8 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbtautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to \tau\tau b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1\to \tau\tau b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b b\bar b@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbbb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1\to b\bar b b\bar b)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi1phi1_bbbb_1_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b b\bar b@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbbb_1_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbbb_1_CMS13 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbbb_1_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi3_phi1phi1_bbbb_2_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$gg\to phi3\to phi1phi1\to b\bar b b\bar b@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbbb_2_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbbb_2_CMS13 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbbb_2_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{gg\to phi3}\cdot BR(phi3\to phi1phi1\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbgaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi1phi1_bbgaga_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to bb \gamma \gamma@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbgaga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbgaga_CMS13 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbgaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to bb \gamma \gamma)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1\to bb \gamma \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi1phi1_bbtautau_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbtautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbtautau_ATLAS13 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbtautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b \tau \tau)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1\to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi3_phi1phi1_bbtautau_1_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbtautau_1_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_pp_phi3_phi1phi1_bbtautau_1_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi1phi1_bbtautau_2_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbtautau_2_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_pp_phi3_phi1phi1_bbtautau_2_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi3_phi1phi1_bbVV_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b VV\to b\bar b \ell \ell \nu \nu@f$.
 */
class Hobs_pp_phi3_phi1phi1_bbVV_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbVV_CMS13 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbVV_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b VV\to b\bar b \ell \ell \nu \nu)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1\to b\bar b VV\to b\bar b \ell \ell \nu \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_phi3_phi1phi1_bbWW_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1 [\to b\bar b WW]\f$.
 */
class Hobs_pp_phi3_phi1phi1_bbWW_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi1phi1_bbWW_ATLAS13 constructor.
     */
    Hobs_pp_phi3_phi1phi1_bbWW_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b WW)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi1phi1\to b\bar b WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class Hobs_gg_phi3_phi1phi1_gagaWW_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi3_phi1phi1_gagaWW_ATLAS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi3_phi1phi1_gagaWW_ATLAS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_phi1Z_bbZ_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi3_phi1Z_bbZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi3_phi1Z_bbZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_phi1Z_bbll_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi3_phi1Z_bbll_CMS8: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi3_phi1Z_bbll_CMS8(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_phi1Z_tautauZ_ATLAS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi3_phi1Z_tautauZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi3_phi1Z_tautauZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};
/**
 * @class Hobs_gg_phi3_phi1Z_tautaull_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi3_phi1Z_tautaull_CMS8: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi3_phi1Z_tautaull_CMS8(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_gg_phi3_phi1Z_bbZ_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi3_phi1Z_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi3_phi1Z_bbZ_ATLAS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_phi1Z_bbZ_1_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi3_phi1Z_bbZ_1_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi3_phi1Z_bbZ_1_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_phi1Z_bbZ_2_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_gg_phi3_phi1Z_bbZ_2_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_gg_phi3_phi1Z_bbZ_2_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi3_phi1Z_bbZ_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_bb_phi3_phi1Z_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_bb_phi3_phi1Z_bbZ_ATLAS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi3_phi1Z_bbZ_1_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_bb_phi3_phi1Z_bbZ_1_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_bb_phi3_phi1Z_bbZ_1_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi3_phi1Z_bbZ_2_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$.
 */
class Hobs_bb_phi3_phi1Z_bbZ_2_CMS13: public ThObservable {
public:

    /**
     * @brief
     */
    Hobs_bb_phi3_phi1Z_bbZ_2_CMS13(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi2Z_bbll_1_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to b\bar b \ell\ell@f$.
 */
class Hobs_pp_phi3_phi2Z_bbll_1_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi2Z_bbll_1_CMS8 constructor.
     */
    Hobs_pp_phi3_phi2Z_bbll_1_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi2Z\to b\bar b\ell \ell)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to phi2Z\to b\bar b\ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi3Z_bbll_1_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to b\bar b \ell\ell@f$.
 */
class Hobs_pp_phi2_phi3Z_bbll_1_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi3Z_bbll_1_CMS8 constructor.
     */
    Hobs_pp_phi2_phi3Z_bbll_1_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi2Z\to b\bar b\ell \ell)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi3\to phi3Z\to b\bar b\ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi2Z_bbll_2_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to b\bar b \ell\ell@f$.
 */
class Hobs_pp_phi3_phi2Z_bbll_2_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi2Z_bbll_2_CMS8 constructor.
     */
    Hobs_pp_phi3_phi2Z_bbll_2_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi2Z\to b\bar b\ell \ell)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi3Z\to b\bar b\ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi3Z_bbll_2_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to b\bar b \ell\ell@f$.
 */
class Hobs_pp_phi2_phi3Z_bbll_2_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi3Z_bbll_2_CMS8 constructor.
     */
    Hobs_pp_phi2_phi3Z_bbll_2_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi3Z\to b\bar b\ell \ell)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi3Z\to b\bar b\ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi2Z_tautaull_1_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi2Z\to \tau\tau\ell\ell@f$.
 */
class Hobs_pp_phi3_phi2Z_tautaull_1_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi2Z_tautaull_1_CMS8 constructor.
     */
    Hobs_pp_phi3_phi2Z_tautaull_1_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi2Z\to \tau\tau\ell \ell)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to AZ\to \tau\tau\ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi3Z_tautaull_1_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi2\to phi3Z\to \tau\tau\ell\ell@f$.
 */
class Hobs_pp_phi2_phi3Z_tautaull_1_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi3Z_tautaull_1_CMS8 constructor.
     */
    Hobs_pp_phi2_phi3Z_tautaull_1_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi3Z\to \tau\tau\ell \ell)]_{\text{theo}} / [\sigma_{pp\to phi2}\cdot BR(phi2\to phi3Z\to \tau\tau\ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi3_phi2Z_tautaull_2_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to phi2Z\to \tau\tau\ell\ell@f$.
 */
class Hobs_pp_phi3_phi2Z_tautaull_2_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi3_phi2Z_tautaull_2_CMS8 constructor.
     */
    Hobs_pp_phi3_phi2Z_tautaull_2_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi2Z\to \tau\tau\ell \ell)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to AZ\to \tau\tau\ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_phi2_phi3Z_tautaull_2_CMS8
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to \tau\tau\ell\ell@f$.
 */
class Hobs_pp_phi2_phi3Z_tautaull_2_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_phi2_phi3Z_tautaull_2_CMS8 constructor.
     */
    Hobs_pp_phi2_phi3Z_tautaull_2_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi3Z\to \tau\tau\ell \ell)]_{\text{theo}} / [\sigma_{pp\to phi3}\cdot BR(phi3\to AZ\to \tau\tau\ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi3_phi2Z_bbZ_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio 
 */
class Hobs_gg_phi3_phi2Z_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi3_phi2Z_bbZ_ATLAS13 constructor.
     */
    Hobs_gg_phi3_phi2Z_bbZ_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_phi2_phi3Z_bbZ_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio 
 */
class Hobs_gg_phi2_phi3Z_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_phi2_phi3Z_bbZ_ATLAS13 constructor.
     */
    Hobs_gg_phi2_phi3Z_bbZ_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi3_phi2Z_bbZ_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio 
 */
class Hobs_bb_phi3_phi2Z_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi3_phi2Z_bbZ_ATLAS13 constructor.
     */
    Hobs_bb_phi3_phi2Z_bbZ_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_bb_phi2_phi3Z_bbZ_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio 
 */
class Hobs_bb_phi2_phi3Z_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_phi2_phi3Z_bbZ_ATLAS13 constructor.
     */
    Hobs_bb_phi2_phi3Z_bbZ_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_Hpm_taunu_ATLAS8_GTHDM
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio 
 */
class Hobs_pp_Hpm_taunu_ATLAS8_GTHDM: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hpm_taunu_ATLAS8_GTHDM constructor.
     */
    Hobs_pp_Hpm_taunu_ATLAS8_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_Hp_taunu_CMS8_GTHDM
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio 
 */
class Hobs_pp_Hp_taunu_CMS8_GTHDM: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hp_taunu_CMS8_GTHDM constructor.
     */
    Hobs_pp_Hp_taunu_CMS8_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_Hpm_taunu_ATLAS13_GTHDM
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to \tau\tau\ell\ell@f$.
 */
class Hobs_pp_Hpm_taunu_ATLAS13_GTHDM: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hpm_taunu_ATLAS13_GTHDM constructor.
     */
    Hobs_pp_Hpm_taunu_ATLAS13_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_Hpm_taunu_CMS13_GTHDM
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to \tau\tau\ell\ell@f$.
 */
class Hobs_pp_Hpm_taunu_CMS13_GTHDM: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hpm_taunu_CMS13_GTHDM constructor.
     */
    Hobs_pp_Hpm_taunu_CMS13_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_Hpm_tb_ATLAS8_GTHDM
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to \tau\tau\ell\ell@f$.
 */
class Hobs_pp_Hpm_tb_ATLAS8_GTHDM: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hpm_tb_ATLAS8_GTHDM constructor.
     */
    Hobs_pp_Hpm_tb_ATLAS8_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_Hp_tb_CMS8_GTHDM
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to \tau\tau\ell\ell@f$.
 */
class Hobs_pp_Hp_tb_CMS8_GTHDM: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hp_tb_CMS8_GTHDM constructor.
     */
    Hobs_pp_Hp_tb_CMS8_GTHDM(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_Hpm_tb_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to phi3\to AZ\to \tau\tau\ell\ell@f$.
 */
class Hobs_pp_Hpm_tb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_Hpm_tb_ATLAS13 constructor.
     */
    Hobs_pp_Hpm_tb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_tt_phi2_tt_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$t\bar t\to phi2\to t\bar t@f$ at 13 TeV.
 */
class log10_tt_phi2_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_tt_phi2_tt_TH13 constructor.
     */
    log10_tt_phi2_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{t\bar t\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_tt_phi3_tt_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$t\bar t\to phi3\to t\bar t@f$ at 13 TeV.
 */
class log10_tt_phi3_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_tt_phi3_tt_TH13 constructor.
     */
    log10_tt_phi3_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{t\bar t\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_bb_phi2_tt_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi2\to t\bar t@f$ at 13 TeV.
 */
class log10_bb_phi2_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_phi2_tt_TH13 constructor.
     */
    log10_bb_phi2_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_bb_phi3_tt_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi3\to t\bar t@f$ at 13 TeV.
 */
class log10_bb_phi3_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_phi3_tt_TH13 constructor.
     */
    log10_bb_phi3_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to t\bar t)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi2_bb_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi2\to b\bar b@f$ at 8 TeV.
 */
class log10_bb_phi2_bb_TH8: public ThObservable {
public:

    /**
     * @brief log10_bb_phi2_bb_TH8 constructor.
     */
    log10_bb_phi2_bb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi3_bb_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi3\to b\bar b@f$ at 8 TeV.
 */
class log10_bb_phi3_bb_TH8: public ThObservable {
public:

    /**
     * @brief log10_bb_phi3_bb_TH8 constructor.
     */
    log10_bb_phi3_bb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_bb_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi2\to b\bar b@f$ at 8 TeV.
 */
class log10_gg_phi2_bb_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_bb_TH8 constructor.
     */
    log10_gg_phi2_bb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{g\bar g\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi3_bb_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi3\to b\bar b@f$ at 8 TeV.
 */
class log10_gg_phi3_bb_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_bb_TH8 constructor.
     */
    log10_gg_phi3_bb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{g\bar g\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi2_bb_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to b\bar b@f$ at 13 TeV.
 */
class log10_pp_phi2_bb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_bb_TH13 constructor.
     */
    log10_pp_phi2_bb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi3_bb_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to b\bar b@f$ at 13 TeV.
 */
class log10_pp_phi3_bb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_bb_TH13 constructor.
     */
    log10_pp_phi3_bb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi2_bb_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to b\bar b@f$ at 13 TeV.
 */
class log10_bb_phi2_bb_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_phi2_bb_TH13 constructor.
     */
    log10_bb_phi2_bb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar{b} \to phi2}\cdot BR^{\text{GTHDM}}(phi2\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi3_bb_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to b\bar b@f$ at 13 TeV.
 */
class log10_bb_phi3_bb_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_phi3_bb_TH13 constructor.
     */
    log10_bb_phi3_bb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar{b} \to phi3}\cdot BR^{\text{GTHDM}}(phi3\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_tautau_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to \tau\tau@f$ at 8 TeV.
 */
class log10_gg_phi2_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_tautau_TH8 constructor.
     */
    log10_gg_phi2_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
    
};


/**
 * @class log10_gg_phi3_tautau_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to \tau\tau@f$ at 8 TeV.
 */
class log10_gg_phi3_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_tautau_TH8 constructor.
     */
    log10_gg_phi3_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
    
};

/**
 * @class log10_bb_phi2_tautau_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi2\to \tau\tau@f$ at 8 TeV.
 */
class log10_bb_phi2_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_bb_phi2_tautau_TH8 constructor.
     */
    log10_bb_phi2_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi3_tautau_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi3\to \tau\tau@f$ at 8 TeV.
 */
class log10_bb_phi3_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_bb_phi3_tautau_TH8 constructor.
     */
    log10_bb_phi3_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_tautau_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to \tau \tau@f$ at 13 TeV.
 */
class log10_gg_phi2_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_tautau_TH13 constructor.
     */
    log10_gg_phi2_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau \tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class log10_gg_phi3_tautau_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to \tau \tau@f$ at 13 TeV.
 */
class log10_gg_phi3_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_tautau_TH13 constructor.
     */
    log10_gg_phi3_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau \tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi2_tautau_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi2\to \tau \tau@f$ at 13 TeV.
 */
class log10_bb_phi2_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_phi2_tautau_TH13 constructor.
     */
    log10_bb_phi2_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar b\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \tau \tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi3_tautau_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$b\bar b\to phi3\to \tau \tau@f$ at 13 TeV.
 */
class log10_bb_phi3_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_phi3_tautau_TH13 constructor.
     */
    log10_bb_phi3_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{b\bar b\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \tau \tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi2_gaga_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to \gamma\gamma@f$ at 8 TeV.
 */
class log10_gg_phi2_gaga_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_gaga_TH8 constructor.
     */
    log10_gg_phi2_gaga_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \gamma\gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi3_gaga_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to \gamma\gamma@f$ at 8 TeV.
 */
class log10_gg_phi3_gaga_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_gaga_TH8 constructor.
     */
    log10_gg_phi3_gaga_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \gamma\gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi2_gaga_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to \gamma \gamma@f$ at 13 TeV.
 */
class log10_pp_phi2_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_gaga_TH13 constructor.
     */
    log10_pp_phi2_gaga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi3_gaga_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to \gamma \gamma@f$ at 13 TeV.
 */
class log10_pp_phi3_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_gaga_TH13 constructor.
     */
    log10_pp_phi3_gaga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_gaga_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to \gamma \gamma@f$ at 13 TeV.
 */
class log10_gg_phi2_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_gaga_TH13 constructor.
     */
    log10_gg_phi2_gaga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi3_gaga_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to \gamma \gamma@f$ at 13 TeV.
 */
class log10_gg_phi3_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_gaga_TH13 constructor.
     */
    log10_gg_phi3_gaga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi2_Zga_llga_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to Z\gamma \to \ell \ell \gamma@f$ at 8 TeV.
 */
class log10_pp_phi2_Zga_llga_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_Zga_llga_TH8 constructor.
     */
    log10_pp_phi2_Zga_llga_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to Z\gamma \to \ell \ell \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi3_Zga_llga_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to Z\gamma \to \ell \ell \gamma@f$ at 8 TeV.
 */
class log10_pp_phi3_Zga_llga_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_Zga_llga_TH8 constructor.
     */
    log10_pp_phi3_Zga_llga_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to Z\gamma \to \ell \ell \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_Zga_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to Z\gamma@f$ at 13 TeV.
 */
class log10_gg_phi2_Zga_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_Zga_TH13 constructor.
     */
    log10_gg_phi2_Zga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to Z\gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi3_Zga_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to Z\gamma@f$ at 13 TeV.
 */
class log10_gg_phi3_Zga_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_Zga_TH13 constructor.
     */
    log10_gg_phi3_Zga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to Z\gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_ZZ_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to ZZ@f$ at 8 TeV.
 */
class log10_gg_phi2_ZZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_ZZ_TH8 constructor.
     */
    log10_gg_phi2_ZZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi3_ZZ_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to ZZ@f$ at 8 TeV.
 */
class log10_gg_phi3_ZZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_ZZ_TH8 constructor.
     */
    log10_gg_phi3_ZZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_VV_phi2_ZZ_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to phi2\to ZZ@f$ at 8 TeV.
 */
class log10_VV_phi2_ZZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_VV_phi2_ZZ_TH8 constructor.
     */
    log10_VV_phi2_ZZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{VV\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_VV_phi3_ZZ_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to phi3\to ZZ@f$ at 8 TeV.
 */
class log10_VV_phi3_ZZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_VV_phi3_ZZ_TH8 constructor.
     */
    log10_VV_phi3_ZZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{VV\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_ZZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to ZZ@f$ at 13 TeV.
 */
class log10_gg_phi2_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_ZZ_TH13 constructor.
     */
    log10_gg_phi2_ZZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi3_ZZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to ZZ@f$ at 13 TeV.
 */
class log10_gg_phi3_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_ZZ_TH13 constructor.
     */
    log10_gg_phi3_ZZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class log10_VV_phi2_ZZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to phi2\to ZZ@f$ at 13 TeV.
 */
class log10_VV_phi2_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_VV_phi2_ZZ_TH13 constructor.
     */
    log10_VV_phi2_ZZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{VV\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_VV_phi3_ZZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to phi3\to ZZ@f$ at 13 TeV.
 */
class log10_VV_phi3_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_VV_phi3_ZZ_TH13 constructor.
     */
    log10_VV_phi3_ZZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{VV\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi2_ZZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to ZZ@f$ at 13 TeV.
 */
class log10_pp_phi2_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_ZZ_TH13 constructor.
     */
    log10_pp_phi2_ZZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_ZZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to ZZ@f$ at 13 TeV.
 */
class log10_pp_phi3_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_ZZ_TH13 constructor.
     */
    log10_pp_phi3_ZZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi2_WW_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to WW@f$ at 8 TeV.
 */
class log10_gg_phi2_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_WW_TH8 constructor.
     */
    log10_gg_phi2_WW_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi3_WW_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to WW@f$ at 8 TeV.
 */
class log10_gg_phi3_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_WW_TH8 constructor.
     */
    log10_gg_phi3_WW_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_VV_phi2_WW_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to phi2\to WW@f$ at 8 TeV.
 */
class log10_VV_phi2_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_VV_phi2_WW_TH8 constructor.
     */
    log10_VV_phi2_WW_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{VV\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_VV_phi3_WW_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to phi3\to WW@f$ at 8 TeV.
 */
class log10_VV_phi3_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_VV_phi3_WW_TH8 constructor.
     */
    log10_VV_phi3_WW_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{VV\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_WW_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to WW@f$ at 13 TeV.
 */
class log10_gg_phi2_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_WW_TH13 constructor.
     */
    log10_gg_phi2_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi3_WW_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to WW@f$ at 13 TeV.
 */
class log10_gg_phi3_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_WW_TH13 constructor.
     */
    log10_gg_phi3_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_VV_phi2_WW_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to phi2\to WW@f$ at 13 TeV.
 */
class log10_VV_phi2_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_VV_phi2_WW_TH13 constructor.
     */
    log10_VV_phi2_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{VV\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_VV_phi3_WW_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$VV\to phi3\to WW@f$ at 13 TeV.
 */
class log10_VV_phi3_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_VV_phi3_WW_TH13 constructor.
     */
    log10_VV_phi3_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{VV\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_ggVV_phi2_WW_lnulnu_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$(gg+VV)\to phi2\to WW\to \ell \nu \ell \nu@f$ at 13 TeV.
 */
class log10_ggVV_phi2_WW_lnulnu_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggVV_phi2_WW_lnulnu_TH13 constructor.
     */
    log10_ggVV_phi2_WW_lnulnu_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{(gg+VV)\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to WW\to \ell \nu \ell \nu)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_ggVV_phi3_WW_lnulnu_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$(gg+VV)\to phi3\to WW\to \ell \nu \ell \nu@f$ at 13 TeV.
 */
class log10_ggVV_phi3_WW_lnulnu_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggVV_phi3_WW_lnulnu_TH13 constructor.
     */
    log10_ggVV_phi3_WW_lnulnu_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{(gg+VV)\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to WW\to \ell \nu \ell \nu)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi2_WW_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$(pp)\to phi2\to WWf$ at 13 TeV.
 */
class log10_pp_phi2_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_WW_TH13 constructor.
     */
    log10_pp_phi2_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_WW_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$(pp)\to phi3\to WWf$ at 13 TeV.
 */
class log10_pp_phi3_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_WW_TH13 constructor.
     */
    log10_pp_phi3_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to WW)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_mu_pp_phi2_VV_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the signal strength of the process @f$pp\to phi2\to VV@f$ at 8 TeV.
 */
class log10_mu_pp_phi2_VV_TH8: public ThObservable {
public:

    /**
     * @brief log10_mu_pp_phi2_VV_TH8 constructor.
     */
    log10_mu_pp_phi2_VV_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\mu_H^{\text{GTHDM}}(phi2\to VV)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_mu_pp_phi3_VV_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the signal strength of the process @f$pp\to phi3\to VV@f$ at 8 TeV.
 */
class log10_mu_pp_phi3_VV_TH8: public ThObservable {
public:

    /**
     * @brief log10_mu_pp_phi3_VV_TH8 constructor.
     */
    log10_mu_pp_phi3_VV_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\mu_H^{\text{GTHDM}}(phi3\to VV)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi2_VV_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to WW,ZZ@f$ at 13 TeV.
 */
class log10_pp_phi2_VV_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_VV_TH13 constructor.
     */
    log10_pp_phi2_VV_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot [BR^{\text{GTHDM}}(phi2\to WW)+BR^{\text{GTHDM}}(phi2\to ZZ)]]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_VV_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to WW,ZZ@f$ at 13 TeV.
 */
class log10_pp_phi3_VV_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_VV_TH13 constructor.
     */
    log10_pp_phi3_VV_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot [BR^{\text{GTHDM}}(phi3\to WW)+BR^{\text{GTHDM}}(phi3\to ZZ)]]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_phi1phi1_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to phi1phi1@f$ at 8 TeV.
 */
class log10_gg_phi2_phi1phi1_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_phi1phi1_TH8 constructor.
     */
    log10_gg_phi2_phi1phi1_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi3_phi1phi1_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to phi1phi1@f$ at 8 TeV.
 */
class log10_gg_phi3_phi1phi1_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_phi1phi1_TH8 constructor.
     */
    log10_gg_phi3_phi1phi1_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi2_phi1phi1_bbbb_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b b\bar b@f$ at 8 TeV.
 */
class log10_pp_phi2_phi1phi1_bbbb_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_phi1phi1_bbbb_TH8 constructor.
     */
    log10_pp_phi2_phi1phi1_bbbb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi3_phi1phi1_bbbb_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b b\bar b@f$ at 8 TeV.
 */
class log10_pp_phi3_phi1phi1_bbbb_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_phi1phi1_bbbb_TH8 constructor.
     */
    log10_pp_phi3_phi1phi1_bbbb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class log10_pp_phi2_phi1phi1_bbgaga_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \gamma \gamma@f$ at 8 TeV.
 */
class log10_pp_phi2_phi1phi1_bbgaga_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_phi1phi1_bbgaga_TH8 constructor.
     */
    log10_pp_phi2_phi1phi1_bbgaga_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi3_phi1phi1_bbgaga_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \gamma \gamma@f$ at 8 TeV.
 */
class log10_pp_phi3_phi1phi1_bbgaga_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_phi1phi1_bbgaga_TH8 constructor.
     */
    log10_pp_phi3_phi1phi1_bbgaga_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi2_phi1phi1_bbtautau_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi2\to phi1phi1\to b\bar b \tau\tau@f$ at 8 TeV.
 */
class log10_gg_phi2_phi1phi1_bbtautau_TH8: public ThObservable {
 public:

  /**                                                                                                                                         
   * @brief log10_gg_phi2_phi1phi1_bbtautau_TH8 constructor.                                                                                                                      
   */
  log10_gg_phi2_phi1phi1_bbtautau_TH8(const StandardModel& SM_i);

  /**                                                                                                                                         
   * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi2phi2\to b\bar b \tau\tau)]@f$
   */
  double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi3_phi1phi1_bbtautau_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to phi3\to phi1phi1\to b\bar b \tau\tau@f$ at 8 TeV.
 */
class log10_gg_phi3_phi1phi1_bbtautau_TH8: public ThObservable {
 public:

  /**                                                                                                                                         
   * @brief log10_gg_phi3_phi1phi1_bbtautau_TH8 constructor.                                                                                                                      
   */
  log10_gg_phi3_phi1phi1_bbtautau_TH8(const StandardModel& SM_i);

  /**                                                                                                                                         
   * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{gg\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi2phi2\to b\bar b \tau\tau)]@f$
   */
  double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi2_phi1phi1_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1@f$ at 8 TeV.
 */
class log10_pp_phi2_phi1phi1_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_phi1phi1_TH8 constructor.
     */
    log10_pp_phi2_phi1phi1_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_phi1phi1_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1@f$ at 8 TeV.
 */
class log10_pp_phi3_phi1phi1_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_phi1phi1_TH8 constructor.
     */
    log10_pp_phi3_phi1phi1_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi2_phi1phi1_bbbb_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b b\bar b@f$ at 13 TeV.
 */
class log10_pp_phi2_phi1phi1_bbbb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_phi1phi1_bbbb_TH13 constructor.
     */
    log10_pp_phi2_phi1phi1_bbbb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_phi1phi1_bbbb_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b b\bar b@f$ at 13 TeV.
 */
class log10_pp_phi3_phi1phi1_bbbb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_phi1phi1_bbbb_TH13 constructor.
     */
    log10_pp_phi3_phi1phi1_bbbb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi2_phi1phi1_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1@f$ at 13 TeV.
 */
class log10_pp_phi2_phi1phi1_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_phi1phi1_TH13 constructor.
     */
    log10_pp_phi2_phi1phi1_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_phi1phi1_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1@f$ at 13 TeV.
 */
class log10_pp_phi3_phi1phi1_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_phi1phi1_TH13 constructor.
     */
    log10_pp_phi3_phi1phi1_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class log10_pp_phi2_phi1phi1_bbgaga_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1 \to b \bar b \gamma \gamma @f$ at 13 TeV.
 */
class log10_pp_phi2_phi1phi1_bbgaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_phi1phi1_bbgaga_TH13 constructor.
     */
    log10_pp_phi2_phi1phi1_bbgaga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1 \to b \bar b \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_phi1phi1_bbgaga_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1 \to b \bar b \gamma \gamma @f$ at 13 TeV.
 */
class log10_pp_phi3_phi1phi1_bbgaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_phi1phi1_bbgaga_TH13 constructor.
     */
    log10_pp_phi3_phi1phi1_bbgaga_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1 \to b \bar b \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi2_phi1phi1_bbtautau_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b \tau \tau@f$ at 13 TeV.
 */
class log10_pp_phi2_phi1phi1_bbtautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_phi1phi1_bbtautau_TH13 constructor.
     */
    log10_pp_phi2_phi1phi1_bbtautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b \tau \tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_phi1phi1_bbtautau_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b \tau \tau@f$ at 13 TeV.
 */
class log10_pp_phi3_phi1phi1_bbtautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_phi1phi1_bbtautau_TH13 constructor.
     */
    log10_pp_phi3_phi1phi1_bbtautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b \tau \tau)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi2_phi1phi1_bbVV_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to b\bar b VV\to b\bar b \ell \ell \nu \nu@f$ at 13 TeV.
 */
class log10_pp_phi2_phi1phi1_bbVV_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi2_phi1phi1_bbVV_TH13 constructor.
     */
    log10_pp_phi2_phi1phi1_bbVV_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1\to b\bar b WW\to b\bar b \ell \ell \nu \nu)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_phi1phi1_bbVV_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to b\bar b VV\to b\bar b \ell \ell \nu \nu@f$ at 13 TeV.
 */
class log10_pp_phi3_phi1phi1_bbVV_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_phi1phi1_bbVV_TH13 constructor.
     */
    log10_pp_phi3_phi1phi1_bbVV_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1\to b\bar b WW\to b\bar b \ell \ell \nu \nu)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi2_phi1phi1_gagaWW_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1phi1\to \gamma \gamma WW @f$ at 13 TeV.
 */
class log10_gg_phi2_phi1phi1_gagaWW_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_phi1phi1_gagaWW_TH13 constructor.
     */
    log10_gg_phi2_phi1phi1_gagaWW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1phi1 \to \gamma \gamma WW )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi3_phi1phi1_gagaWW_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1phi1\to \gamma \gamma WW @f$ at 13 TeV.
 */
class log10_gg_phi3_phi1phi1_gagaWW_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_phi1phi1_gagaWW_TH13 constructor.
     */
    log10_gg_phi3_phi1phi1_gagaWW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1phi1 \to \gamma \gamma WW )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_phi1Z_bbZ_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1Z\to bbZ @f$ at 8 TeV.
 */
class log10_gg_phi2_phi1Z_bbZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_phi1Z_bbZ_TH8 constructor.
     */
    log10_gg_phi2_phi1Z_bbZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1Z\to bbZ )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi3_phi1Z_bbZ_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1Z\to bbZ @f$ at 8 TeV.
 */
class log10_gg_phi3_phi1Z_bbZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_phi1Z_bbZ_TH8 constructor.
     */
    log10_gg_phi3_phi1Z_bbZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi3\to phi1Z \to bbZ )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_phi1Z_bbll_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1Z \to bbll @f$ at 8 TeV.
 */
class log10_gg_phi2_phi1Z_bbll_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_phi1Z_bbll_TH8 constructor.
     */
    log10_gg_phi2_phi1Z_bbll_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1Z \to bbll )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi3_phi1Z_bbll_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1Z \to bbll @f$ at 8 TeV.
 */
class log10_gg_phi3_phi1Z_bbll_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_phi1Z_bbll_TH8 constructor.
     */
    log10_gg_phi3_phi1Z_bbll_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1Z \to bbll )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_phi1Z_tautauZ_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1Z \to tautauZ @f$ at 8 TeV.
 */
class log10_gg_phi2_phi1Z_tautauZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_phi1Z_tautauZ_TH8 constructor.
     */
    log10_gg_phi2_phi1Z_tautauZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi2\to phi1Z \to tautauZ )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi3_phi1Z_tautauZ_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1Z \to tautauZ @f$ at 8 TeV.
 */
class log10_gg_phi3_phi1Z_tautauZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_phi1Z_tautauZ_TH8 constructor.
     */
    log10_gg_phi3_phi1Z_tautauZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1Z \to tautauZ )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_gg_phi2_phi1Z_tautaull_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1Z \to tautaull @f$ at 8 TeV.
 */
class log10_gg_phi2_phi1Z_tautaull_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_phi1Z_tautaull_TH8 constructor.
     */
    log10_gg_phi2_phi1Z_tautaull_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1Z \to tautaull )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi3_phi1Z_tautaull_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1Z \to tautaull @f$ at 13 TeV.
 */
class log10_gg_phi3_phi1Z_tautaull_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_phi1Z_tautaull_TH8 constructor.
     */
    log10_gg_phi3_phi1Z_tautaull_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1Z \to tautaull )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi2_phi1Z_bbZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi2\to phi1Z \to bbZ @f$ at 13 TeV.
 */
class log10_gg_phi2_phi1Z_bbZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi2_phi1Z_bbZ_TH13 constructor.
     */
    log10_gg_phi2_phi1Z_bbZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1Z \to bbZ )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class log10_gg_phi3_phi1Z_bbZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi1Z \to bbZ @f$ at 13 TeV.
 */
class log10_gg_phi3_phi1Z_bbZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_phi1Z_bbZ_TH13 constructor.
     */
    log10_gg_phi3_phi1Z_bbZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1Z \to bbZ )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi2_phi1Z_bbZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$bb\to phi2\to phi1Z \to bbZ @f$ at 13 TeV.
 */
class log10_bb_phi2_phi1Z_bbZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_phi2_phi1Z_bbZ_TH13 constructor.
     */
    log10_bb_phi2_phi1Z_bbZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{bb\to phi2}\cdot BR^{\text{GTHDM}}(phi2\to phi1Z \to bbZ )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi3_phi1Z_bbZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$bb\to phi3\to phi1Z \to bbZ @f$ at 13 TeV.
 */
class log10_bb_phi3_phi1Z_bbZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_phi2_phi1Z_bbZ_TH13 constructor.
     */
    log10_bb_phi3_phi1Z_bbZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{bb\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi1Z \to bbZ )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_phi3_AZ_bbll_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi2Z\to b\bar b \ell \ell@f$ at 8 TeV.
 */
class log10_pp_phi3_phi2Z_bbll_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_AZ_bbll_TH8 constructor.
     */
    log10_pp_phi3_phi2Z_bbll_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to AZ\to b\bar b \ell \ell)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_phi3_phi2Z_tautaull_TH8
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to phi2Z \to \tau \tau \ell \ell@f$ at 8 TeV.
 */
class log10_pp_phi3_phi2Z_tautaull_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_phi3_phi2Z_tautaull_TH8 constructor.
     */
    log10_pp_phi3_phi2Z_tautaull_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to phi2 Z\to \tau tau \ell \ell)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_gg_phi3_phi2Z_bbZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to phi3\to \phi2Z \to bbZ @f$ at 13 TeV.
 */
class log10_gg_phi3_phi2Z_bbZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_phi3_phi2Z_bbZ_TH13 constructor.
     */
    log10_gg_phi3_phi2Z_bbZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \phi2 Z\to b\bar bZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_bb_phi3_phi2Z_bbZ_TH13
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$bb\to phi3\to \phi2Z \to bbZ @f$ at 13 TeV.
 */
class log10_bb_phi3_phi2Z_bbZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_phi3_phi2Z_bbZ_TH13 constructor.
     */
    log10_bb_phi3_phi2Z_bbZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{bb\to phi3}\cdot BR^{\text{GTHDM}}(phi3\to \phi2 Z\to b\bar bZ)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_Hpm_taunu_TH8_GTHDM
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^{\pm} \to \tau \nu @f$ at 8 TeV.
 */
class log10_pp_Hpm_taunu_TH8_GTHDM: public ThObservable {
public:

    /**
     * @brief log10_pp_Hpm_taunu_TH8_GTHDM constructor.
     */
    log10_pp_Hpm_taunu_TH8_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to H^{\pm}}\cdot BR^{\text{GTHDM}}(H^{\pm}\to \phi2 \tau \nu)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_Hp_taunu_TH8_GTHDM
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^{\+} \to \tau \nu @f$ at 8 TeV.
 */
class log10_pp_Hp_taunu_TH8_GTHDM: public ThObservable {
public:

    /**
     * @brief log10_pp_Hp_taunu_TH8_GTHDM constructor.
     */
    log10_pp_Hp_taunu_TH8_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to H^{+}}\cdot BR^{\text{GTHDM}}(H^+\to \phi2 \tau \nu)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_Hpm_taunu_TH13_GTHDM
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^{\pm} \to \tau \nu @f$ at 13 TeV.
 */
class log10_pp_Hpm_taunu_TH13_GTHDM: public ThObservable {
public:

    /**
     * @brief log10_pp_Hpm_taunu_TH13 constructor.
     */
    log10_pp_Hpm_taunu_TH13_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to H^{\pm}}\cdot BR^{\text{GTHDM}}(H\to \phi2 \tau \nu)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_Hpm_tb_TH8_GTHDM
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^{\pm} \to t b  @f$ at 8 TeV.
 */
class log10_pp_Hpm_tb_TH8_GTHDM: public ThObservable {
public:

    /**
     * @brief log10_pp_Hpm_tb_TH8_GTHDM constructor.
     */
    log10_pp_Hpm_tb_TH8_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to H^{\pm}}\cdot BR^{\text{GTHDM}}(H^{\pm}\to \phi2 t )]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class log10_pp_Hp_tb_TH8_GTHDM
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^{\+} \to t b @f$ at 8 TeV.
 */
class log10_pp_Hp_tb_TH8_GTHDM: public ThObservable {
public:

    /**
     * @brief log10_pp_Hp_tb_TH8_GTHDM constructor.
     */
    log10_pp_Hp_tb_TH8_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to H^{+}}\cdot BR^{\text{GTHDM}}(H^+\to \phi2 t b)]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class log10_pp_Hpm_tb_TH13_GTHDM
 * @ingroup GeneralTHDM
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$pp\to H^{\pm} \to t b @f$ at 13 TeV.
 */
class log10_pp_Hpm_tb_TH13_GTHDM: public ThObservable {
public:

    /**
     * @brief log10_pp_Hpm_taunu_TH13_GTHDM constructor.
     */
    log10_pp_Hpm_tb_TH13_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GTHDM}}_{pp\to H^{\pm}}\cdot BR^{\text{GTHDM}}(H\to \phi2 t b]@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

//NOT CLEANED YET


/**
 * @class Gamma_phi3_GTHDM
 * @ingroup GeneralTHDM
 * @brief Total phi3 decay rate in the %GTHDM.
 */
class Gamma_phi3_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Gamma_phi3_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$\Gamma_phi3@f$ in units of GeV
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

///**
// * @class rHH_gaga_GTHDM
// * @ingroup GeneralTHDM
// * @brief Squared relative coupling of @f$H@f$ to two photons.
// */
//class rHH_gaga_GTHDM : public ThObservable {
//public:
//    
//    /**
//     * @brief Constructor.
//     */
//    rHH_gaga_GTHDM(const StandardModel& SM_i);
//    
//    /**
//     * @return @f$r^{(H)}_{\gamma \gamma}@f$
//     */
//    double computeThValue ();
//};

/**
 * @class rHH_gg_GTHDM
 * @ingroup GeneralTHDM
 * @brief Squared relative coupling of @f$H@f$ to two gluons.
 */
class rphi3_ggE_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor for the squared relative coupling of @f$H@f$ to two gluons.
     */
    rphi3_ggE_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(H)}_{gg}@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

class rphi3_ggO_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor for the squared relative coupling of @f$H@f$ to two gluons.
     */
    rphi3_ggO_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(H)}_{gg}@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_HH_phi1phi1_GTHDM
 * @ingroup GeneralTHDM
 * @brief %GTHDM branching ratio of @f$H@f$ to two @f$h@f$.
 */
class BR_phi3_phi1phi1_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi3_phi1phi1_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi3\to phi1phi1)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_HH_phi2phi2_GTHDM
 * @ingroup GeneralTHDM
 * @brief %GTHDM branching ratio of @f$H@f$ to two @f$h@f$.
 */
class BR_phi3_phi2phi2_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi3_phi2phi2_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi3\to phi2phi2)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_HH_HpHm_GTHDM
 * @ingroup GeneralTHDM
 * @brief %GTHDM branching ratio of @f$H@f$ to charged Higgs bosons.
 */
class BR_phi3_HpHm_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi3_HpHm_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi3\to phi3^\pm H^\mp)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_HH_AZ_GTHDM
 * @ingroup GeneralTHDM
 * @brief %GTHDM branching ratio of @f$H@f$ to an @f$A@f$ and a @f$Z@f$ boson.
 */
class BR_phi3_phi1Z_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi3_phi1Z_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi3\to AZ)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_HH_AZ_GTHDM
 * @ingroup GeneralTHDM
 * @brief %THDM branching ratio of @f$H@f$ to an @f$A@f$ and a @f$Z@f$ boson.
 */
class BR_phi3_phi2Z_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi3_phi2Z_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi3\to AZ)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_HH_HpW_GTHDM
 * @ingroup GeneralTHDM
 * @brief %THDM branching ratio of @f$H@f$ to a charged Higgs boson and a @f$W@f$ boson.
 */
class BR_phi3_HpW_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi3_HpW_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi3\to phi3^\pm W^\mp)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Gamma_phi2_GTHDM
 * @ingroup GeneralTHDM
 * @brief Total phi2 decay rate in the %GTHDM.
 */
class Gamma_phi2_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Gamma_phi2_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$\Gamma_phi2@f$ in units of GeV
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

///**
// * @class rHH_gaga_GTHDM
// * @ingroup GeneralTHDM
// * @brief Squared relative coupling of @f$H@f$ to two photons.
// */
//class rHH_gaga_GTHDM : public ThObservable {
//public:
//    
//    /**
//     * @brief Constructor.
//     */
//    rHH_gaga_GTHDM(const StandardModel& SM_i);
//    
//    /**
//     * @return @f$r^{(H)}_{\gamma \gamma}@f$
//     */
//    double computeThValue ();
//};

/**
 * @class rHH_gg_GTHDM
 * @ingroup GeneralTHDM
 * @brief Squared relative coupling of @f$H@f$ to two gluons.
 */
class rphi2_ggE_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor for the squared relative coupling of @f$H@f$ to two gluons.
     */
    rphi2_ggE_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(H)}_{gg}@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

class rphi2_ggO_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor for the squared relative coupling of @f$H@f$ to two gluons.
     */
    rphi2_ggO_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(H)}_{gg}@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_HH_phi1phi1_GTHDM
 * @ingroup GeneralTHDM
 * @brief %GTHDM branching ratio of @f$H@f$ to two @f$h@f$.
 */
class BR_phi2_phi1phi1_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi2_phi1phi1_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi2\to phi1phi1)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class BR_HH_HpHm_GTHDM
 * @ingroup GeneralTHDM
 * @brief %GTHDM branching ratio of @f$H@f$ to charged Higgs bosons.
 */
class BR_phi2_HpHm_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi2_HpHm_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi2\to phi2^\pm H^\mp)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_HH_AZ_GTHDM
 * @ingroup GeneralTHDM
 * @brief %GTHDM branching ratio of @f$H@f$ to an @f$A@f$ and a @f$Z@f$ boson.
 */
class BR_phi2_phi1Z_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi2_phi1Z_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi2\to AZ)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_HH_HpW_GTHDM
 * @ingroup GeneralTHDM
 * @brief %THDM branching ratio of @f$H@f$ to a charged Higgs boson and a @f$W@f$ boson.
 */
class BR_phi2_HpW_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_phi2_HpW_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(phi2\to phi2^\pm W^\mp)@f$
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};


#endif	/*GENERALTHDMHEAVYHIGGS_H*/