/*
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMLOWMASS_H
#define GENERALTHDMLOWMASS_H

#include "ThObservable.h"

class GeneralTHDM;
class GeneralTHDMcache;

/**
 * @class GeneralTHDMLowMass
 * @ingroup GeneralTHDM 
 * @brief Base class for low mass Higgs direct search observables.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theory over experiment ratios of
 * several phenomena that include an extra scalar lighter than the SM Higgs
 */

/*********************************/
/* Observables with phi_3, i.e A */
/*********************************/

/**
 * @class Hobs_pp_h_phi3phi3_mumutautau_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to h_{125}\to AA\to \mu\mu\tau\tau@f$.
 */
class Hobs_pp_h_phi3phi3_mumutautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_mumutautau_CMS13 constructor.
     */
    Hobs_pp_h_phi3phi3_mumutautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM/SM}}_{pp\to h_{125}}\cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to \mu\mu) \cdot BR^{\text{GTHDM}}(A\to \tau\tau)]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi3phi3_bbtautau_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to h_{125}\to AA\to bb\tau\tau@f$.
 */
class Hobs_pp_h_phi3phi3_bbtautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_bbtautau_CMS13 constructor.
     */
    Hobs_pp_h_phi3phi3_bbtautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM/SM}}_{pp\to h_{125}}\cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to bb) \cdot BR^{\text{GTHDM}}(A\to \tau\tau)]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi3phi3_bbmumu_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to h_{125}\to AA\to bb\tau\tau@f$.
 */
class Hobs_pp_h_phi3phi3_bbmumu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_bbmumu_CMS13 constructor.
     */
    Hobs_pp_h_phi3phi3_bbmumu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM/SM}}_{pp\to h_{125}}\cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to bb) \cdot BR^{\text{GTHDM}}(A\to \mu\mu)]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi3Z_mumull_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the branching ratio of the process @f$ h_{125}\to AZ\to \mu\mu ll@f$.
 */
class Hobs_pp_h_phi3Z_mumull_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3Z_mumull_CMS13 constructor.
     */
    Hobs_pp_h_phi3Z_mumull_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[BR^{\text{GTHDM}}(h_{125}\to AZ) \cdot BR^{\text{GTHDM}}(A\to \mu\mu)]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi3phi3_mumumumu_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the branching ratio of the process @f$h_{125}\to AA\to \mu\mu\mu\mu@f$.
 */
class Hobs_pp_h_phi3phi3_mumumumu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_mumumumu_CMS13 constructor.
     */
    Hobs_pp_h_phi3phi3_mumumumu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to \mu\mu) \cdot BR^{\text{GTHDM}}(A\to \mu\mu)]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi3phi3_gagagaga_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to h_{125}\to AA\to \gamma\gamma\gamma\gamma@f$.
 */
class Hobs_pp_h_phi3phi3_gagagaga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_gagagaga_CMS13 constructor.
     */
    Hobs_pp_h_phi3phi3_gagagaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to h_{125}}\cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to \gamma\gamma) \cdot BR^{\text{GTHDM}}(A\to \gamma\gamma)]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi3phi3_tautautautau_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross section times branching ratio of the process @f$pp\to h_{125}\to AA\to \tau\tau\tau\tau@f$.
 */
class Hobs_pp_h_phi3phi3_tautautautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_tautautautau_CMS13 constructor.
     */
    Hobs_pp_h_phi3phi3_tautautautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM/SM}}_{pp\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to \tau\tau) \cdot BR^{\text{GTHDM}}(A\to \tau\tau)]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi3phi3_bbmumu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the branching ratio of the process @f$h_{125}\to AA\to bb\mu\mu@f$.
 */
class Hobs_pp_h_phi3phi3_bbmumu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_bbmumu_ATLAS13 constructor.
     */
    Hobs_pp_h_phi3phi3_bbmumu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to bb) \cdot BR^{\text{GTHDM}}(A\to \mu\mu)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_h_phi3phi3_mumumumu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the branching ratio of the process @f$gg \to h_{125}\to AA\to \mu\mu\mu\mu@f$.
 */
class Hobs_gg_h_phi3phi3_mumumumu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_h_phi3phi3_mumumumu_ATLAS13 constructor.
     */
    Hobs_gg_h_phi3phi3_mumumumu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to \mu\mu) \cdot BR^{\text{GTHDM}}(A\to \mu\mu)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_h_phi3Z_mumull_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the branching ratio of the process @f$gg \to h_{125}\to AZ\to \mu\mull@f$.
 */
class Hobs_gg_h_phi3Z_mumull_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_h_phi3Z_mumull_ATLAS13 constructor.
     */
    Hobs_gg_h_phi3Z_mumull_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to AZ) \cdot BR^{\text{GTHDM}}(A\to \mu\mu) \cdot BR^{\text{GTHDM}}(Z\to ll)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_Vh_h_phi3phi3_bbbb_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the branching ratio of the process @f$Vh \to h_{125}\to AA\to bbbb@f$.
 */
class Hobs_Vh_h_phi3phi3_bbbb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_Vh_h_phi3phi3_bbbb_ATLAS13 constructor.
     */
    Hobs_Vh_h_phi3phi3_bbbb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{Vh\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to bb) \cdot BR^{\text{GTHDM}}(A\to bb)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_Zh_h_phi3phi3_bbbb_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the branching ratio of the process @f$Zh \to h_{125}\to AA\to bbbb@f$.
 */
class Hobs_Zh_h_phi3phi3_bbbb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_Zh_h_phi3phi3_bbbb_ATLAS13 constructor.
     */
    Hobs_Zh_h_phi3phi3_bbbb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{Zh\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to bb) \cdot BR^{\text{GTHDM}}(A\to bb)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi3phi3_bbmumu_ATLAS13_old
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to h_{125}\to AA\to bb\mu\mu@f$ at 13 TeV.
 */
class Hobs_pp_h_phi3phi3_bbmumu_ATLAS13_old: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_bbmumu_ATLAS13_old constructor.
     */
    Hobs_pp_h_phi3phi3_bbmumu_ATLAS13_old(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to bb) \cdot BR^{\text{GTHDM}}(A\to \mu\mu)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi3phi3_gagagg_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section  times branching ratio of the process @f$pp \to h_{125}\to AA\to \gamma\gamma gg@f$.
 */
class Hobs_pp_h_phi3phi3_gagagg_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi3phi3_gagagg_ATLAS13 constructor.
     */
    Hobs_pp_h_phi3phi3_gagagg_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to AA) \cdot BR^{\text{GTHDM}}(A\to \gamma\gamma) \cdot BR^{\text{GTHDM}}(A\to gg)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/*********************************/
/* Observables with phi_2, i.e H */
/*********************************/

/**
 * @class Hobs_pp_h_phi2Z_mumull_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the branching ratio of the process @f$ h_{125}\to HZ\to \mu\mu ll@f$.
 */
class Hobs_pp_h_phi2Z_mumull_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi2Z_mumull_CMS13 constructor.
     */
    Hobs_pp_h_phi2Z_mumull_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[BR^{\text{GTHDM}}(h_{125}\to HZ) \cdot BR^{\text{GTHDM}}(H\to \mu\mu) ]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Hobs_pp_h_phi2phi2_mumumumu_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the branching ratio of the process @f$ h_{125}\to HH\to \mu\mu\mu\mu@f$.
 */
class Hobs_pp_h_phi2phi2_mumumumu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi2phi2_mumumumu_CMS13 constructor.
     */
    Hobs_pp_h_phi2phi2_mumumumu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[BR^{\text{GTHDM}}(h_{125}\to HH) \cdot BR^{\text{GTHDM}}(H\to \mu\mu) \cdot BR^{\text{GTHDM}}(H\to \mu\mu)]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_h_phi2phi2_mumumumu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the branching ratio of the process @f$gg \to h_{125}\to HH\to \mu\mu\mu\mu@f$.
 */
class Hobs_gg_h_phi2phi2_mumumumu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_h_phi2phi2_mumumumu_ATLAS13 constructor.
     */
    Hobs_gg_h_phi2phi2_mumumumu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to HH) \cdot BR^{\text{GTHDM}}(H\to \mu\mu) \cdot BR^{\text{GTHDM}}(H\to \mu\mu)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_gg_h_phi3Z_mumull_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the branching ratio of the process @f$gg \to h_{125}\to AZ\to \mu\mull@f$.
 */
class Hobs_gg_h_phi2Z_mumull_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_h_phi2Z_mumull_ATLAS13 constructor.
     */
    Hobs_gg_h_phi2Z_mumull_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{gg\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to HZ) \cdot BR^{\text{GTHDM}}(H\to \mu\mu) \cdot BR^{\text{GTHDM}}(Z\to ll)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_Vh_h_phi2phi2_bbbb_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the branching ratio of the process @f$Vh \to h_{125}\to HH\to bbbb@f$.
 */
class Hobs_Vh_h_phi2phi2_bbbb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_Vh_h_phi2phi2_bbbb_ATLAS13 constructor.
     */
    Hobs_Vh_h_phi2phi2_bbbb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{Vh\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to HH) \cdot BR^{\text{GTHDM}}(H\to bb) \cdot BR^{\text{GTHDM}}(H\to bb)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_Zh_h_phi2phi2_bbbb_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the branching ratio of the process @f$Zh \to h_{125}\to HH\to bbbb@f$.
 */
class Hobs_Zh_h_phi2phi2_bbbb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_Zh_h_phi2phi2_bbbb_ATLAS13 constructor.
     */
    Hobs_Zh_h_phi2phi2_bbbb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{Zh\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to HH) \cdot BR^{\text{GTHDM}}(H\to bb) \cdot BR^{\text{GTHDM}}(H\to bb)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi2phi2_bbmumu_ATLAS13_old
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section times branching ratio of the process @f$pp\to h_{125}\to HH\to bb\mu\mu@f$ at 13 TeV.
 */
class Hobs_pp_h_phi2phi2_bbmumu_ATLAS13_old: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi2phi2_bbmumu_ATLAS13_old constructor.
     */
    Hobs_pp_h_phi2phi2_bbmumu_ATLAS13_old(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to HH) \cdot BR^{\text{GTHDM}}(H\to bb) \cdot BR^{\text{GTHDM}}(H\to \mu\mu)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_h_phi2phi2_gagagg_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section  times branching ratio of the process @f$pp \to h_{125}\to HH\to \gamma\gamma gg@f$.
 */
class Hobs_pp_h_phi2phi2_gagagg_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_h_phi2phi2_gagagg_ATLAS13 constructor.
     */
    Hobs_pp_h_phi2phi2_gagagg_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GTHDM}}_{pp\to h_{125}} \cdot BR^{\text{GTHDM}}(h_{125}\to HH) \cdot BR^{\text{GTHDM}}(H\to \gamma\gamma) \cdot BR^{\text{GTHDM}}(H\to gg)]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

#endif /* GENERALTHDMLOWMASS_H */