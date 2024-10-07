/*
 * Copyright (C) 2024 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMSUSYOBS_H
#define GENERALTHDMSUSYOBS_H

#include "ThObservable.h"

class GeneralTHDM;
class GeneralTHDMcache;

/**
 * @class GeneralTHDMSusyObs
 * @ingroup GeneralTHDM 
 * @brief Base class for SUSY observables that are interpretable as THDM constraints.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theory over experiment ratios of
 * several phenomena containing sparticles which can also describe searches for
 * generic extra scalars within a multi-Higgs model context.
 */


/**
 * @class Hobs_pp_HpHm_taunutaunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross-section times branching ratios
 * associated with the process @f$p p\to H^{+} H^{-} \to \tau^{+} \nu_{\tau} \tau^{-} \bar\nu_{\tau}@f$,
 * as interpreted from a search for pair-produced staus decaying into taus and neutralinos
 */
class Hobs_pp_HpHm_taunutaunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_HpHm_taunutaunu_ATLAS13 constructor.
     */
    Hobs_pp_HpHm_taunutaunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text GTHDM}(p p \to H^+ H^-) \cdot BR^{\text{GTHDM}}(H^{+} \to \tau^{+} \nu_{\tau}) \cdot BR^{\text{GTHDM}}(H^{-} \to \tau^{-} \bar\nu_{\tau})]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_HpHm_taunutaunu_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross-section times branching ratios
 * associated with the process @f$p p\to H^{+} H^{-} \to \tau^{+} \nu_{\tau} \tau^{-} \bar\nu_{\tau}@f$,
 * as interpreted from a search for pair-produced staus decaying into taus and neutralinos
 */
class Hobs_pp_HpHm_taunutaunu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_HpHm_taunutaunu_CMS13 constructor.
     */
    Hobs_pp_HpHm_taunutaunu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text GTHDM}(p p \to H^+ H^-) \cdot BR^{\text{GTHDM}}(H^{+} \to \tau^{+} \nu_{\tau}) \cdot BR^{\text{GTHDM}}(H^{-} \to \tau^{-} \bar\nu_{\tau})]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_HpHm_munumunu_ATLAS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and ATLAS upper limit for the cross-section times branching ratios
 * associated with the process @f$p p\to H^{+} H^{-} \to \mu^{+} \nu_{\mu} \mu^{-} \bar\nu_{\mu}@f$,
 * as interpreted from a search for pair-produced sleptons decaying into charged leptons and neutralinos
 */
class Hobs_pp_HpHm_munumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_HpHm_munumunu_ATLAS13 constructor.
     */
    Hobs_pp_HpHm_munumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text GTHDM}(p p \to H^+ H^-) \cdot BR^{\text{GTHDM}}(H^{+} \to \mu^{+} \nu_{\mu}) \cdot BR^{\text{GTHDM}}(H^{-} \to \mu^{-} \bar\nu_{\mu})]_{\frac{\text{theo}}{\text{ATLAS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_pp_HpHm_munumunu_CMS13
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and CMS upper limit for the cross-section times branching ratios
 * associated with the process @f$p p\to H^{+} H^{-} \to \mu^{+} \nu_{\mu} \mu^{-} \bar\nu_{\mu}@f$,
 * as interpreted from a search for pair-produced sleptons decaying into charged leptons and neutralinos
 */
class Hobs_pp_HpHm_munumunu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_HpHm_munumunu_CMS13 constructor.
     */
    Hobs_pp_HpHm_munumunu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text GTHDM}(p p \to H^+ H^-) \cdot BR^{\text{GTHDM}}(H^{+} \to \mu^{+} \nu_{\mu}) \cdot BR^{\text{GTHDM}}(H^{-} \to \mu^{-} \bar\nu_{\mu})]_{\frac{\text{theo}}{\text{CMS,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Hobs_HpHm_munumunu_LEP208
 * @ingroup GeneralTHDM
 * @brief Ratio of the prediction and combined LEP upper limit for the cross-section times branching ratios
 * associated with the process @f$e^{+} e^{-}\to H^{+} H^{-} \to \mu^{+} \nu_{\mu} \mu^{-} \bar\nu_{\mu}@f$,
 * as interpreted from a search for pair-produced sleptons decaying into charged leptons and neutralinos
 */
class Hobs_HpHm_munumunu_LEP208: public ThObservable {
public:

    /**
     * @brief Hobs_HpHm_munumunu_LEP208 constructor.
     */
    Hobs_HpHm_munumunu_LEP208(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text GTHDM}(e^+ e^- \to H^+ H^-) \cdot BR^{\text{GTHDM}}(H^{+} \to \mu^{+} \nu_{\mu}) \cdot BR^{\text{GTHDM}}(H^{-} \to \mu^{-} \bar\nu_{\mu})]_{\frac{\text{theo}}{\text{LEP,95\%}}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};


#endif /* GENERALTHDMSUSYOBS_H */
