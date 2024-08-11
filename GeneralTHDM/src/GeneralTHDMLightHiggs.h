/* 
 * Copyright (C) 2018 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMLIGHTHIGGS_H
#define GENERALTHDMLIGHTHIGGS_H

#include "ThObservable.h"

class GeneralTHDM;
class GeneralTHDMcache;

/**
 * @class Gamma_h_GTHDM
 * @ingroup GeneralTHDM
 * @brief Total @f$h@f$ decay rate in the %GeneralTHDM.
 */
class Gamma_h_GTHDM: public ThObservable {
public:

    /**
     * @brief Gamma_h_GTHDM constructor.
     */
    Gamma_h_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\Gamma_{h,\text{GTHDM}}@f$ in units of GeV
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/////// IMPORTANT!!!
/////// So far these branching ratios don't seem to be used since we use the signal strenghts as observables
/////// It seems that the code could be optimised and at least use this definitions when computing the signal
/////// strengths, which is done in the GeneralTHDM.cpp. We'll leave it as it is for the moment though.

/**
 * @class BR_h_bb_GTHDM
 * @ingroup GTHDM
 * @brief %GTHDM branching ratio of @f$h\to b \bar b@f$.
 */
class BR_h_bb_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_h_bb_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @return @f$BR{\text GTHDM}(h\to b \bar b)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_h_gaga_GTHDM
 * @ingroup GTHDM
 * @brief %GTHDM branching ratio of @f$h\to \gamma \gamma@f$.
 */
class BR_h_gaga_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_h_gaga_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR{\text GTHDM}(h\to \gamma \gamma)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class BR_h_tautau_GTHDM
 * @ingroup GTHDM
 * @brief %GTHDM branching ratio of @f$h\to \tau \tau@f$.
 */
class BR_h_tautau_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    BR_h_tautau_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR{\text GTHDM}(h\to \tau \tau)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_tth_htobb8
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two b quarks at 8 TeV.
 */
class GTHDM_ggF_tth_htobb8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htobb8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class GTHDM_ggF_tth_htoWW8
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$W@f$ bosons at 8 TeV.
 */
class GTHDM_ggF_tth_htoWW8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htoWW8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to WW)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class GTHDM_ggF_tth_htotautau8
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$\tau@f$ leptons at 8 TeV.
 */
class GTHDM_ggF_tth_htotautau8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htotautau8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \tau\tau)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class GTHDM_ggF_tth_htoZZ8
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$Z@f$ bosons at 8 TeV.
 */
class GTHDM_ggF_tth_htoZZ8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htoZZ8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_tth_htogaga8
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two photons at 8 TeV.
 */
class GTHDM_ggF_tth_htogaga8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htogaga8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \gamma\gamma)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_tth_htobb13
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two b quarks at 13 TeV.
 */
class GTHDM_ggF_tth_htobb13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htobb13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_tth_htoWW13
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$W@f$ bosons at 13 TeV.
 */
class GTHDM_ggF_tth_htoWW13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htoWW13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to WW)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_tth_htotautau13
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$\tau@f$ leptons at 13 TeV.
 */
class GTHDM_ggF_tth_htotautau13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htotautau13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \tau\tau)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_tth_htoZZ13
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$Z@f$ bosons at 13 TeV.
 */
class GTHDM_ggF_tth_htoZZ13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htoZZ13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_tth_htogaga13
 * @ingroup GTHDM
 * @brief Signal strength of a ggF or tth produced h decaying to two photons at 13 TeV.
 */
class GTHDM_ggF_tth_htogaga13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_tth_htogaga13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \gamma\gamma)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_VBF_Vh_htobb
 * @ingroup GTHDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$b@f$ quarks.
 */
class GTHDM_VBF_Vh_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_VBF_Vh_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_VBF_Vh_htoWW
 * @ingroup GTHDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$W@f$ bosons.
 */
class GTHDM_VBF_Vh_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_VBF_Vh_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to WW)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_VBF_Vh_htotautau
 * @ingroup GTHDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$\tau@f$ leptons.
 */
class GTHDM_VBF_Vh_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_VBF_Vh_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to \tau\tau)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_VBF_Vh_htoZZ
 * @ingroup GTHDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$Z@f$ bosons.
 */
class GTHDM_VBF_Vh_htoZZ : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_VBF_Vh_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_VBF_Vh_htogaga
 * @ingroup GTHDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two photons.
 */
class GTHDM_VBF_Vh_htogaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_VBF_Vh_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to \gamma\gamma)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_VBF_Vh_htogg
 * @ingroup GTHDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two gluons.
 */
class GTHDM_VBF_Vh_htogg : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_VBF_Vh_htogg(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to gg)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_VBF_Vh_htocc
 * @ingroup GTHDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$c@f$ quarks.
 */
class GTHDM_VBF_Vh_htocc : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_VBF_Vh_htocc(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to c\bar c)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class ggF_htobb
 * @ingroup GTHDM
 * @brief Signal strength of a ggF produced h decaying to two b quarks at 13 TeV.
 */
class GTHDM_ggF_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_htoWW
 * @ingroup GTHDM
 * @brief Signal strength of a ggF produced h decaying to two @f$W@f$ bosons at 13 TeV.
 */
class GTHDM_ggF_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF}(h\to WW)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_htotautau
 * @ingroup GTHDM
 * @brief Signal strength of a ggF produced h decaying to two @f$\tau@f$ leptons at 13 TeV.
 */
class GTHDM_ggF_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF}(h\to \tau \tau)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_htoZZ
 * @ingroup GTHDM
 * @brief Signal strength of a ggF produced h decaying to two @f$Z@f$ bosons at 13 TeV.
 */
class GTHDM_ggF_htoZZ : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_ggF_htogaga
 * @ingroup GTHDM
 * @brief Signal strength of a ggF produced h decaying to two photons at 13 TeV.
 */
class GTHDM_ggF_htogaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_ggF_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text _ggF}(h\to \gamma \gamma)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_tth_htobb
 * @ingroup GTHDM
 * @brief Signal strength of a tth produced h decaying to two b quarks at 13 TeV.
 */
class GTHDM_tth_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_tth_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_tth_htoWW
 * @ingroup GTHDM
 * @brief Signal strength of a tth produced h decaying to two @f$W@f$ bosons at 13 TeV.
 */
class GTHDM_tth_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_tth_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to WW)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_tth_htotautau
 * @ingroup GTHDM
 * @brief Signal strength of a tth produced h decaying to two @f$\tau@f$ leptons at 13 TeV.
 */
class GTHDM_tth_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_tth_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to \tau \tau)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_tth_htoZZ
 * @ingroup GTHDM
 * @brief Signal strength of a tth produced h decaying to two @f$Z@f$ bosons at 13 TeV.
 */
class GTHDM_tth_htoZZ : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_tth_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_tth_htogaga
 * @ingroup GTHDM
 * @brief Signal strength of a tth produced h decaying to two photons at 13 TeV.
 */
class GTHDM_tth_htogaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_tth_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to \gamma \gamma)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_mu_htobb
 * @ingroup GTHDM
 * @brief Signal strength of an h decaying to two b quarks at 13 TeV.
 */
class GTHDM_mu_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_mu_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text pp}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_mu_htoWW
 * @ingroup GTHDM
 * @brief Signal strength of an h decaying to two @f$W@f$ bosons at 13 TeV.
 */
class GTHDM_mu_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_mu_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text pp}(h\to WW)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_mu_htotautau
 * @ingroup GTHDM
 * @brief Signal strength of an h decaying to two @f$\tau@f$ leptons at 13 TeV.
 */
class GTHDM_mu_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_mu_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text pp}(h\to \tau\tau)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class GTHDM_mu_htoZga
 * @ingroup GTHDM
 * @brief Signal strength of an h decaying to a @f$Z@f$ boson and a photon at 13 TeV.
 */
class GTHDM_mu_htoZga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    GTHDM_mu_htoZga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text pp}(h\to Zga)@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class yu1r_GTHDM
 * @ingroup GeneralTHDM
 * @brief  Coupling of the SM-Higgs to up quarks real part
 */
class yu1r_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    yu1r_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$yu1r_GTHDM@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
}; 

/**
 * @class yd1r_GTHDM
 * @ingroup GeneralTHDM
 * @brief  Coupling of the SM-Higgs to down quarks  real part
 */
class yd1r_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    yd1r_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$yd1r_GTHDM@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class yl1r_GTHDM
 * @ingroup GeneralTHDM
 * @brief  Coupling of the SM-Higgs to leptons  real part
 */
class yl1r_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    yl1r_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$yl1r_GTHDM@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class suR_GTHDM
 * @ingroup GeneralTHDM
 * @brief  Real part Yukawa up coupling
 */
class suR_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    suR_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$suR_GTHDM@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class sdR_GTHDM
 * @ingroup GeneralTHDM
 * @brief  Real part of Yukawa down coupling
 */
class sdR_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    sdR_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$sdR_GTHDM@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class slR_GTHDM
 * @ingroup GeneralTHDM
 * @brief  Real part of Yukawa leptonic coupling
 */
class slR_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    slR_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$slR_GTHDM@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class rh_gg_GTHDM
 * @ingroup GeneralTHDM
 * @brief Squared relative coupling of @f$h@f$ to two gluons.
 */
class rh_gg_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    rh_gg_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(h)}_{gg}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
}; 

    /**
 * @class rh_gaga_GTHDM
 * @ingroup GeneralTHDM
 * @brief Squared relative coupling of @f$h@f$ to two photons.
 */
class rh_gaga_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    rh_gaga_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(h)}_{\gamma \gamma}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
}; 
    
        /**
 * @class rh_Zga_GTHDM
 * @ingroup GeneralTHDM
 * @brief Squared relative coupling of @f$h@f$ to a Z and a photon.
 */
class rh_Zga_GTHDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    rh_Zga_GTHDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(h)}_{Z\gamma}@f$
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
}; 




#endif /* GENERALTHDMLIGHTHIGGS_H */
