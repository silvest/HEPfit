/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LIGHTHIGGS_H
#define	LIGHTHIGGS_H

#include "ThObservable.h"

class THDM;
class THDMcache;

/**
 * @class THDM_BR_h_bb
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$h\to b \bar b@f$.
 */
class THDM_BR_h_bb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    THDM_BR_h_bb(const StandardModel& SM_i);
    
    /**
     * @return @return @f$BR{\text THDM}(h\to b \bar b)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class THDM_BR_h_gaga
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$h\to \gamma \gamma@f$.
 */
class THDM_BR_h_gaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    THDM_BR_h_gaga(const StandardModel& SM_i);
    
    /**
     * @return @f$BR{\text THDM}(h\to \gamma \gamma)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class THDM_BR_h_tautau
 * @ingroup THDM
 * @brief %THDM branching ratio of @f$h\to \tau \tau@f$.
 */
class THDM_BR_h_tautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    THDM_BR_h_tautau(const StandardModel& SM_i);
    
    /**
     * @return @f$BR{\text THDM}(h\to \tau \tau)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htobb8
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two b quarks at 8 TeV.
 */
class ggF_tth_htobb8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htobb8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htoWW8
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$W@f$ bosons at 8 TeV.
 */
class ggF_tth_htoWW8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoWW8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to WW)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htotautau8
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$\tau@f$ leptons at 8 TeV.
 */
class ggF_tth_htotautau8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htotautau8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \tau\tau)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htoZZ8
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$Z@f$ bosons at 8 TeV.
 */
class ggF_tth_htoZZ8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoZZ8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htogaga8
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two photons at 8 TeV.
 */
class ggF_tth_htogaga8 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htogaga8(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \gamma\gamma)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htobb13
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two b quarks at 13 TeV.
 */
class ggF_tth_htobb13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htobb13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htoWW13
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$W@f$ bosons at 13 TeV.
 */
class ggF_tth_htoWW13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoWW13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to WW)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htotautau13
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$\tau@f$ leptons at 13 TeV.
 */
class ggF_tth_htotautau13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htotautau13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \tau\tau)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htoZZ13
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$Z@f$ bosons at 13 TeV.
 */
class ggF_tth_htoZZ13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoZZ13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htogaga13
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two photons at 13 TeV.
 */
class ggF_tth_htogaga13 : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htogaga13(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \gamma\gamma)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htobb
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$b@f$ quarks.
 */
class VBF_Vh_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htoWW
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$W@f$ bosons.
 */
class VBF_Vh_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to WW)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htotautau
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$\tau@f$ leptons.
 */
class VBF_Vh_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to \tau\tau)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htoZZ
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$Z@f$ bosons.
 */
class VBF_Vh_htoZZ : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htogaga
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two photons.
 */
class VBF_Vh_htogaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to \gamma\gamma)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htogg
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two gluons.
 */
class VBF_Vh_htogg : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htogg(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to gg)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htocc
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$c@f$ quarks.
 */
class VBF_Vh_htocc : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htocc(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to c\bar c)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_htobb
 * @ingroup THDM
 * @brief Signal strength of a ggF produced h decaying to two b quarks at 13 TeV.
 */
class ggF_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_htoWW
 * @ingroup THDM
 * @brief Signal strength of a ggF produced h decaying to two @f$W@f$ bosons at 13 TeV.
 */
class ggF_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF}(h\to WW)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_htotautau
 * @ingroup THDM
 * @brief Signal strength of a ggF produced h decaying to two @f$\tau@f$ leptons at 13 TeV.
 */
class ggF_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF}(h\to \tau \tau)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_htoZZ
 * @ingroup THDM
 * @brief Signal strength of a ggF produced h decaying to two @f$Z@f$ bosons at 13 TeV.
 */
class ggF_htoZZ : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_htogaga
 * @ingroup THDM
 * @brief Signal strength of a ggF produced h decaying to two photons at 13 TeV.
 */
class ggF_htogaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF}(h\to \gamma \gamma)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class tth_htobb
 * @ingroup THDM
 * @brief Signal strength of a tth produced h decaying to two b quarks at 13 TeV.
 */
class tth_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    tth_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class tth_htoWW
 * @ingroup THDM
 * @brief Signal strength of a tth produced h decaying to two @f$W@f$ bosons at 13 TeV.
 */
class tth_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    tth_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to WW)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class tth_htotautau
 * @ingroup THDM
 * @brief Signal strength of a tth produced h decaying to two @f$\tau@f$ leptons at 13 TeV.
 */
class tth_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    tth_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to \tau \tau)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class tth_htoZZ
 * @ingroup THDM
 * @brief Signal strength of a tth produced h decaying to two @f$Z@f$ bosons at 13 TeV.
 */
class tth_htoZZ : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    tth_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to ZZ)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class tth_htogaga
 * @ingroup THDM
 * @brief Signal strength of a tth produced h decaying to two photons at 13 TeV.
 */
class tth_htogaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    tth_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text tth}(h\to \gamma \gamma)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class mu_htobb
 * @ingroup THDM
 * @brief Signal strength of an h decaying to two b quarks at 13 TeV.
 */
class mu_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    mu_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text pp}(h\to b\bar b)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class mu_htoWW
 * @ingroup THDM
 * @brief Signal strength of an h decaying to two @f$W@f$ bosons at 13 TeV.
 */
class mu_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    mu_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text pp}(h\to WW)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class mu_htotautau
 * @ingroup THDM
 * @brief Signal strength of an h decaying to two @f$\tau@f$ leptons at 13 TeV.
 */
class mu_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    mu_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text pp}(h\to \tau\tau)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class mu_htoZga
 * @ingroup THDM
 * @brief Signal strength of an h decaying to a @f$Z@f$ boson and a photon at 13 TeV.
 */
class mu_htoZga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    mu_htoZga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text pp}(h\to Zga)@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};






/**
 * @class Gamma_h_THDM
 * @ingroup THDM
 * @brief Total h decay rate in the %THDM.
 */
class Gamma_h_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Gamma_h_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$\Gamma_h@f$ in units of GeV
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class rh_gaga_THDM
 * @ingroup THDM
 * @brief Squared relative coupling of @f$h@f$ to two photons.
 */
class rh_gaga_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    rh_gaga_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(h)}_{\gamma \gamma}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class rh_gg_THDM
 * @ingroup THDM
 * @brief Squared relative coupling of @f$h@f$ to two gluons.
 */
class rh_gg_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    rh_gg_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(h)}_{gg}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};

/**
 * @class rh_Zga_THDM
 * @ingroup THDM
 * @brief Squared relative coupling of @f$h@f$ to a @f$Z@f$ boson and a photon.
 */
class rh_Zga_THDM : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    rh_Zga_THDM(const StandardModel& SM_i);
    
    /**
     * @return @f$r^{(h)}_{Z\gamma}@f$
     */
    double computeThValue();
private:
    const THDM& myTHDM;
};


#endif	/* HIGGSSIGSTR_H */
