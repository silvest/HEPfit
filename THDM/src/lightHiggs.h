/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LIGHTHIGGS_H
#define	LIGHTHIGGS_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"
#include "THDMcache.h"

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
    double computeThValue ();
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
    double computeThValue ();
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
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htobb
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two b quarks.
 */
class ggF_tth_htobb : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to b\bar b)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htoWW
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$W@f$ bosons.
 */
class ggF_tth_htoWW : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to WW)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htotautau
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$\tau@f$ leptons.
 */
class ggF_tth_htotautau : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \tau\tau)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htoZZ
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$Z@f$ bosons.
 */
class ggF_tth_htoZZ : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to ZZ)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class ggF_tth_htogaga
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two photons.
 */
class ggF_tth_htogaga : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \gamma\gamma)@f$
     */
    double computeThValue ();
private:
    const THDM& myTHDM;
};

/**
 * @class VBF_Vh_htobb
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two b quarks.
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
    double computeThValue ();
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
    double computeThValue ();
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
    double computeThValue ();
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
    double computeThValue ();
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
    double computeThValue ();
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
    double computeThValue ();
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
    double computeThValue ();
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
    double computeThValue ();
private:
    const THDM& myTHDM;
};


#endif	/* HIGGSSIGSTR_H */
