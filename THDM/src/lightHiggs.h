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
#include "THDMfunctions.h"

/**
 * @class lightHiggs
 * @ingroup THDM
 * @brief .
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The squared relative @f$\gamma \gamma@f$ and @f$Z\gamma@f$ couplings are calculated at one-loop
 * following @cite Gunion:1989we.
 */
class lightHiggs : public ThObservable {
public:
    lightHiggs(const StandardModel& SM_i);
    virtual ~lightHiggs();
    void computeSignalStrengthQuantities();

    /**
     * @brief Empty function
     */
    double computeThValue();

    /**
     * @brief THDM branching ratio of @f$h\to b \bar b@f$.
     * @return @f$BR{\text THDM}(h\to b \bar b)@f$
     * @details This is also needed for the @f$H\to hh@f$ decay.
     */
    double THDM_BR_h_bb();

    /**
     * @brief THDM branching ratio of @f$h\to \gamma \gamma@f$.
     * @return @f$BR{\text THDM}(h\to \gamma \gamma)@f$
     * @details This is also needed for the @f$H\to hh@f$ decay.
     */
    double THDM_BR_h_gaga();

    /**
     * @brief THDM branching ratio of @f$h\to \tau \tau@f$.
     * @return @f$BR{\text THDM}(h\to \tau \tau)@f$
     * @details This is also needed for the @f$H\to hh@f$ decay.
     */
    double THDM_BR_h_tautau();

protected:

    /**
     * @brief SM branching ratio of @f$h\to b \bar b@f$.
     * @return @f$BR{\text SM}(h\to b \bar b)@f$
     */
    double BrSM_htobb;

    /**
     * @brief SM branching ratio of @f$h\to \gamma \gamma@f$.
     * @return @f$BR{\text SM}(h\to \gamma \gamma)@f$
     */
    double BrSM_htogaga;

    /**
     * @brief SM branching ratio of @f$h\to \tau \tau@f$.
     * @return @f$BR{\text SM}(h\to \tau \tau)@f$
     */
    double BrSM_htotautau;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two down type quarks.
     * @return @f$r^{(h)}_{Q_dQ_d}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_QdQd;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two massive vector bosons.
     * @return @f$r^{(h)}_{WW}=r^{(h)}_{ZZ}@f$
     */
    double rh_VV;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two charged leptons.
     * @return @f$r^{(h)}_{\ell \ell}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_ll;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two photons.
     * @return @f$r^{(h)}_{\gamma \gamma}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_gaga;

    /**
     * @brief Ratio of THDM and SM cross sections for ggF and tth production of h.
     * @return @f$\sigma^{\text THDM}_{\text ggF+tth}/\sigma^{\text SM}_{\text ggF+tth}@f$
     */
    double ggF_tth;

    /**
     * @brief Ratio of THDM and SM cross sections for VBF and Vh production of h.
     * @return @f$\sigma^{\text THDM}_{\text VBF+Vh}/\sigma^{\text SM}_{\text VBF+Vh}@f$
     */
    double VBF_Vh;

    /**
     * @brief Sum of the modified branching ratios.
     * @return @f$\sum _i r^{(h)}_{i} BR^{\text SM}(h\to i)@f$
     */
    double sumModBRs;

private:
    const THDM * myTHDM;
    const StandardModel& mySM;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two gluons.
     * @return @f$r^{(h)}_{gg}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_gg;

    /**
     * @brief Squared relative coupling of @f$h@f$ to a @f$Z@f$ boson and a photon.
     * @return @f$r^{(h)}_{Z\gamma}@f$
     * @details Depends on the type of @f$Z_2@f$ symmetry.
     */
    double rh_Zga;

    /**
     * @brief Squared relative coupling of @f$h@f$ to two up type quarks.
     * @return @f$r^{(h)}_{Q_uQ_u}@f$
     */
    double rh_QuQu;
};

/**
 * @class ggF_tth_htobb
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two b quarks.
 */
class ggF_tth_htobb : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to b\bar b)@f$
     */
    double computeThValue ();
};

/**
 * @class ggF_tth_htoWW
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$W@f$ bosons.
 */
class ggF_tth_htoWW : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to WW)@f$
     */
    double computeThValue ();
};

/**
 * @class ggF_tth_htotautau
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$\tau@f$ leptons.
 */
class ggF_tth_htotautau : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \tau\tau)@f$
     */
    double computeThValue ();
};

/**
 * @class ggF_tth_htoZZ
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two @f$Z@f$ bosons.
 */
class ggF_tth_htoZZ : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to ZZ)@f$
     */
    double computeThValue ();
};

/**
 * @class ggF_tth_htogaga
 * @ingroup THDM
 * @brief Signal strength of a ggF or tth produced h decaying to two photons.
 */
class ggF_tth_htogaga : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    ggF_tth_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text ggF+tth}(h\to \gamma\gamma)@f$
     */
    double computeThValue ();
};

/**
 * @class VBF_Vh_htobb
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two b quarks.
 */
class VBF_Vh_htobb : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htobb(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to b\bar b)@f$
     */
    double computeThValue ();
};

/**
 * @class VBF_Vh_htoWW
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$W@f$ bosons.
 */
class VBF_Vh_htoWW : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htoWW(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to WW)@f$
     */
    double computeThValue ();
};

/**
 * @class VBF_Vh_htotautau
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$\tau@f$ leptons.
 */
class VBF_Vh_htotautau : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htotautau(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to \tau\tau)@f$
     */
    double computeThValue ();
};

/**
 * @class VBF_Vh_htoZZ
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two @f$Z@f$ bosons.
 */
class VBF_Vh_htoZZ : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htoZZ(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to ZZ)@f$
     */
    double computeThValue ();
};

/**
 * @class VBF_Vh_htogaga
 * @ingroup THDM
 * @brief Signal strength of a VBF or Vh produced h decaying to two photons.
 */
class VBF_Vh_htogaga : public lightHiggs {
public:
    
    /**
     * @brief Constructor.
     */
    VBF_Vh_htogaga(const StandardModel& SM_i);
    
    /**
     * @return @f$\mu_{\text VBF+Vh}(h\to \gamma\gamma)@f$
     */
    double computeThValue ();
};


#endif	/* HIGGSSIGSTR_H */
