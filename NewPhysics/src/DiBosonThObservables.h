/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DIBOSONTHOBSERVABLES_H
#define	DIBOSONTHOBSERVABLES_H

#include "ThObservable.h"

class NPbase;

/**
 * @addtogroup NewPhysics
 * @brief A module for model-independent studies of new physics.
 * @{
 */

/**
 * @class mueeWW
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to W^+ W^- }@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to W^+ W^-}@f$ between the 
 * @f$e^+e^- \to W^+ W^-@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV
     */
    mueeWW(const StandardModel& SM_i, const double sqrt_s_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to W^+ W^-}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to W^+ W^-}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeWWPol
 * @ingroup NewPhysics
 * @brief A class for computing the ratio @f$\mu_{e^+e^- \to W^+ W^-}@f$.
 * @author HEPfit CollaborationH
 * @copyright GNU General Public License
 * @details A class for computing the ratio @f$\mu_{e^+e^- \to W^+ W^-}@f$ between the 
 * @f$e^+e^- \to W^+ W^-@f$ 
 * associated production cross-section in the current model and in the Standard Model.
 */
class mueeWWPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in TeV, Pol_em_i and Pol_ep_i
     * are the polarization of electrons and positrons, respectively
     */
    mueeWWPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i);

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to W^+ W^-}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to W^+ W^-}@f$
     */
    double computeThValue();

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

/**
 * @}
 */

#endif	/* DIBOSONTHOBSERVABLES_H */

