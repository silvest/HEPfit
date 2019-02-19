/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DIBOSONTHOBSERVABLES_H
#define	DIBOSONTHOBSERVABLES_H

#include "ThObservable.h"
#include "NPbase.h"

/**
 * @class mueeWW
 * @ingroup HiggsExtensions
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
    mueeWW(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWW called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to W^+ W^-}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to W^+ W^-}@f$
     */
    double computeThValue()
    {
        return myNPbase->mueeWW(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};

/**
 * @class mueeWWPol
 * @ingroup HiggsExtensions
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
    mueeWWPol(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mueeWWPol called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of @f$\mu_{e^+e^- \to W^+ W^-}@f$ in the current model.
     * @return @f$\mu_{e^+e^- \to W^+ W^-}@f$
     */
    double computeThValue()
    {
        return myNPbase->mueeWWPol(sqrt_s,Pol_em, Pol_ep);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s, Pol_em, Pol_ep;
};

#endif	/* DIBOSONTHOBSERVABLES_H */

