/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PP_DIBOSON_H
#define PP_DIBOSON_H

#include <stdexcept>
#include <ThObservable.h>
#include "NPbase.h"

/**
 * @addtogroup EW
 * @brief A module for electroweak precision observables.
 * @details 
 * @{
 */


/**
 * @class ppZHprobe
 * @brief A class for implementing the direction constrained by @f$ p p \to Z H@f$ in the boosted regime, @f$g_p^Z@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the cross section @f$e^+ e^- \to W^+ W^-@f$.
 */
class ppZHprobe : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    ppZHprobe(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("ppZHprobe called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the direction constrained by @f$ p p \to Z H@f$ in the boosted regime in the current model.
     * @return @f$g_p^Z@f$
     */
    double computeThValue()
    {
        return myNPbase->ppZHprobe(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class NpTVppWZ
 * @brief A class for computing the number of events in  @f$ p p \to WZ@f$
 * in a given @f$p_{TV}@f$ bin, normalized to the SM prediction.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the number of events in  @f$ p p \to WZ@f$
 * in a given @f$p_{TV}@f$ bin, normalized to the SM prediction.
 */
class mupTVppWZ : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    mupTVppWZ(const StandardModel& SM_i, const double sqrt_s_i, const double pTV1_i, const double pTV2_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), pTV1(pTV1_i), pTV2(pTV2_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("mupTVppWZ called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the number of events in  @f$ p p \to WZ@f$
     * in a given @f$p_{TV}@f$ bin, normalized to the SM prediction, in the current model.
     * @return @f$N_{ev}^{p_{TV}}/N_{ev,SM}^{p_{TV}}@f$
     */
    double computeThValue()
    {
        return myNPbase->mupTVppWZ(sqrt_s, pTV1, pTV2);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    const double pTV1, pTV2;
};

/** 
 * @}
 */

#endif	/* PP_DIBOSON_H */

