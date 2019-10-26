/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EEWW_H
#define EEWW_H

#include <stdexcept>
#include <ThObservable.h>
#include "NPbase.h"


/**
 * @class xseeWW
 * @brief A class for computing the cross section @f$e^+ e^- \to W^+ W^-@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the cross section @f$e^+ e^- \to W^+ W^-@f$.
 */
class xseeWW : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    xseeWW(const StandardModel& SM_i, const double sqrt_s_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("xseeWW called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the value of the cross section @f$e^+ e^- \to W^+ W^-@f$ in the current model.
     * @return @f$\sigma(e^+ e^- \to W^+ W^-)@f$
     */
    double computeThValue()
    {
        return myNPbase->xseeWW(sqrt_s);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
};


/**
 * @class dxseeWWdcosBin
 * @brief A class for computing the integral of the differential cross section 
 * for @f$e^+ e^- \to W^+ W^-@f$ in a given @f$\cos{\theta}@f$ bin.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details A class for computing the integral of the differential cross section
 * for @f$e^+ e^- \to W^+ W^-@f$ in a given @f$\cos{\theta}@f$ bin..
 */
class dxseeWWdcosBin : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to a StandardModel object or to any extension of it
     * @param[in] sqrt_s_i the center-of-mass energy in GeV
     */
    dxseeWWdcosBin(const StandardModel& SM_i, const double sqrt_s_i, const double cos1_i, const double cos2_i)
    : ThObservable(SM_i), sqrt_s(sqrt_s_i), cos1(cos1_i), cos2(cos2_i)
    {
        if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
            throw std::runtime_error("dxseeWWdcosBin called with a class whose parent is not NPbase");
    }

    /**
     * @brief A method to compute the integral of the differential cross section 
     * for @f$e^+ e^- \to W^+ W^-@f$ in a given @f$\cos{\theta}@f$ bin in the current model.
     * @return @f$\int_{\cos{\theta_1}}^{\cos{\theta_2}} d\sigma(e^+ e^- \to W^+ W^-)/d\cos{\theta}@f$
     */
    double computeThValue()
    {
        return myNPbase->dxseeWWdcosBin(sqrt_s, cos1, cos2);
    }

private:
    const NPbase* myNPbase;
    const double sqrt_s;
    const double cos1, cos2;
};

#endif	/* EEWW_H */

