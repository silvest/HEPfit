/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FUNCTIONS_H
#define	FUNCTIONS_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class THDMfunctions
 * @ingroup THDM 
 * @brief 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class THDMfunctions : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    THDMfunctions(const StandardModel& SM_i) 
    : ThObservable(SM_i), myTHDM(static_cast<const THDM*> (&SM_i))
    {
    };

    /**
     * @brief The quartic coupling @f$\lambda_5@f$.
     * @return @f$\lambda_5@f$
     */
    double computeThValue();
    
    gslpp::complex f_func(const double x)const;
    gslpp::complex g_func(const double x)const;
    gslpp::complex Int1(const double tau, const double lambda)const;
    gslpp::complex Int2(const double tau, const double lambda)const;
    int HSTheta (const double x) const;
    double KaellenFunction (const double a, const double b, const double c) const;

    private:
        const THDM * myTHDM;
};

#endif	/* FUNCTIONS_H */





