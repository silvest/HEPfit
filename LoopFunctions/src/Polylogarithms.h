/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef POLYLOGARITHMS_H
#define	POLYLOGARITHMS_H

#include <gslpp.h>
#include "BernoulliNumbers.h"
using namespace gslpp;

/**
 * @class Polylogarithms
 * @ingroup LoopFunctions 
 * @brief A class for the polylogarithms.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class handles the polylogarithms 
 * @f$\mathrm{Li}_2(x)@f$, @f$\mathrm{Li}_2(z)@f$ and @f$\mathrm{Li}_3(x)@f$,
 * where @f$x@f$ and @f$z@f$ are real and complex variables, respectively. 
 */
class Polylogarithms : public BernoulliNumbers {
public:

    /**
     * @brief The default constructor.
     */
    Polylogarithms();

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The dilogarithm with a real argument, @f$\mathrm{Li}_2(x)@f$.
     * @details This function calls the GSL function gsl_sf_complex_dilog_xy_e(). 
     * @param[in] x a real variable.
     * @return @f$\mathrm{Li}_2(x)@f$
     */
    complex Li2(const double x) const;    
    
    /**
     * @brief The dilogarithm with a complex argument, @f$\mathrm{Li}_2(z)@f$.
     * @details This function calls the GSL function gsl_sf_complex_dilog_xy_e(). 
     * @param[in] z a complex variable.
     * @return @f$\mathrm{Li}_2(z)@f$
     */
    complex Li2(const complex z) const;
    
    /**
     * @brief The trilogarithm @f$\mathrm{Li}_3(x)@f$
     * @param[in] x a real variable.
     * @return @f$\mathrm{Li}_3(x)@f$
     *
     * @attention This function is applicable for real @f$x \leq 1@f$.
     */
    double Li3(const double x) const;
    
};

#endif	/* POLYLOGARITHMS_H */

