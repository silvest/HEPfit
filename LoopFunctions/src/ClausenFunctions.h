/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CLAUSENFUNCTIONS_H
#define	CLAUSENFUNCTIONS_H

#include "BernoulliNumbers.h"

/**
 * @class ClausenFunctions
 * @ingroup LoopFunctions 
 * @brief A class for the Clausen functions.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class handles the Clausen functions
 * @f$\mathrm{Cl}_2(\phi)@f$ and @f$\mathrm{Cl}_3(\phi)@f$. 
 */
class ClausenFunctions : public BernoulliNumbers {
public:
    
    /**
     * @brief The default constructor.
     */
    ClausenFunctions();
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The Clausen function of index 2, @f$\mathrm{Cl}_2(\phi)@f$.
     * @details This function calls the GSL function gsl_sf_clausen().
     * @param[in] phi a real variable
     * @return @f$\mathrm{Cl}_2(\phi)@f$
     */
    double Cl2(const double phi) const;    
    
    /**
     * @brief The Clausen function of index 3, @f$\mathrm{Cl}_3(\phi)@f$.
     * @details The function @f$\mathrm{Cl}_3(\phi)@f$ is computed with the help
     * of the GSL function gsl_sf_zeta_int(3). See @cite Kniehl:1989qu. 
     * @param[in] phi a real variable
     * @return @f$\mathrm{Cl}_3(\phi) = \zeta(3) - \phi^2 \left( \frac{3}{4} - \frac{\ln(\phi)}{2} + \sum_{n=1}^{\infty} \frac{B_n \phi^{2n}}{2n(2n+1)(2n+2)(2n)!} \right)@f$
     * @note The function is even, so the argument is internally reduced to the range [0, pi]. The sum is truncated at n=18.
     * @note The Bernoulli numbers are stored in the base class BernoulliNumbers.
     */
    double Cl3(const double phi) const;

};

#endif	/* CLAUSENFUNCTIONS_H */

