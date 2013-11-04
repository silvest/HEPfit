/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
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
 * @brief A class for polylogarithms. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class Polylogarithms : public BernoulliNumbers {
public:

    /**
     * @brief Polylogarithms constructor. 
     */
    Polylogarithms();

    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Dilogarithm function. 
     * @param[in] x Real variable. 
     * @return Dilogarithm Li_2(x).
     */
    complex Li2(const double x) const;    
    
    /**
     * @brief Dilogarithm function.
     * @param[in] z Complex variable.
     * @return Dilogarithm Li_2(z).
     */
    complex Li2(const complex z) const;
    
    /**
     * @brief Trilogarithm function.
     * @param[in] x Real variable.
     * @return Trilogarithm Li_3(x). 
     * @attention applicable only for real x and x <= 1. 
     */
    double Li3(const double x) const;

    ////////////////////////////////////////////////////////////////////////

private:
    
};

#endif	/* POLYLOGARITHMS_H */

