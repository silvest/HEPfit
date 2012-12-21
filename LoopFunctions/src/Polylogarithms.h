/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef POLYLOGARITHMS_H
#define	POLYLOGARITHMS_H

#include <gslpp.h>
#include "BernoulliNumbers.h"
using namespace gslpp;


class Polylogarithms : public BernoulliNumbers {
public:

    /**
     * @brief Polylogarithms constructor
     */
    Polylogarithms();


    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] x real variable
     * @return dilogarithm Li_2(x)
     */
    complex Li2(const double x) const;    
    
    /**
     * @param[in] z complex variable
     * @return dilogarithm Li_2(z)
     */
    complex Li2(const complex z) const;
    
    /**
     * @param[in] x real variable
     * @return trilogarithm Li_3(x)
     * @attention applicable only for real x and x <= 1. 
     */
    double Li3(const double x) const;


    ////////////////////////////////////////////////////////////////////////

private:

  
    
};

#endif	/* POLYLOGARITHMS_H */

