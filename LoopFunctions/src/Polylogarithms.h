/* 
 * File:   Polylogarithms.h
 * Author: mishima
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
     * @param[in] z
     * @return dilogarithm Li_2(z)
     */
    complex Li2(const complex z) const;
    
    /**
     * @param[in] x 
     * @return trilogarithm Li_3(x)
     * @attention applicable only for real x and x <= 1. 
     */
    double Li3(const double x) const;


    ////////////////////////////////////////////////////////////////////////

private:

  
    
};

#endif	/* POLYLOGARITHMS_H */

