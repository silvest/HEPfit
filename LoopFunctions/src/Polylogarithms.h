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
    
    /**
     * @brief Polylogarithms copy constructor
     * @param[in] orig reference to a Polylogarithms object
     */
    //Polylogarithms(const Polylogarithms& orig);
    
    /**
     * @brief Polylogarithms destructor
     */
    virtual ~Polylogarithms();


    ////////////////////////////////////////////////////////////////////////

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

