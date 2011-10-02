/* 
 * File:   Polylogarithms.h
 * Author: mishima
 */

#ifndef POLYLOGARITHMS_H
#define	POLYLOGARITHMS_H

#include <gsl/gsl_math.h>
#include <gslpp.h>
using namespace gslpp;


class Polylogarithms {
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
     * @param[in] k the value for the argument x
     * @param[in] dim dim=1
     * @param[in] params this parameter is unnecessary. 
     * @return Li_2(x)/x
     */
    static double integrand_for_Li3(double *k, size_t dim, void *params);
    
    /**
     * @param[in] x 
     * @return Polylogarithm of order three Li_3(x)
     * @attention valid only for x <= 1. 
     */
    double Li3(const double x) const;
    
    
    
    

    ////////////////////////////////////////////////////////////////////////

private:
    double B0[19]; /* Bernoulli numbers */
  
    
    
};

#endif	/* POLYLOGARITHMS_H */

