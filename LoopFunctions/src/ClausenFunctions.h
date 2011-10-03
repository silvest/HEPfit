/* 
 * File:   ClausenFunctions.h
 * Author: mishima
 */

#ifndef CLAUSENFUNCTIONS_H
#define	CLAUSENFUNCTIONS_H

#include <gsl/gsl_math.h>
#include <gslpp.h>
using namespace gslpp;


class ClausenFunctions {
public:
    
    /**
     * @brief  ClausenFunctions constructor
     */
    ClausenFunctions();
    
    /**
     * @brief ClausenFunctions copy constructor
     * @param[in] orig reference to a ClausenFunctions object
     */    
    //ClausenFunctions(const ClausenFunctions& orig);

    /**
     * @brief ClausenFunctions destructor
     */
    virtual ~ClausenFunctions();

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] phi
     * @return Clausen function of index 2, Cl_2(phi)
     */
    double Cl2(const double phi) const;    
    
    /**
     * @param[in] phi
     * @return Clausen function of index 3, Cl_3(phi)
     * @attention phi=0.0 is not allowed. 
     */
    double Cl3(const double phi) const;
    
    
    ////////////////////////////////////////////////////////////////////////    
    
private:

    struct my_f_params { double phi; };
    
    static double integrand_for_Li3_imag(double *k, size_t dim, void *phi);
    
    
};

#endif	/* CLAUSENFUNCTIONS_H */

