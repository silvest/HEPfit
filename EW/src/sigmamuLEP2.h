/* 
 * File:   sigmamuLEP2.h
 * Author: giovannigrilli
 *
 * Created on 30 agosto 2012, 11.40
 */

#ifndef SIGMAMULEP2_H
#define	SIGMAMULEP2_H

#include <ThObservable.h>
#include "EW.h"
#include <iostream>
#include <cmath>
#include <cuba.h>


//static int sqrt_s_global;


class sigmamuLEP2 : public ThObservable {
public:
    /**
     * @brief sigmamuLEP2 constructor
     * @param[in] EW_i an object of EW class
     */
    sigmamuLEP2(const EW& EW_i,double sqrt_s_i) : ThObservable(EW_i), 
            myEW(EW_i), sqrt_s_local(sqrt_s_i){};
    
    
    double sigmaLEP2_Integrand(double x);
    
    static int Integrand(const int *ndim, const double theta[],
                     const int *ncomp, double dsigma[12], void *userdata);
    
    /**
     * @return the muon cross section for LEP2 energies 
     */
    double getThValue();
   
    //static const double sqrt_s[12]; 


private:
    const EW& myEW;
    const double sqrt_s_local;

};

#endif	/* SIGMAMULEP2_H */

