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

class sigmamuLEP2 : public ThObservable {
public:
    /**
     * @brief sigmamuLEP2 constructor
     * @param[in] EW_i an object of EW class
     */
    sigmamuLEP2(const EW& EW_i,const double sqrt_s): ThObservable(EW_i), myEW(EW_i){};
    
    /**
     * @return the muoniv cross section for LEP2 energies 
     */
    double getThValue();
    
    static int Integrand(const int *ndim, const double theta,const double sqrt_s,
                     const int *ncomp, double dsigma[12], void *userdata);


private:
    EW& myEW;

};

#endif	/* SIGMAMULEP2_H */

