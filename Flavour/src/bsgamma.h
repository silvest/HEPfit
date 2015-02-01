/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BSGAMMA_H
#define	BSGAMMA_H

#include <ThObservable.h>
#include "Flavour.h"
#include <StandardModel.h>

/**
 * @class Bsgamma
 * @ingroup flavour
 * @brief A class for the decay b -> s gamma. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Bsgamma : public ThObservable {
public:   
    
    
    Bsgamma(const StandardModel& SM_i, int obsFlag);
    
    double GF;
    double ale;
    double Mb;
    double mu_b;
    double width;
    gslpp::complex lambda_t;
    double coeff;
    gslpp::vector<gslpp::complex> ** allcoeff;
    gslpp::complex C_2;
    gslpp::complex C_7;
    gslpp::complex C_8;
    int obs;
    
    double Gamma;
    double Gamma_conj;
    
    void computeCoeff(orders order);
    void computeGamma(orders order);
    double computeThValue();
    
    
private:
    
    const StandardModel& mySM;
};

#endif	/* BSGAMMA_H */