/* 
 * File:   Resummations.cpp
 * Author: mishima
 */

#include "Resummations.h"


Resummations::Resummations() {
}

//Resummations::Resummations(const Resummations& orig) {
//}

Resummations::~Resummations() {
}


////////////////////////////////////////////////////////////////////////

double Resummations::ONEplusDeltaR(const schemes_EW schemeMw, 
                                   const double cW2_to_sW2, 
                                   const double DeltaAlpha_l5q, 
                                   const double DeltaRho, 
                                   const double DeltaR_rem) {
    double X;

    switch (schemeMw) {
        case NORESUM: 
            X = DeltaAlpha_l5q - cW2_to_sW2*DeltaRho + DeltaR_rem;
            break;
        case OMSI:
            throw "Write codes for resummationMw()";
            break;
        case INTERMEDIATE:
            throw "Write codes for resummationMw()";            
            break;        
        case OMSII:
            throw "Write codes for resummationMw()";            
            break;
        default:
            throw "Error in resummationMw()";            
            break;
    }   
    
    return X;
}

complex Resummations::rhoZ(const schemes_EW schemeRhoZ, 
                           const double cW2_to_sW2, const double DeltaRho,
                           const complex deltaRho_rem) {
    complex rhoZ;
    
    /* Real parts */
    switch (schemeRhoZ) {
        case NORESUM: 
            rhoZ.real() = 1.0 + DeltaRho + deltaRho_rem.real();
            break;
        case OMSI:
            throw "Write codes for resummationRhoZ()";
            break;
        case INTERMEDIATE:
            throw "Write codes for resummationRhoZ()";
            break;        
        case OMSII:
            throw "Write codes for resummationRhoZ()";
            break;
        default:
            throw "Error in resummationRhoZ()";
            break;
    }

    /* Imaginary parts */
    rhoZ.imag() = deltaRho_rem.imag();
    
    return rhoZ;    
}

complex Resummations::kappaZ(const schemes_EW schemeKappaZ, 
                             const double cW2_to_sW2, const double DeltaRho, 
                             const complex deltaKappa_rem) {
    complex kappaZ;
    
    /* Real parts */
    switch (schemeKappaZ) {
        case NORESUM: 
            kappaZ.real() = 1.0 + cW2_to_sW2*DeltaRho + deltaKappa_rem.real();  
            break;
        case OMSI:
            throw "Write codes for resummationKappaZ()";
            break;
        case INTERMEDIATE:
            throw "Write codes for resummationKappaZ()";
            break;        
        case OMSII:
            throw "Write codes for resummationKappaZ()";
            break;
        case APPROXIMATEFORMULA:
            /* The real parts are given by the approximate formulae. 
             * See ComputeKappaZ() */
            break;
        default:
            throw "Error in resummationKappaZ()";
            break;
    }

    /* Imaginary parts */
    kappaZ.imag() = deltaKappa_rem.imag();

    return kappaZ;
}

