/* 
 * File:   Resummations.h
 * Author: mishima
 */

#ifndef RESUMMATIONS_H
#define	RESUMMATIONS_H

#include <gslpp.h>
#include "Schemes.h"

using namespace gslpp;


class Resummations {
public:

    /**
     * @brief Resummations constructor
     */
    Resummations();
    
    /**
     * @brief Resummations copy constructor
     * @param[in] orig reference to a Resummations object
     */
    //Resummations(const Resummations& orig);
    
    /**
     * @brief Resummations destructor
     */
    virtual ~Resummations();

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] schemeMw
     * @param[in] cW2_to_sW2
     * @param[in] DeltaAlpha_l5q
     * @param[in] DeltaRho
     * @param[in] DeltaR_rem
     * @return resummed 1+DeltaR ==> 1/(1-DeltaR)
     */
    double ONEplusDeltaR(const schemes_EW schemeMw, 
                         const double cW2_to_sW2, const double DeltaAlpha_l5q, 
                         const double DeltaRho, const double DeltaR_rem);
    
    /**
     * @param[in] schemeRhoZ
     * @param[in] cW2_to_sW2
     * @param[in] DeltaRho
     * @param[in] deltaRho_rem
     * @return resummed rho_Z^f
     */
    complex rhoZ(const schemes_EW schemeRhoZ, 
                 const double cW2_to_sW2, const double DeltaRho,
                 const complex deltaRho_rem);        
    
    /**
     * @param[in] schemeKappaZ
     * @param[in] cW2_to_sW2
     * @param[in] DeltaRho
     * @param[in] deltaKappa_rem
     * @return resummed kappa_Z^f
     */
    complex kappaZ(const schemes_EW schemeKappaZ, 
                   const double cW2_to_sW2, const double DeltaRho, 
                   const complex deltaKappa_rem);     
    
    
    ////////////////////////////////////////////////////////////////////////
    
private:

};

#endif	/* RESUMMATIONS_H */

