/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CPENGUINBOXMU_H
#define	CPENGUINBOXMU_H

#include <StandardModel.h>
#include <sstream>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_clausen.h>

class CPenguinBoxMu {
    
public:
    
    /**
     * 
     * @brief constructor
     */
    CPenguinBoxMu(const StandardModel& model_i);
    
    /**
     * 
     * @brief destructor
     */
    ~CPenguinBoxMu();
    
    /**
     * 
     * @brief hep-ph/9512380v1, page 92.
     * @return Z-penguin charm contribution to the decay K -> mu mubar. 
     */
    double C_NL();
    
     /**
     * 
     * @brief hep-ph/9512380v2, page 99.
     * @return EW box charm contribution to the decay K -> mu mubar, aleays NLO
     */
    double B_NL();
    
    /**
     * 
     * @return the charm contribution to to the decay K -> mu mubar, aleays NLO
     */
    double X_ch();
    
private:
    const StandardModel& model;
    const StandardModelMatching& modelmatching; 
    
};

#endif	/* CPENGUINBOXMU_H */
