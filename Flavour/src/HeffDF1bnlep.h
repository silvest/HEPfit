/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDF1BNLEP_H
#define	HEFFDF1BNLEP_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>
#include "EvolDF1nlep.h"
#include <sstream>

class HeffDF1bnlep {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDF1bnlep(const StandardModel & SM, StandardModelMatching& modelmatching);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDF1bnlep();
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the reonrmalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 0 | \f$, \f$ |\Delta S = 0 | \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBnlep00(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 1 | \f$, \f$ |\Delta S = 0 | \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBnlep10(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 0 | \f$, \f$ |\Delta S = 1 | \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBnlep01(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 1 | \f$, \f$ |\Delta S = 1 | \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBnlep11(double mu, schemes scheme = NDR);
    
    WilsonCoefficient getCoeffbnlep00() const {
        return coeffbnlep00;
    }
    
    WilsonCoefficient getCoeffbnlep10() const {
        return coeffbnlep01;
    }
    
    WilsonCoefficient getCoeffbnlep01() const {
        return coeffbnlep10;
    }
    
    WilsonCoefficient getCoeffbnlep11() const {
        return coeffbnlep11;
    }
    
    EvolDF1nlep getUDF1() const {
        return u;
    }

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;
    
    WilsonCoefficient coeffbnlep00qcd, coeffbnlep00;
    WilsonCoefficient coeffbnlep10qcd, coeffbnlep10;
    WilsonCoefficient coeffbnlep01, coeffbnlep01A, coeffbnlep01B, coeffbnlep00CC;
    WilsonCoefficient coeffbnlep11, coeffbnlep11A, coeffbnlep11B, coeffbnlep10CC;
    EvolDF1nlep u;
    
    gslpp::vector<gslpp::complex> bnlep, bnlep2, bnlepCC;
};

#endif	/* HEFFDF1BNLEP_H */


