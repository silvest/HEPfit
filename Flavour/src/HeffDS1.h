/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDS1_H
#define	HEFFDS1_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>
#include "EvolDF1nlep.h"
#include "EvolDB1Mll.h"
#include <sstream>

class HeffDS1{
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDS1(const StandardModel & SM);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDS1();
    
    /**
     * 
     * @brief the effective hamiltonian at the scale mu for \f$ K \rightarrow \pi \pi \f$
     * @brief with the current current open-charm operator contribution. 
     * @param mu is the low energy scale
     * @param scheme indicates the reonrmalization scheme
     * @return the effective hamiltonian at the scale mu for \f$ K \rightarrow \pi \pi \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffDS1PP(double mu, schemes scheme = NDR);
    
    gslpp::vector<gslpp::complex>** ComputeCoeffDS1pnunu();
    
    gslpp::vector<gslpp::complex>** ComputeCoeffDS1mumu();
    
    WilsonCoefficient getCoeffDS1PP() const {
        return coeffds1;
    }
    
    WilsonCoefficient getCoeffDS1cc() const {
        return coeffds1cc;
    }
    
    WilsonCoefficient getCoeffDS1pnunu() const {
        return coeffds1pnunu;
    }
    
    WilsonCoefficient getCoeffDS1mumu() const {
        return coeffds1mumu;
    }
    
    EvolDF1nlep getUDF1B() const {
        return u;
    }
    
    EvolDB1Mll getUDF1M() const {
        return uM;
    }

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;
    
    WilsonCoefficient coeffds1, coeffds1cc, coeffds1pnunu, coeffds1mumu;
    EvolDF1nlep u;
    EvolDB1Mll uM;
    
    gslpp::vector<gslpp::complex> DS1cce, DS1cc;
    
    /**
     * 
     * @brief compute the matching at the charm threshold within the SM, NDR scheme implemented
     * @return the SM effective hamiltonian below the charm threshold for \f$ K \rightarrow \pi \pi \f$
     */
    void CharmMatch();
};


#endif	/* HEFFDS1_H */

