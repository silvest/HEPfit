/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDF1_H
#define	HEFFDF1_H

#include "StandardModel.h"
#include "StandardModelMatching.h"
#include "WilsonCoefficient.h"
#include "EvolDF1.h"
#include <map>

#define N_OPS 15 // number of operators in the basis

class HeffDF1 {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDF1(unsigned int nops, std::string blocks, const StandardModel & SM, orders order = NLO, orders_qed order_qed = NO_QED);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDF1();
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> K^*ll decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    gslpp::vector<gslpp::complex>** ComputeCoeff_s(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> K^*ll decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    //gslpp::vector<gslpp::complex>** ComputeCoeffprime(double mu, schemes scheme = NDR);
    
    
    EvolDF1 getEvol() const {
        return evolDF1;
    }
    
    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;
    
    WilsonCoefficient coeff;
    EvolDF1 evolDF1;
    
    std::string blocks;
    unsigned int nops;
    double mu_cache;
    schemes scheme_cache;
    std::vector<double> Vmu_cache;
    std::vector<WilsonCoefficient> WC_cache;

};

#endif	/* HEFFDF1_H */
