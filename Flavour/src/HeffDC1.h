/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDC1_H
#define	HEFFDC1_H

#include <StandardModel.h>
/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>
#include "EvolDC1.h"
#include "EvolDC1Buras.h"
#include <sstream>

/**
 * @class HeffDC1
 * @brief \f$ |\Delta C = 1 | \f$ Hamiltonian Class
 * @details The Hamiltonian for \f$ |\Delta C = 1 | \f$ processes as \f$ D^{0} \, \rightarrow \, \pi^{+} \, \pi^{-} \f$
 * and \f$ D^{0} \, \rightarrow \, K^{+} \, K^{-} \f$ 
 */
class HeffDC1 {
public:
/**
 * @brief HeffDC1 constructor 
 * @param SM an object of StandardModel class
 * @param SM_Matching an object of StandardModelMatching class
 */
    HeffDC1(const StandardModel & SM, StandardModelMatching & SM_Matching);
    /**
     * @brief HeffDC1 destructor
     */
    virtual ~HeffDC1();
    /**
     * @brief a method returning the evolved Wilson related to \f$ D^{0} \, \rightarrow \, \pi^{+} \, \pi^{-} \f$ 
     * @details it returns the Wilson coefficients of  the process \f$ D^{0} \, \rightarrow \, \pi^{+} \, \pi^{-} \f$
     * evolved to the low scale \f$ \mu \f$ with the associated coupling constants and CKM factors
     * @param mu a double for the low scale of the evolution 
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @return a vector<complex> pointer to pointer method for the evolved Wilson coefficients 
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffDC1_pi(double mu, schemes scheme = NDR);
    /**
     * @brief a method returning the evolved Wilson related to \f$ D^{0} \, \rightarrow \, K^{+} \, K^{-} \f$ 
     * @details it returns the Wilson coefficients of the process \f$ D^{0} \, \rightarrow \, K^{+} \, K^{-} \f$
     * evolved down to the low scale \f$ \mu \f$ with the associated coupling constants and CKM factors
     * @param mu a double for the low scale of the evolution 
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @return a vector<complex> pointer to pointer method for the evolved Wilson coefficients 
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffDC1_K(double mu, schemes scheme = NDR);
    /**
     *  
     * @return an object of WilsonCoefficient class
     */
    WilsonCoefficient getCoeffDC1() const {
        return coeffdc1;
    }
    /**
     *  
     * @return an object of EvolDC1 class
     */
    EvolDC1 getUDC1() const {
        return ug;
    }
    /**
     * 
     * @return an object of EvolDC1Buras class
     */
    EvolDC1Buras getUDC1Buras() const {
        return u;
    }
    /**
     * 
     * @return an object of StandardModel class
     */
    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;
    WilsonCoefficient coeffdc1, coeffdc1g;
    EvolDC1 ug;
    EvolDC1Buras u;
    gslpp::matrix<gslpp::complex> ckm, COEFF_pi, COEFF_K;
};

#endif	/* HEFFDC1_H */

