/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef HEFFDF1_DIUJLKNU_H
#define HEFFDF1_DIUJLKNU_H


class StandardModel;
class EvolDF1_diujlknu;

#include "QCD.h"
#include "WilsonCoefficient.h"
#include "gslpp.h"
#include <sstream>
#include <memory>



/**
 * @class HeffDF1_diujlknu
 * @brief \f$ |\Delta F = 1 | \f$ Hamiltonian Class for \f$ d_i \rightarrow u_j \ell_k \nu \f$ transitions (leptonic and semileptonic charged-current meson decays)
 * @details The \f$ |\Delta F = 1 | \f$ Hamiltonian Class for \f$ d_i \rightarrow u_j \ell_k \nu \f$ transitions in the JMS basis ordered as CnueduVLLkkij, CnueduVLRkkij, CnueduSRRkkij, CnueduSRLkkij, CnueduTRRkkij
 * 
 */
class HeffDF1_diujlknu {
public:
/**
 * @brief HeffDF1_Plepnu constructor 
 * @param SM an object of StandardModel class
 * @param SM_Matching an object of StandardModelMatching class
 */
    HeffDF1_diujlknu(const StandardModel & SM);
    /**
     * @brief HeffDC1 destructor
     */
    virtual ~HeffDF1_diujlknu();

    /**
     * 
     * @param i the flavour of the down-type quark
     * @param j the flavour of the up-type quark
     * @param k the flavour of the charged lepton
     * @param mu the scale at which the coefficients should be evaluated
     * @return short distance contribution to \f$ \bar{d_i} u_j \bar{\ell_k} \nu \f$ transitions in the JMS basis ordered as CnueduVLLkkij, CnueduVLRkkij, CnueduSRRkkij, CnueduSRLkkij, CnueduTRRkkij
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffdiujleptonknu(int i, int j, int k, double mu);

    
    /**
     * 
     * @return an object of StandardModel class
     */
    const StandardModel& GetModel() const {
        return model;
    }
    
    
private:
    
    const StandardModel& model;
    WilsonCoefficient coeffdiujleptonknu;
    std::unique_ptr<EvolDF1_diujlknu> evolDF1;
    
};

#endif /* HEFFDF1_DIUJLKNU_H */

