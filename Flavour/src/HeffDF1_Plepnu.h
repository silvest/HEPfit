/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef HEFFDF1_PLEPNU_H
#define HEFFDF1_PLEPNU_H


class StandardModel;

#include "QCD.h"
#include "WilsonCoefficient.h"
#include "gslpp.h"
#include <sstream>
#include <memory>



/**
 * @class HeffDF1_Plepnu
 * @brief \f$ |\Delta F = 1 | \f$ Hamiltonian Class for Pion leptonic decay
 * @details The Hamiltonian for \f$ |\Delta F = 1 | \f$ processes as \f$ \pi \, \rightarrow \, \mu \nu \f$
 * 
 */
class HeffDF1_Plepnu {
public:
/**
 * @brief HeffDF1_Plepnu constructor 
 * @param SM an object of StandardModel class
 * @param SM_Matching an object of StandardModelMatching class
 */
    HeffDF1_Plepnu(const StandardModel & SM);
    /**
     * @brief HeffDC1 destructor
     */
    virtual ~HeffDF1_Plepnu();

    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ D \rightarrow \lepton \nu \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffuleptonnu(QCD::meson meson_i, QCD::lepton lepton_i);

    
    /**
     * 
     * @return an object of StandardModel class
     */
    const StandardModel& GetModel() const {
        return model;
    }
    
    
private:
    
    const StandardModel& model;
    WilsonCoefficient coeffuleptonnu;
    
};

#endif /* HEFFDF1_PLEPNU_H */

