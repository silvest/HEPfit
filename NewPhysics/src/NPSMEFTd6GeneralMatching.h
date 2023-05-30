/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSMEFTD6GENERALMATCHING_H
#define NPSMEFTD6GENERALMATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"

class NPSMEFTd6General;

/**
 * @class NPSMEFTd6GeneralMatching
 * @ingroup NewPhysics
 * @brief A class for the matching in the NPSMEFTd6_General model at the scale @f$ \mu_W @f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  This class, after update, contains the SMEFT coefficients at the scale @f$ \mu_W @f$ defined in the SMEFT model
 */
class NPSMEFTd6GeneralMatching : public StandardModelMatching {
public:
    NPSMEFTd6GeneralMatching(const NPSMEFTd6General & NPSMEFTd6General_i);
    
    virtual ~NPSMEFTd6GeneralMatching();
    
    /**
     *
     * @brief Updates to new FlavourWilsonCoefficient parameter sets.
     * @return
     */
    
    void updateNPSMEFTd6GeneralParameters();
    
private:
    const NPSMEFTd6General & mySMEFT;
    
};

#endif /* NPSMEFTD6GENERALMATCHING_H */


