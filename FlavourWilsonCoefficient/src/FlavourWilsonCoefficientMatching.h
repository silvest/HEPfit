/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOURWILSONCOEFFICIENTMATCHING_H
#define	FLAVOURWILSONCOEFFICIENTMATCHING_H

#include <gslpp.h>
#include "StandardModelMatching.h"

class FlavourWilsonCoefficient;

/**
 * @class FlavourWilsonCoefficientMatching
 * @ingroup FlavourWilsonCoefficient
 * @brief A class for the matching in the FlavourWilsonCoefficient. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class FlavourWilsonCoefficientMatching : public StandardModelMatching {
public:
    FlavourWilsonCoefficientMatching(const FlavourWilsonCoefficient & FlavourWilsonCoefficient_i);
    
    /**
     *
     * @brief Updates to new FlavourWilsonCoefficient parameter sets.
     * @return
     */
    
    void updateFlavourWilsonCoefficientParameters();

private:
    const FlavourWilsonCoefficient & myFlavourWilsonCoefficient;
    gslpp::matrix<complex> myCKM;
};

#endif	/* FLAVOURWILSONCOEFFICIENTMATCHING_H */

