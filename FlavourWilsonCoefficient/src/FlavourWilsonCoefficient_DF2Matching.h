/* 
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOURWILSONCOEFFICIENT_DF2MATCHING_H
#define FLAVOURWILSONCOEFFICIENT_DF2MATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"

class FlavourWilsonCoefficient_DF2;

/**
 * @class FlavourWilsonCoefficient_DF2Matching
 * @ingroup FlavourWilsonCoefficient
 * @brief A class for the matching in the FlavourWilsonCoefficient_DF2Matching. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class FlavourWilsonCoefficient_DF2Matching : public StandardModelMatching {
public:
    /**
     * @brief FlavourWilsonCoefficient_DF2Matching constructor
     * @param[in] An object of the FlavourWilsonCoefficient_DF2 class
     */
    FlavourWilsonCoefficient_DF2Matching(const FlavourWilsonCoefficient_DF2& FWC_i);
    
        /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{d} \f$ 
     * @return return the vector of SM Wilson coefficients (SM + NP)
     */
    virtual  std::vector<WilsonCoefficient>& CMdbd2() ;
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{s} \f$ 
     * @return return the vector of SM Wilson coefficients (SM + NP)
     */
    virtual   std::vector<WilsonCoefficient>& CMdbs2() ;

    /**
     * 
     * @brief \f$ \Delta C = 2 \f$,
     * @return return the vector of SM Wilson coefficients (SM + NP)
     */
    virtual   std::vector<WilsonCoefficient>& CMdd2() ;
    
    /**
     * 
     * @brief \f$ \Delta S = 2 \f$ 
     * @return return the vector of SM Wilson coefficients (SM + NP)
     */
    virtual   std::vector<WilsonCoefficient>& CMdk2() ;
       
protected:
    std::vector<WilsonCoefficient> vmcdbd2, vmcdbs2, vmcdc2, vmcds2;


private:
    const FlavourWilsonCoefficient_DF2& FWC;///< An object of the %FLAVOURWILSONCOEFFICIENT_DF2 class.
    WilsonCoefficient mcdbd2, mcdbs2, mcdc2, mcds2;

};

#endif /* FLAVOURWILSONCOEFFICIENT_DF2MATCHING_H */

