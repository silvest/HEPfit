/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOURWILSONCOEFFICIENTMATCHING_H
#define	FLAVOURWILSONCOEFFICIENTMATCHING_H

#include "gslpp.h"
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
    
    std::vector<WilsonCoefficient>& CMbsg();
    
    std::vector<WilsonCoefficient>& CMprimebsg();
    
    std::vector<WilsonCoefficient>& CMBMll();
    
    std::vector<WilsonCoefficient>& CMprimeBMll();

private:
    const FlavourWilsonCoefficient & myFlavourWilsonCoefficient;
    gslpp::matrix<gslpp::complex> myCKM;
    gslpp::complex DC1;
    gslpp::complex DC2;
    gslpp::complex DC3;
    gslpp::complex DC4;
    gslpp::complex DC5;
    gslpp::complex DC6;
    gslpp::complex DC7;
    gslpp::complex DC8;
    gslpp::complex DC9;
    gslpp::complex DC10;
    gslpp::complex DC7g;
    gslpp::complex DC8g;
    
    gslpp::complex DC7p;
    gslpp::complex DC8p;
    gslpp::complex DC9p;
    gslpp::complex DC10p;
    gslpp::complex DC7gp;
    gslpp::complex DC8gp;
    
    WilsonCoefficient mcbsg, mcprimebsg;
    WilsonCoefficient mcBMll, mcprimeBMll;
    std::vector<WilsonCoefficient> vmcbsg, vmcprimebsg;
    std::vector<WilsonCoefficient> vmcBMll, vmcprimeBMll;
};

#endif	/* FLAVOURWILSONCOEFFICIENTMATCHING_H */

