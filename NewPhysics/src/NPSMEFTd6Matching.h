/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSMEFTD6MATCHING_H
#define NPSMEFTD6MATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"

class NPSMEFTd6;

/**
 * @class NPSMEFTd6Matching
 * @ingroup NewPhysics
 * @brief A class for the matching in the NPSMEFTd6.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPSMEFTd6Matching : public StandardModelMatching {
public:
    NPSMEFTd6Matching(const NPSMEFTd6 & NPSMEFTd6_i);
    
    virtual ~NPSMEFTd6Matching();
    
    /**
     *
     * @brief Updates to new FlavourWilsonCoefficient parameter sets.
     * @return
     */
    
    void updateNPSMEFTd6Parameters();
    
    std::vector<WilsonCoefficient>& CMbsg();
    
    std::vector<WilsonCoefficient>& CMbsmm();
    
    std::vector<WilsonCoefficient>& CMbdmm();
    
    std::vector<WilsonCoefficient>& CMprimebsg();
    
    std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton);
    
    std::vector<WilsonCoefficient>& CMprimeBMll(QCD::lepton lepton);

private:
    const NPSMEFTd6 & myNPSMEFTd6;
    gslpp::matrix<gslpp::complex> myCKM;
    
    gslpp::complex C9;
    gslpp::complex C10;
    gslpp::complex C7g;
    gslpp::complex C8g;
    
    gslpp::complex C9p;
    gslpp::complex C10p;
    gslpp::complex C7gp;
    gslpp::complex C8gp;
    
    double Muw;
    
    WilsonCoefficient mcbsg, mcprimebsg;
    WilsonCoefficient mcbsmm;
    WilsonCoefficient mcbdmm;
    WilsonCoefficient mcBMll, mcprimeBMll;
    std::vector<WilsonCoefficient> vmcbsg, vmcprimebsg;
    std::vector<WilsonCoefficient> vmcbsmm;
    std::vector<WilsonCoefficient> vmcbdmm;
    std::vector<WilsonCoefficient> vmcBMll, vmcprimeBMll;
};

#endif /* NPSMEFTD6MATCHING_H */

