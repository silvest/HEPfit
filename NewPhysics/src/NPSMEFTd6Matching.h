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
    
    gslpp::complex SMEFT_factor, loop_factor;
    
    gslpp::complex C7NP;
    gslpp::complex C8NP;
    gslpp::complex C7pNP;
    gslpp::complex C8pNP;
    
    gslpp::complex C9NPmu;
    gslpp::complex C10NPmu;
    gslpp::complex C9pNPmu;
    gslpp::complex C10pNPmu;
    gslpp::complex CSNPmu;
    gslpp::complex CSpNPmu;
    gslpp::complex CPNPmu;
    gslpp::complex CPpNPmu;
    
    gslpp::complex C9NPe;
    gslpp::complex C10NPe;
    gslpp::complex C9pNPe;
    gslpp::complex C10pNPe;
    gslpp::complex CSNPe;
    gslpp::complex CSpNPe;
    gslpp::complex CPNPe;
    gslpp::complex CPpNPe;
    
    gslpp::complex CLQ1_1123_tot, CLQ1_2223_tot;
    gslpp::complex CLQ3_1123_tot, CLQ3_2223_tot; 
    gslpp::complex CQe_2311_tot, CQe_2322_tot;
    
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

