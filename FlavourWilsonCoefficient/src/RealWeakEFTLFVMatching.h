/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef REALWEAKEFTLFVMATCHING_H
#define REALWEAKEFTLFVMATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"

class RealWeakEFTLFV;

/**
 * @class NPSMEFTd6Matching
 * @ingroup NewPhysics
 * @brief A class for the matching in the NPSMEFTd6.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class RealWeakEFTLFVMatching : public StandardModelMatching {
public:
    RealWeakEFTLFVMatching(const RealWeakEFTLFV & RealWeakEFTLFV_i);
    
    virtual ~RealWeakEFTLFVMatching();
    
    /**
     *
     * @brief Updates to new FlavourWilsonCoefficient parameter sets.
     * @return
     */
    
    void updateRealWeakEFTLFVParameters();
    
    std::vector<WilsonCoefficient>& CMbsg();
    
    std::vector<WilsonCoefficient>& CMbsmm();
    
    std::vector<WilsonCoefficient>& CMbdmm();
    
    std::vector<WilsonCoefficient>& CMprimebsg();
    
    std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton);
    
    std::vector<WilsonCoefficient>& CMprimeBMll(QCD::lepton lepton);

private:
    const RealWeakEFTLFV & myRealWeakEFTLFV;
    
    gslpp::complex SMEFT_factor, loop_factor;
    
    gslpp::complex C7NP;
    gslpp::complex C7pNP;
    gslpp::complex C8NP;
    gslpp::complex C8pNP;
    
    gslpp::complex C9NPe;
    gslpp::complex C9pNPe;
    gslpp::complex C10NPe;
    gslpp::complex C10pNPe;
    gslpp::complex CSNPe;
    gslpp::complex CSpNPe;
    gslpp::complex CPNPe;
    gslpp::complex CPpNPe;
    
    gslpp::complex C9NPmu;
    gslpp::complex C9pNPmu;
    gslpp::complex C10NPmu;
    gslpp::complex C10pNPmu;
    gslpp::complex CSNPmu;
    gslpp::complex CSpNPmu;
    gslpp::complex CPNPmu;
    gslpp::complex CPpNPmu;
    
    double WCscale;
    
    WilsonCoefficient mcbsg, mcprimebsg;
    WilsonCoefficient mcbsmm;
    WilsonCoefficient mcbdmm;
    WilsonCoefficient mcBMll, mcprimeBMll;
    std::vector<WilsonCoefficient> vmcbsg, vmcprimebsg;
    std::vector<WilsonCoefficient> vmcbsmm;
    std::vector<WilsonCoefficient> vmcbdmm;
    std::vector<WilsonCoefficient> vmcBMll, vmcprimeBMll;
};

#endif /* REALWEAKEFTLFVMATCHING_H */

