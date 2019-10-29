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
 * @class RealWeakEFTLFVMatching
 * @ingroup FlavourWilsonCoefficient
 * @brief A class for the matching in the RealWeakEFTLFVMatching. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class RealWeakEFTLFVMatching : public StandardModelMatching {
public:
    /**
     * @brief FlavourWilsonCoefficient_DF2Matching constructor
     * @param[in] RealWeakEFTLFV_i An object of the RealWeakEFTLFV class
     */
    RealWeakEFTLFVMatching(const RealWeakEFTLFV & RealWeakEFTLFV_i);
    
    /**
     * @brief FlavourWilsonCoefficient_DF2Matching destructor
     */
    virtual ~RealWeakEFTLFVMatching();
    
    /**
     *
     * @brief Updates to new RealWeakEFTLFVMatching parameter sets.
     */
    void updateRealWeakEFTLFVParameters();
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    std::vector<WilsonCoefficient>& CMbsg();
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    std::vector<WilsonCoefficient>& CMprimebsg();
    
    /**
     * @brief Wilson coefficient for \f$ B_{s} \rightarrow \mu \mu \f$
     * @return Wilson coefficient for \f$ B_{s} \rightarrow \mu \mu \f$
     */
    std::vector<WilsonCoefficient>& CMbsmm();
    
    /**
     * @brief Wilson coefficient for \f$ B_{d} \rightarrow \mu \mu \f$
     * @return Wilson coefficient for \f$ B_{d} \rightarrow \mu \mu \f$
     */
    std::vector<WilsonCoefficient>& CMbdmm();
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow K^* l^{+} l^{-} \f$
     */
    std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton);
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow K^* l^{+} l^{-} \f$
     */
    std::vector<WilsonCoefficient>& CMprimeBMll(QCD::lepton lepton);

private:
    const RealWeakEFTLFV & myRealWeakEFTLFV;///< An object of the %RealWeakEFTLFV class.
    
    gslpp::complex SMEFT_factor, loop_factor;
    
    gslpp::complex C7NP;///< the Wilson Coefficient \f$C_{7}^{NP}\f$.
    gslpp::complex C7pNP;///< the Wilson Coefficient \f$C_{7}^{\prime NP}\f$.
    gslpp::complex C8NP;///< the Wilson Coefficient \f$C_{8}^{NP}\f$.
    gslpp::complex C8pNP;///< the Wilson Coefficient \f$C_{8}^{\prime NP}\f$.
    
    gslpp::complex C9NPe;///< the Wilson Coefficient \f$C_{9,e}^{NP}\f$.
    gslpp::complex C9pNPe;///< the Wilson Coefficient \f$C_{9,e}^{\prime NP}\f$.
    gslpp::complex C10NPe;///< the Wilson Coefficient \f$C_{10,e}^{NP}\f$.
    gslpp::complex C10pNPe;///< the Wilson Coefficient \f$C_{10,e}^{\prime NP}\f$.
    gslpp::complex CSNPe;///< the Wilson Coefficient \f$C_{S,e}^{NP}\f$.
    gslpp::complex CSpNPe;///< the Wilson Coefficient \f$C_{S,e}^{\prime NP}\f$.
    gslpp::complex CPNPe;///< the Wilson Coefficient \f$C_{P,e}^{NP}\f$.
    gslpp::complex CPpNPe;///< the Wilson Coefficient \f$C_{P,e}^{\prime NP}\f$.
    
    gslpp::complex C9NPmu;///< the Wilson Coefficient \f$C_{9,\mu}^{NP}\f$.
    gslpp::complex C9pNPmu;///< the Wilson Coefficient \f$C_{9,\mu}^{\prime NP}\f$.
    gslpp::complex C10NPmu;///< the Wilson Coefficient \f$C_{10,\mu}^{NP}\f$.
    gslpp::complex C10pNPmu;///< the Wilson Coefficient \f$C_{10,\mu}^{\prime NP}\f$.
    gslpp::complex CSNPmu;///< the Wilson Coefficient \f$C_{S,\mu}^{NP}\f$.
    gslpp::complex CSpNPmu;///< the Wilson Coefficient \f$C_{S,\mu}^{\prime NP}\f$.
    gslpp::complex CPNPmu;///< the Wilson Coefficient \f$C_{P,\mu}^{NP}\f$.
    gslpp::complex CPpNPmu;///< the Wilson Coefficient \f$C_{P,\mu}^{\prime NP}\f$.
    
    double WCscale;///< The scale of the Wilson coefficients.
    
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

