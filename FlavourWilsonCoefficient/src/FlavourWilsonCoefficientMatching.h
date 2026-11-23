/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
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
    /**
     * @brief FlavourWilsonCoefficient constructor
     * @param[in] FlavourWilsonCoefficient_i An object of the FlavourWilsonCoefficient class
     */
    FlavourWilsonCoefficientMatching(const FlavourWilsonCoefficient & FlavourWilsonCoefficient_i);
    
    /**
     * @brief FlavourWilsonCoefficientMatching destructor
     */
    virtual ~FlavourWilsonCoefficientMatching();
    
    /**
     *
     * @brief Updates to new FlavourWilsonCoefficient parameter sets.
     * @return
     */
    void updateFlavourWilsonCoefficientParameters();
    
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
    const FlavourWilsonCoefficient & myFlavourWilsonCoefficient;///< An object of the %FlavourWilsonCoefficient class.
    gslpp::matrix<gslpp::complex> myCKM;///< The CKM matrix.
    gslpp::complex DC1;///< the modification Wilson Coefficient \f$C_{1}\f$.
    gslpp::complex DC2;///< the modification Wilson Coefficient \f$C_{2}\f$.
    gslpp::complex DC3;///< the modification Wilson Coefficient \f$C_{3}\f$.
    gslpp::complex DC4;///< the modification Wilson Coefficient \f$C_{4}\f$.
    gslpp::complex DC5;///< the modification Wilson Coefficient \f$C_{5}\f$.
    gslpp::complex DC6;///< the modification Wilson Coefficient \f$C_{6}\f$.
    gslpp::complex DC7;///< the modification Wilson Coefficient \f$C_{7}\f$.
    gslpp::complex DC8;///< the modification Wilson Coefficient \f$C_{8}\f$.
    gslpp::complex DC9;///< the modification Wilson Coefficient \f$C_{9}\f$.
    gslpp::complex DC10;///< the modification Wilson Coefficient \f$C_{10}\f$.
    gslpp::complex DC7g;///< the modification Wilson Coefficient \f$C_{7g}\f$.
    gslpp::complex DC8g;///< the modification Wilson Coefficient \f$C_{8g}\f$.
    
    gslpp::complex DC7p;///< the modification Wilson Coefficient \f$C_{7}^{\prime}\f$.
    gslpp::complex DC8p;///< the modification Wilson Coefficient \f$C_{8}^{\prime}\f$.
    gslpp::complex DC9p;///< the modification Wilson Coefficient \f$C_{9}^{\prime}\f$.
    gslpp::complex DC10p;///< the modification Wilson Coefficient \f$C_{10}^{\prime}\f$.
    gslpp::complex DC7gp;///< the modification Wilson Coefficient \f$C_{7g}^{\prime}\f$.
    gslpp::complex DC8gp;///< the modification Wilson Coefficient \f$C_{8g}^{\prime}\f$.
    
    WilsonCoefficient mcbsg, mcprimebsg;
    WilsonCoefficient mcbsmm;
    WilsonCoefficient mcbdmm;
    WilsonCoefficient mcBMll, mcprimeBMll;
    std::vector<WilsonCoefficient> vmcbsg, vmcprimebsg;
    std::vector<WilsonCoefficient> vmcbsmm;
    std::vector<WilsonCoefficient> vmcbdmm;
    std::vector<WilsonCoefficient> vmcBMll, vmcprimeBMll;
};

#endif	/* FLAVOURWILSONCOEFFICIENTMATCHING_H */

