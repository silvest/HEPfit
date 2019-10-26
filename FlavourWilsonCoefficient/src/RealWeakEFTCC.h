/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef REALWEAKEFTCC_H
#define REALWEAKEFTCC_H

#include "StandardModel.h"
#include "gslpp.h"

/**
 * @class RealWeakEFTCC
 * @ingroup FlavourWilsonCoefficent
 * @brief %Model for NP contributions to charged current
 * processes like @f$b\to c@f$ decays.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class RealWeakEFTCC: public StandardModel {
public:

    static const int NRealWeakEFTCCvars = 5;

    static const std::string RealWeakEFTCCvars[NRealWeakEFTCCvars];
    static const std::string RealWeakEFTCCPMvars[NRealWeakEFTCCvars];
    
    /**
     * @brief RealWeakEFTCC constructor
     */
    RealWeakEFTCC(const unsigned int basis);
    
    /**
     * @brief RealWeakEFTCC destructor
     */
    ~RealWeakEFTCC();
    
    /**
     * @brief The post-update method for %RealWeakEFTCC.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief Initializes the %RealWeakEFTCC parameters found in the argument.
     * @param[in] DPars a map containing the parameters (all as double) to be used in Monte Carlo
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The pre-update method for %RealWeakEFTCC
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();
    
    /**
     * @brief The update method for %RealWeakEFTCC.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The post-update method for %RealWeakEFTCC.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %RealWeakEFTCC
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to set a flag of %RealWeakEFTCC.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     *
     * @return \f$ C_S or C_{S}^{L}\f$
     */
    double getCCC1() const {
        return C1;
    }
    
    /**
     *
     * @return \f$ C_P or C_{S}^{R}\f$
     */
    double getCCC2() const {
        return C2;
    }
    
    /**
     *
     * @return \f$ C_V or C_{V}^{L}\f$
     */
    double getCCC3() const {
        return C3;
    }
    
        /**
     *
     * @return \f$ C_A or or C_{V}^{R}$\f
     */
    double getCCC4() const {
        return C4;
    }
    
    /**
     *
     * @return \f$ C_T\f$
     */
    double getC5() const {
        return C5;
    }
    
protected: 
    
    /**
     * @brief A method to set the value of a parameter of %FlavourWilsonCoefficient.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string, const double&);

private:

    double C1, C2, C3, C4, C5;///< The Wilson Coefficients.
    const unsigned int basis;///< The choice of basis.
};

#endif /* REALWEAKEFTCC_H */

