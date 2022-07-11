/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef REALWETHEFFDF1_H
#define REALWETHEFFDF1_H

#include "StandardModel.h"
#include "gslpp.h"

/**
 * @class RealWETHeffDF1
 * @ingroup FlavourWilsonCoefficent
 * @brief %Model for NP contributions to HeffDF1
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 * 
 * 
 * @anchor RealWETHeffDF1Parameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %RealWETHeffDF1 model are summarised below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CSL, %CSR, %CVL, %CVR</td>
 *   <td class="mod_symb">@f$C_S^L@f$, @f$C_S^R@f$, @f$C_V^L@f$
 *   @f$C_V^R@f$</td>
 *   <td class="mod_desc">The left and right handed scalar and
 *   vector couplings for the RealWeakEFTCCPM %Model</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CS, %CP, %CV, %CA</td>
 *   <td class="mod_symb">@f$C_S@f$, @f$C_P@f$, @f$C_V@f$
 *   @f$C_A@f$</td>
 *   <td class="mod_desc">The scalar, pseudo-scalar, vector,
 *   and axialvector couplingsfor the %RealWeakEFTCC %Model</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CT</td>
 *   <td class="mod_symb">@f$C_T@f$</td>
 *   <td class="mod_desc">The tensor coupling for both
 *   the %RealWeakEFTCC and RealWeakEFTCCPM %Model</td>
 * </tr>
 * </table>
 * 
 */
class RealWETHeffDF1: public StandardModel {
public:

    static const int NRealWETHeffDF1vars = 8;

    static const std::string RealWETHeffDF1vars[NRealWETHeffDF1vars];
    
    /**
     * @brief RealWETHeffDF1 constructor
     */
    RealWETHeffDF1();
    
    /**
     * @brief RealWETHeffDF1 destructor
     */
    ~RealWETHeffDF1();
    
    /**
     * @brief The post-update method for %RealWETHeffDF1.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief Initializes the %RealWETHeffDF1 parameters found in the argument.
     * @param[in] DPars a map containing the parameters (all as double) to be used in Monte Carlo
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The pre-update method for %RealWETHeffDF1
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();
    
    /**
     * @brief The update method for %RealWETHeffDF1.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The post-update method for %RealWETHeffDF1.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %RealWETHeffDF1
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to set a flag of %RealWETHeffDF1.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     *
     * @return \f$ C_9^{11}\f$
     */
    double getC9_11() const {
        return C9_11;
    }
    
    /**
     *
     * @return \f$ C_9^{22}\f$
     */
    double getC9_22() const {
        return C9_22;
    }
    
    /**
     *
     * @return \f$ C_{10}^{11}\f$
     */
    double getC10_11() const {
        return C10_11;
    }
    
    /**
     *
     * @return \f$ C_{10}^{22}\f$
     */
    double getC10_22() const {
        return C10_22;
    }
    
    /**
     *
     * @return \f$ C_9^{\prime 11}\f$
     */
    double getC9p_11() const {
        return C9p_11;
    }
    
    /**
     *
     * @return \f$ C_9^{\prime 22}\f$
     */
    double getC9p_22() const {
        return C9p_22;
    }
    
    /**
     *
     * @return \f$ C_{10}^{\prime 11}\f$
     */
    double getC10p_11() const {
        return C10p_11;
    }
    
    /**
     *
     * @return \f$ C_{10}^{\prime 22}\f$
     */
    double getC10p_22() const {
        return C10p_22;
    }
    
protected: 
    
    /**
     * @brief A method to set the value of a parameter of %FlavourWilsonCoefficient.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string, const double&);

private:

    double C9_11, C10_11, C9p_11, C10p_11;///< The Wilson Coefficients for electric currents
    double C9_22, C10_22, C9p_22, C10p_22;///< The Wilson Coefficients for muonic currents
};

#endif /* REALWETHEFFDF1_H */
