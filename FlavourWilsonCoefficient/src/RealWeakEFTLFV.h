/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef REALWEAKEFTLFV_H
#define REALWEAKEFTLFV_H

#include "StandardModel.h"
#include "gslpp.h"
#include "RealWeakEFTLFVMatching.h"

/**
 * @class RealWeakEFTLFV
 * @ingroup FlavourWilsonCoefficent
 * @brief %Model for WEFT LFV contributions to @f$\Delta F=1@f$
 * processes like @f$ b\to s@f$ decays.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * 
 * @anchor RealWeakEFTLFVParameters
 * <h3>%Model parameters</h3>
 *
 * 
 * 
 * The model parameters of %FlavourWilsonCoefficient model are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C7, %C7p</td>
 *   <td class="mod_symb">@f$C_7@f$ and @f$C_7^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the dipole operator
 *   @f$O_7@f$ and the chirality flipped @f$O_7^\prime@f$</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C8, %C8p</td>
 *   <td class="mod_symb">@f$C_8@f$ and @f$C_8^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the chromomagnetic operator
 *   @f$O_8@f$ and the chirality flipped @f$O_8^\prime@f$</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C9_11, %C9p_11</td>
 *   <td class="mod_symb">@f$C_{9,e}@f$ and @f$C_{9,e}^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the semileptonic operator
 *   @f$O_{9,e}@f$ and the chirality flipped @f$O_{9,e}^\prime@f$ coupling to electrons</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C10_11, %C10p_11</td>
 *   <td class="mod_symb">@f$C_{10,e}@f$ and @f$C_{10,e}^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the semileptonic operator
 *   @f$O_{10,e}@f$ and the chirality flipped @f$O_{10,e}^\prime@f$ coupling to electrons</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CS_11, %CSp_11</td>
 *   <td class="mod_symb">@f$C_{S,e}@f$ and @f$C_{S,e}^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the scalar operator
 *   @f$O_{S,e}@f$ and the chirality flipped @f$O_{S,e}^\prime@f$ coupling to electrons</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CP_11, %CPp_11</td>
 *   <td class="mod_symb">@f$C_{P,e}@f$ and @f$C_{P,e}^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the pseudo-scalar operator
 *   @f$O_{P,e}@f$ and the chirality flipped @f$O_{P,e}^\prime@f$ coupling to electrons</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C9_22, %C9p_22</td>
 *   <td class="mod_symb">@f$C_{9,\mu}@f$ and @f$C_{9,\mu}^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the semileptonic operator
 *   @f$O_{9,\mu}@f$ and the chirality flipped @f$O_{9,\mu}^\prime@f$ coupling to muons</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%C10_22, %C10p_22</td>
 *   <td class="mod_symb">@f$C_{10,\mu}@f$ and @f$C_{10,\mu}^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the semileptonic operator
 *   @f$O_{10,\mu}@f$ and the chirality flipped @f$O_{10,\mu}^\prime@f$ coupling to muons</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CS_22, %CSp_22</td>
 *   <td class="mod_symb">@f$C_{S,\mu}@f$ and @f$C_{S,\mu}^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the scalar operator
 *   @f$O_{S,\mu}@f$ and the chirality flipped @f$O_{S,\mu}^\prime@f$ coupling to muons</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%CP_22, %CPp_22</td>
 *   <td class="mod_symb">@f$C_{P,\mu}@f$ and @f$C_{P,\mu}^\prime@f$</td>
 *   <td class="mod_desc">The Wilson coefficient of the pseudo-scalar operator
 *   @f$O_{P,\mu}@f$ and the chirality flipped @f$O_{P,\mu}^\prime@f$ coupling to muons</td>
 * </tr>
 * </table>
 */
 
class RealWeakEFTLFV: public StandardModel {
public:

    static const int NRealWeakEFTLFVvars = 21;

    static const std::string RealWeakEFTLFVvars[NRealWeakEFTLFVvars];
    
    /**
     * @brief RealWeakEFTLFV constructor
     */
    RealWeakEFTLFV();
    
    /**
     * @brief RealWeakEFTLFV destructor
     */
    ~RealWeakEFTLFV();
    
    /**
     * @brief The post-update method for %RealWeakEFTLFV.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief Initializes the %RealWeakEFTLFV parameters found in the argument.
     * @param[in] DPars a map containing the parameters (all as double) to be used in Monte Carlo
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The pre-update method for %RealWeakEFTLFV
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();
    
    /**
     * @brief The update method for %RealWeakEFTLFV.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The post-update method for %RealWeakEFTLFV.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %RealWeakEFTLFV
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to set a flag of %RealWeakEFTLFV.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     * @brief A get method to access the member reference of type %RealWeakEFTLFVMatching.
     * @return a reference to a %RealWeakEFTLFVMatching object
     */
    virtual RealWeakEFTLFVMatching& getMatching() const
    {
        return ReWEFTM.getObj();
    }
    
    /**
     *
     * @return \f$ C_7\f$
     */
    double getC7() const {
        return C7;
    }
    
    /**
     *
     * @return \f$ C_7'\f$
     */
    double getC7p() const {
        return C7p;
    }
    
    /**
     *
     * @return \f$ C_8\f$
     */
    double getC8() const {
        return C8;
    }
    
        /**
     *
     * @return \f$ C_8'$\f
     */
    double getC8p() const {
        return C8p;
    }
    
    /**
     *
     * @return \f$ C_9^{11}\f$
     */
    double getC9_11() const {
        return C9_11;
    }
    
        /**
     *
     * @return \f$ C_9'^{11}\f$
     */
    double getC9p_11() const {
        return C9p_11;
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
     * @return \f$ C_{10}'^{11}\f$
     */
    double getC10p_11() const {
        return C10p_11;
    }
    
    /**
     *
     * @return \f$ C_S^{11}\f$
     */
    double getCS_11() const {
        return CS_11;
    }
    
        /**
     *
     * @return \f$ C_S'^{11}\f$
     */
    double getCSp_11() const {
        return CSp_11;
    }
    
    /**
     *
     * @return \f$ C_P^{11}\f$
     */
    double getCP_11() const {
        return CP_11;
    }
    
        /**
     *
     * @return \f$ C_P'^{11}\f$
     */
    double getCPp_11() const {
        return CPp_11;
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
     * @return \f$ C_9'^{22}\f$
     */
    double getC9p_22() const {
        return C9p_22;
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
     * @return \f$ C_{10}'^{22}\f$
     */
    double getC10p_22() const {
        return C10p_22;
    }
    
    /**
     *
     * @return \f$ C_S^{22}\f$
     */
    double getCS_22() const {
        return CS_22;
    }
    
        /**
     *
     * @return \f$ C_S'^{22}\f$
     */
    double getCSp_22() const {
        return CSp_22;
    }
    
    /**
     *
     * @return \f$ C_P^{22}\f$
     */
    double getCP_22() const {
        return CP_22;
    }
    
        /**
     *
     * @return \f$ C_P'^{22}\f$
     */
    double getCPp_22() const {
        return CPp_22;
    }
    
    /**
     *
     * @return the matching scale of the Wilson coefficients
     */
    double getWCscale() const {
        return WCscale;
    }
    
protected: 
    
    /**
     * @brief A method to set the value of a parameter of %RealWeakEFTLFV.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string, const double&);
    mutable Matching<RealWeakEFTLFVMatching,RealWeakEFTLFV> ReWEFTM;///< The RealWeakEFTLFVMatching object.

private:

    double C7, C7p, C8, C8p;///< The real Wilson Coefficients.
    double C9_11, C9p_11, C10_11, C10p_11, CS_11, CSp_11, CP_11, CPp_11;///< The real Wilson Coefficients for electrons.
    double C9_22, C9p_22, C10_22, C10p_22, CS_22, CSp_22, CP_22, CPp_22;///< The real Wilson Coefficients for muons.
    
    double WCscale;///< The scale of the  Wilson Coefficients.
    
      
};

#endif /* REALWEAKEFTLFV_H */

