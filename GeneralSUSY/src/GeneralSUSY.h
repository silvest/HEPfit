/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALSUSY_H
#define	GENERALSUSY_H

#include "SUSY.h"

/**
 * @addtogroup GeneralSUSY
 * @brief A module for the general MSSM.
 */

/**
 * @class GeneralSUSY
 * @brief A class for the general MSSM. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Sets all soft masses and trilinear couplings for squarks and sleptons.
 * 
 *
 * @anchor GeneralSUSYParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %GeneralSUSY model are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msQhat2_11r, %msQhat2_12r, %msQhat2_13r
 *                    <br>%msQhat2_21r, %msQhat2_22r, %msQhat2_23r</td>
 *                    <br>%msQhat2_31r, %msQhat2_32r, %msQhat2_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(m_Q^2)@f$</td>
 *   <td class="mod_desc">Real part of the squared soft mass matrix
 *   of the @f$SU(2)_L@f$ doublet squarks</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msQhat2_11i, %msQhat2_12i, %msQhat2_13i
 *                    <br>%msQhat2_21i, %msQhat2_22i, %msQhat2_23i</td>
 *                    <br>%msQhat2_31i, %msQhat2_32i, %msQhat2_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(m_Q^2)@f$</td>
 *   <td class="mod_desc">Imaginary part of the squared soft mass matrix
 *   of the @f$SU(2)_L@f$ doublet squarks</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msUhat2_11r, %msUhat2_12r, %msUhat2_13r
 *                    <br>%msUhat2_21r, %msUhat2_22r, %msUhat2_23r</td>
 *                    <br>%msUhat2_31r, %msUhat2_32r, %msUhat2_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(m_U^2)@f$</td>
 *   <td class="mod_desc">Real part of the squared soft mass matrix
 *   of the up-type  @f$SU(2)_L@f$ singlet squarks</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msUhat2_11i, %msUhat2_12i, %msUhat2_13i
 *                    <br>%msUhat2_21i, %msUhat2_22i, %msUhat2_23i</td>
 *                    <br>%msUhat2_31i, %msUhat2_32i, %msUhat2_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(m_U^2)@f$</td>
 *   <td class="mod_desc">Imaginary part of the squared soft mass matrix
 *   of the up-type  @f$SU(2)_L@f$ singlet squarks</td>
 * <tr>
 *   <td class="mod_name">%msDhat2_11r, %msDhat2_12r, %msDhat2_13r
 *                    <br>%msDhat2_21r, %msDhat2_22r, %msDhat2_23r</td>
 *                    <br>%msDhat2_31r, %msDhat2_32r, %msDhat2_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(m_D^2)@f$</td>
 *   <td class="mod_desc">Real part of the squared soft mass matrix
 *   of the down-type  @f$SU(2)_L@f$ singlet squarks</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msDhat2_11i, %msDhat2_12i, %msDhat2_13i
 *                    <br>%msDhat2_21i, %msDhat2_22i, %msDhat2_23i</td>
 *                    <br>%msDhat2_31i, %msDhat2_32i, %msDhat2_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(m_D^2)@f$</td>
 *   <td class="mod_desc">Imaginary part of the squared soft mass matrix
 *   of the down-type  @f$SU(2)_L@f$ singlet squarks</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msLhat2_11r, %msLhat2_12r, %msLhat2_13r
 *                    <br>%msLhat2_21r, %msLhat2_22r, %msLhat2_23r</td>
 *                    <br>%msLhat2_31r, %msLhat2_32r, %msLhat2_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(m_L^2)@f$</td>
 *   <td class="mod_desc">Real part of the squared soft mass matrix
 *   of the @f$SU(2)_L@f$ doublet sleptons</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msLhat2_11i, %msLhat2_12i, %msLhat2_13i
 *                    <br>%msLhat2_21i, %msLhat2_22i, %msLhat2_23i</td>
 *                    <br>%msLhat2_31i, %msLhat2_32i, %msLhat2_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(m_L^2)@f$</td>
 *   <td class="mod_desc">Imaginary part of the squared soft mass matrix
 *   of the @f$SU(2)_L@f$ doublet sleptons</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msEhat2_11r, %msEhat2_12r, %msEhat2_13r
 *                    <br>%msEhat2_21r, %msEhat2_22r, %msEhat2_23r</td>
 *                    <br>%msEhat2_31r, %msEhat2_32r, %msEhat2_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(m_E^2)@f$</td>
 *   <td class="mod_desc">Real part of the squared soft mass matrix
 *   of the up-type  @f$SU(2)_L@f$ singlet slepton</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msEhat2_11i, %msEhat2_12i, %msEhat2_13i
 *                    <br>%msEhat2_21i, %msEhat2_22i, %msEhat2_23i</td>
 *                    <br>%msEhat2_31i, %msEhat2_32i, %msEhat2_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(m_E^2)@f$</td>
 *   <td class="mod_desc">Imaginary part of the squared soft mass matrix
 *   of the up-type  @f$SU(2)_L@f$ singlet slepton</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msNhat2_11r, %msNhat2_12r, %msNhat2_13r
 *                    <br>%msNhat2_21r, %msNhat2_22r, %msNhat2_23r</td>
 *                    <br>%msNhat2_31r, %msNhat2_32r, %msNhat2_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(m_N^2)@f$</td>
 *   <td class="mod_desc">Real part of the squared soft mass matrix
 *   of the down-type  @f$SU(2)_L@f$ singlet slepton</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%msNhat2_11i, %msNhat2_12i, %msNhat2_13i
 *                    <br>%msNhat2_21i, %msNhat2_22i, %msNhat2_23i</td>
 *                    <br>%msNhat2_31i, %msNhat2_32i, %msNhat2_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(m_N^2)@f$</td>
 *   <td class="mod_desc">Imaginary part of the squared soft mass matrix
 *   of the down-type  @f$SU(2)_L@f$ singlet slepton</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%TUhat_11r, %TUhat_12r, %TUhat_13r
 *                    <br>%TUhat_21r, %TUhat_22r, %TUhat_23r</td>
 *                    <br>%TUhat_31r, %TUhat_32r, %TUhat_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(T_U)@f$</td>
 *   <td class="mod_desc">Real part of the tri-linear coupling of the up-type 
 *   squarks</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%TUhat_11i, %TUhat_12i, %TUhat_13i
 *                    <br>%TUhat_21i, %TUhat_22i, %TUhat_23i</td>
 *                    <br>%TUhat_31i, %TUhat_32i, %TUhat_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(T_U)@f$</td>
 *   <td class="mod_desc">Imaginary part of the tri-linear coupling of the up-type 
 *   squarks</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%TDhat_11r, %TDhat_12r, %TDhat_13r
 *                    <br>%TDhat_21r, %TDhat_22r, %TDhat_23r</td>
 *                    <br>%TDhat_31r, %TDhat_32r, %TDhat_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(T_D)@f$</td>
 *   <td class="mod_desc">Real part of the tri-linear coupling of the down-type 
 *   squarks</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%TDhat_11i, %TDhat_12i, %TDhat_13i
 *                    <br>%TDhat_21i, %TDhat_22i, %TDhat_23i</td>
 *                    <br>%TDhat_31i, %TDhat_32i, %TDhat_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(T_D)@f$</td>
 *   <td class="mod_desc">Imaginary part of the tri-linear coupling of the down-type 
 *   squarks</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%TEhat_11r, %TEhat_12r, %TEhat_13r
 *                    <br>%TEhat_21r, %TEhat_22r, %TEhat_23r</td>
 *                    <br>%TEhat_31r, %TEhat_32r, %TEhat_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(T_E)@f$</td>
 *   <td class="mod_desc">Real part of the tri-linear coupling of the up-type 
 *   slepton</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%TEhat_11i, %TEhat_12i, %TEhat_13i
 *                    <br>%TEhat_21i, %TEhat_22i, %TEhat_23i</td>
 *                    <br>%TEhat_31i, %TEhat_32i, %TEhat_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(T_E)@f$</td>
 *   <td class="mod_desc">Imaginary part of the tri-linear coupling of the up-type 
 *   slepton</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%TNhat_11r, %TNhat_12r, %TNhat_13r
 *                    <br>%TNhat_21r, %TNhat_22r, %TNhat_23r</td>
 *                    <br>%TNhat_31r, %TNhat_32r, %TNhat_33r</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{R}(T_N)@f$</td>
 *   <td class="mod_desc">Real part of the tri-linear coupling of the down-type 
 *   slepton</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%TNhat_11i, %TNhat_12i, %TNhat_13i
 *                    <br>%TNhat_21i, %TNhat_22i, %TNhat_23i</td>
 *                    <br>%TNhat_31i, %TNhat_32i, %TNhat_33i</td>
 *   <td class="mod_symb">Elements of the matrix @f$\mathcal{I}(T_N)@f$</td>
 *   <td class="mod_desc">Imaginary part of the tri-linear coupling of the down-type 
 *   slepton</td>
 * </tr>
 * </table>
 * 
 * 
 * 
 */
class GeneralSUSY : public SUSY {
public:
    static const int NGeneralSUSYvars = 126;
    static const std::string GeneralSUSYvars[NGeneralSUSYvars];

    /**
     * @brief Friend classes of GeneralSUSY class.
     */
    GeneralSUSY();

    /**
     * @brief A method to initialize the model.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers defined in the current class.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief Initializes the %GeneralSUSY parameters found in the argument.
     * @param[in] DPars a map containing the parameters (all as double) to be used in Monte Carlo
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The pre-update method for %GeneralSUSY
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();
    
    /**
     * @brief The update method for %GeneralSUSY.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The post-update method for %GeneralSUSY.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %GeneralSUSY
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

private:
    /**
     * @brief A method to set the value of a parameter of %GeneralSUSY.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string, const double&);
    virtual void SetSoftTerms();
    double msQhat2_11r, msQhat2_12r, msQhat2_12i, msQhat2_13r, msQhat2_13i, msQhat2_22r, msQhat2_23r, msQhat2_23i, msQhat2_33r,
           msUhat2_11r, msUhat2_12r, msUhat2_12i, msUhat2_13r, msUhat2_13i, msUhat2_22r, msUhat2_23r, msUhat2_23i, msUhat2_33r,
           msDhat2_11r, msDhat2_12r, msDhat2_12i, msDhat2_13r, msDhat2_13i, msDhat2_22r, msDhat2_23r, msDhat2_23i, msDhat2_33r,
           msLhat2_11r, msLhat2_12r, msLhat2_12i, msLhat2_13r, msLhat2_13i, msLhat2_22r, msLhat2_23r, msLhat2_23i, msLhat2_33r,
           msEhat2_11r, msEhat2_12r, msEhat2_12i, msEhat2_13r, msEhat2_13i, msEhat2_22r, msEhat2_23r, msEhat2_23i, msEhat2_33r,
           msNhat2_11r, msNhat2_12r, msNhat2_12i, msNhat2_13r, msNhat2_13i, msNhat2_22r, msNhat2_23r, msNhat2_23i, msNhat2_33r,
           TUhat_11r, TUhat_12r, TUhat_13r, TUhat_21r, TUhat_22r, TUhat_23r, TUhat_31r, TUhat_32r, TUhat_33r,
           TUhat_11i, TUhat_12i, TUhat_13i, TUhat_21i, TUhat_22i, TUhat_23i, TUhat_31i, TUhat_32i, TUhat_33i,
           TDhat_11r, TDhat_12r, TDhat_13r, TDhat_21r, TDhat_22r, TDhat_23r, TDhat_31r, TDhat_32r, TDhat_33r,
           TDhat_11i, TDhat_12i, TDhat_13i, TDhat_21i, TDhat_22i, TDhat_23i, TDhat_31i, TDhat_32i, TDhat_33i,
           TEhat_11r, TEhat_12r, TEhat_13r, TEhat_21r, TEhat_22r, TEhat_23r, TEhat_31r, TEhat_32r, TEhat_33r,
           TEhat_11i, TEhat_12i, TEhat_13i, TEhat_21i, TEhat_22i, TEhat_23i, TEhat_31i, TEhat_32i, TEhat_33i,
           TNhat_11r, TNhat_12r, TNhat_13r, TNhat_21r, TNhat_22r, TNhat_23r, TNhat_31r, TNhat_32r, TNhat_33r,
           TNhat_11i, TNhat_12i, TNhat_13i, TNhat_21i, TNhat_22i, TNhat_23i, TNhat_31i, TNhat_32i, TNhat_33i;
};

/**
 * @}
 */

#endif	/* GeneralSUSY_H */
