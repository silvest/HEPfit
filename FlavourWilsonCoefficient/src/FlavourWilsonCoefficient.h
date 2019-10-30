/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOURWILSONCOEFFICIENT_H
#define	FLAVOURWILSONCOEFFICIENT_H

#include "StandardModel.h"
#include "gslpp.h"
#include "FlavourWilsonCoefficientMatching.h"

/**
 * @class FlavourWilsonCoefficient
 * @ingroup FlavourWilsonCoefficent
 * @brief %Model for NP contributions to flavour
 * through shifts to Standard %Model Wilson coefficients.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 * 
 *
 * @anchor FlavourWilsonCoefficientParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %FlavourWilsonCoefficient model are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%reDC1 -  %reDC10</td>
 *   <td class="mod_symb">@f$ \mathcal{R}(\Delta C_i), (i=1,\ldots,10)  @f$</td>
 *   <td class="mod_desc">The real part of @f$ \Delta C_i, (i=1,\ldots,10) @f$.
 *   (1-4) are current-current operators, (5-8) are electroweak penguin operators
 *   (9-10) are semileptonic operators.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%imDC1 -  %imDC10</td>
 *   <td class="mod_symb">@f$ \mathcal{I}(\Delta C_i), (i=1,\ldots,10)  @f$</td>
 *   <td class="mod_desc">The imaginary part of @f$ \Delta C_i, (i=1,\ldots,10) @f$.
 *   (1-4) are current-current operators, (5-8) are electroweak penguin operators
 *   (9-10) are semileptonic operators.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%reDC7g -  %reDC8g</td>
 *   <td class="mod_symb">@f$ \mathcal{R}(\Delta C_i), (i=7g, 8g)  @f$</td>
 *   <td class="mod_desc">The real part of @f$ \Delta C_i, (i=7g,8g) @f$,
 *   7g being the electric dipole operator and 8g being the chromomagnetic 
 *   dipole operator.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%imDC7g -  %imDC8g</td>
 *   <td class="mod_symb">@f$ \mathcal{I}(\Delta C_i), (i=7g, 8g)  @f$</td>
 *   <td class="mod_desc">The imaginary part of @f$ \Delta C_i, (i=7g,8g) @f$,
 *   7g being the electric dipole operator and 8g being the chromomagnetic 
 *   dipole operator.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%reDC7p -  %reDC10p</td>
 *   <td class="mod_symb">@f$ \mathcal{R}(\Delta C_i^\prime), (i=7,\ldots,10)  @f$</td>
 *   <td class="mod_desc">The real part of @f$ \Delta C_i^\prime, (i=1,\ldots,10) @f$
 *   which are the chirality flipped primed operators.</td>
 * </tr>
 *  <tr>
 *   <td class="mod_name">%imDC7p -  %imDC10p</td>
 *   <td class="mod_symb">@f$ \mathcal{I}(\Delta C_i^\prime), (i=7,\ldots,10)  @f$</td>
 *   <td class="mod_desc">The imaginary part of @f$ \Delta C_i^\prime, (i=1,\ldots,10) @f$
 *   which are the chirality flipped primed operators.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%reDC7gp -  %reDC8gp</td>
 *   <td class="mod_symb">@f$ \mathcal{R}(\Delta C_i^\prime), (i=7g, 8g)  @f$</td>
 *   <td class="mod_desc">The real part of @f$ \Delta C_i^\prime, (i=7g,8g) @f$
 *   which are the chirality flipped primed operators.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%imDC7gp -  %imDC8gp</td>
 *   <td class="mod_symb">@f$ \mathcal{I}(\Delta C_i^\prime), (i=7g, 8g)  @f$</td>
 *   <td class="mod_desc">The imaginary part of @f$ \Delta C_i^\prime, (i=7g,8g) @f$,
 *   which are the chirality flipped primed operators.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%WCscale</td>
 *   <td class="mod_symb">@f$ \mu_{matching} @f$</td>
 *   <td class="mod_desc">The Matching scale of the Wilson coefficients</td>
 * </tr>
 * </table>
 */
class FlavourWilsonCoefficient: public StandardModel {
public:

    static const int NFlavourWilsonCoefficientvars = 37;

    static const std::string FlavourWilsonCoefficientvars[NFlavourWilsonCoefficientvars];
    
    /**
     * @brief %FlavourWilsonCoefficient constructor
     */
    FlavourWilsonCoefficient();
    
    /**
     * @brief %FlavourWilsonCoefficient destructor
     */
    ~FlavourWilsonCoefficient();
    
    /**
     * @brief A method to initialize the model.
     * @details This method, called via InputParser::ReadParameters(), allocates
     * memory to the pointers defined in the current class.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief Initializes the %FlavourWilsonCoefficient parameters found in the argument.
     * @param[in] DPars a map containing the parameters (all as double) to be used in Monte Carlo
     */
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The pre-update method for %FlavourWilsonCoefficient
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PreUpdate();    
    
    /**
     * @brief The update method for %FlavourWilsonCoefficient.
     * @details This method updates all the model parameters with given DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief The post-update method for %FlavourWilsonCoefficient.
     * @details This method runs all the procedures that are need to be executed
     * after the model is successfully updated.
     * @return a boolean that is true if the execution is successful
     */
    virtual bool PostUpdate();
    
    /**
     * @brief A method to check if all the mandatory parameters for %FlavourWilsonCoefficient
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to set a flag of %FlavourWilsonCoefficient.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     * @brief A get method to access the member reference of type %FlavourWilsonCoefficientMatching.
     * @return a reference to a %FlavourWilsonCoefficientMatching object
     */
    virtual FlavourWilsonCoefficientMatching& getMatching() const
    {
        return FWCM.getObj();
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_1\f$
     * @return gslpp::complex \f$\Delta C_1\f$
     */
    gslpp::complex getDC1() const {
        return DC1;
    }
    
    /**
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_2\f$
     * @return gslpp::complex \f$\Delta C_2\f$
     */
    gslpp::complex getDC2() const {
        return DC2;
    }
    
    /**
     * 
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_3\f$
     * @return gslpp::complex \f$\Delta C_3\f$
     */
    gslpp::complex getDC3() const {
        return DC3;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_4\f$
     * @return gslpp::complex \f$\Delta C_4\f$
     */
    gslpp::complex getDC4() const {
        return DC4;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_5\f$
     * @return gslpp::complex \f$\Delta C_5\f$
     */
    gslpp::complex getDC5() const {
        return DC5;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_6\f$
     * @return gslpp::complex \f$\Delta C_6\f$
     */
    gslpp::complex getDC6() const {
        return DC6;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_7\f$
     * @return gslpp::complex \f$\Delta C_7\f$
     */
    gslpp::complex getDC7() const {
        return DC7;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_8\f$
     * @return gslpp::complex \f$\Delta C_8\f$
     */
    gslpp::complex getDC8() const {
        return DC8;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_9\f$
     * @return gslpp::complex \f$\Delta C_9\f$
     */
    gslpp::complex getDC9() const {
        return DC9;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_10\f$
     * @return gslpp::complex \f$\Delta C_10\f$
     */
    gslpp::complex getDC10() const {
        return DC10;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_{7\gamma}\f$
     * @return gslpp::complex \f$\Delta C_{7\gamma}\f$
     */
    gslpp::complex getDC7g() const {
        return DC7g;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C_{8g}\f$
     * @return gslpp::complex \f$\Delta C_{8g}\f$
     */
    gslpp::complex getDC8g() const {
        return DC8g;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C'_7\f$
     * @return gslpp::complex \f$\Delta C'_7\f$
     */
    gslpp::complex getDC7p() const {
        return DC7p;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C'_8\f$
     * @return gslpp::complex \f$\Delta C'_8\f$
     */
    gslpp::complex getDC8p() const {
        return DC8p;
    }
    
    /**
     * 
     * @brief A getter for the complex Wilson coefficient \f$\Delta C'_9\f$
     * @return gslpp::complex \f$\Delta C'_9\f$
     */
    gslpp::complex getDC9p() const {
        return DC9;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C'_10\f$
     * @return gslpp::complex \f$\Delta C'_10\f$
     */
    gslpp::complex getDC10p() const {
        return DC10p;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C'_{7\gamma}\f$
     * @return gslpp::complex \f$\Delta C'_{7\gamma}\f$
     */
    gslpp::complex getDC7gp() const {
        return DC7gp;
    }
    
    /**
     *
     * @brief A getter for the complex Wilson coefficient \f$\Delta C'_{8g}\f$
     * @return gslpp::complex \f$\Delta C'_{8g}\f$
     */
    gslpp::complex getDC8gp() const {
        return DC8gp;
    }
    
    /**
     *
     * @brief A getter for the scale of the Wilson coefficients
     * @return the matching scale of the Wilson coefficients
     */
    double getWCscale() const {
        return WCscale;
    }
    
protected: 
    
    /**
     * @brief A method to set the value of a parameter of %FlavourWilsonCoefficient.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string, const double&);
    mutable Matching<FlavourWilsonCoefficientMatching,FlavourWilsonCoefficient> FWCM; ///< The FlavourWilsonCoefficientMatching object.

private:

    
    gslpp::complex DC1, DC2, DC3, DC4, DC5, DC6, DC7, DC8, DC9, DC10, DC7g, DC8g;///< The complex Wilson Coefficients.
    gslpp::complex DC7p, DC8p, DC9p, DC10p, DC7gp, DC8gp;///< The complex Wilson Coefficients.
    double reDC1, reDC2, reDC3, reDC4, reDC5, reDC6, reDC7, reDC8, reDC9, reDC10, reDC7g, reDC8g;///< The real parts of the  Wilson Coefficients.
    double imDC1, imDC2, imDC3, imDC4, imDC5, imDC6, imDC7, imDC8, imDC9, imDC10, imDC7g, imDC8g;///< The imaginary parts of the  Wilson Coefficients.
    double reDC7p, reDC8p, reDC9p, reDC10p, reDC7gp, reDC8gp;///< The real parts of the  Wilson Coefficients.
    double imDC7p, imDC8p, imDC9p, imDC10p, imDC7gp, imDC8gp;///< The imaginary parts of the  Wilson Coefficients.
    double WCscale;///< The scale of the  Wilson Coefficients.
    
      
};

/**
 * @}
 */

#endif	/* FLAVOURWILSONCOEFFICIENT_H */

