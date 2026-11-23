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
 *   <td class="mod_name">%reDC9, %reDC10</td>
 *   <td class="mod_symb">@f$ \mathcal{R}(\Delta C_i), (i=9,10)  @f$</td>
 *   <td class="mod_desc">The real part of @f$ \Delta C_i, (i=9,10) @f$, the semileptonic operators.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%imDC9, %imDC10</td>
 *   <td class="mod_symb">@f$ \mathcal{I}(\Delta C_i), (i=9,10)  @f$</td>
 *   <td class="mod_desc">The imaginary part of @f$ \Delta C_i, (i=9,10) @f$, the semileptonic operators.</td>
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
 *   <td class="mod_name">%reDC9p, %reDC10p</td>
 *   <td class="mod_symb">@f$ \mathcal{R}(\Delta C_i^\prime), (i=9,10)  @f$</td>
 *   <td class="mod_desc">The real part of @f$ \Delta C_i^\prime, (i=9,10) @f$
 *   which are the chirality flipped primed operators.</td>
 * </tr>
 *  <tr>
 *   <td class="mod_name">%imDC9p, %imDC10p</td>
 *   <td class="mod_symb">@f$ \mathcal{I}(\Delta C_i^\prime), (i=9,10)  @f$</td>
 *   <td class="mod_desc">The imaginary part of @f$ \Delta C_i^\prime, (i=9,10) @f$
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
 * The shifts enter the @f$ b\to s\gamma @f$ (CMbsg, CMprimebsg), @f$ b\to s\ell\ell @f$ (CMBMll, CMprimeBMll) and
 * @f$ B_q\to\mu\mu @f$ (CMbsmm, CMbdmm, C10 only) matchings. The former parameters DC1-DC8, DC7p and DC8p (documented as
 * current-current and electroweak-penguin shifts) had no counterpart in these operator bases and never entered any
 * observable; they have been removed.
 */
class FlavourWilsonCoefficient: public StandardModel {
public:

    static const int NFlavourWilsonCoefficientvars = 17;

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
     * @brief A getter for the complex Wilson coefficient \f$\Delta C'_9\f$
     * @return gslpp::complex \f$\Delta C'_9\f$
     */
    gslpp::complex getDC9p() const {
        return DC9p;
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

    
    gslpp::complex DC9, DC10, DC7g, DC8g;///< The complex Wilson Coefficients.
    gslpp::complex DC9p, DC10p, DC7gp, DC8gp;///< The complex Wilson Coefficients.
    double reDC9, reDC10, reDC7g, reDC8g;///< The real parts of the  Wilson Coefficients.
    double imDC9, imDC10, imDC7g, imDC8g;///< The imaginary parts of the  Wilson Coefficients.
    double reDC9p, reDC10p, reDC7gp, reDC8gp;///< The real parts of the  Wilson Coefficients.
    double imDC9p, imDC10p, imDC7gp, imDC8gp;///< The imaginary parts of the  Wilson Coefficients.
    double WCscale;///< The scale of the  Wilson Coefficients.
    
      
};

/**
 * @}
 */

#endif	/* FLAVOURWILSONCOEFFICIENT_H */

