/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSTUVWXY_H
#define	NPSTUVWXY_H

#include <cmath>
#include "NPbase.h"

/**
 * @class NPSTUVWXY
 * @brief A class for new physics in the form of contributions to the extended oblique 
 * parameters \f$\hat{S},~\hat{T},~\hat{U},~V,~W,~X\f$ and \f$Y\f$. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the necessary functions to compute new physics 
 * tree-level corrections to electroweak precision observables, in the form of 
 * contributions to the extended set of oblique parameters \cite Barbieri:2004qk. These corrections
 * are parameterized in terms of the \f$\hat{S},~\hat{T},~\hat{U},~V,~W,~X\f$ and \f$Y\f$ 
 * contributions to \f$M_W\f$, and to \f$Z\f$-pole observables through the corrections 
 * to the different neutral-current effective couplings to leptons and quarks. The 
 * contributions to the later are implemented in the \b NPbase class. 
 * 
 * @anchor NPSTUVWXYParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of NPSTUVWXY are summarized below: 
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%myObliqueShat </td>
 *   <td class="mod_symb">\f$\hat{S}\f$</td>
 *   <td class="mod_desc">The oblique paramenter \f$\hat{S}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%myObliqueThat </td>
 *   <td class="mod_symb">\f$\hat{T} \f$</td>
 *   <td class="mod_desc">The oblique paramenter \f$\hat{T}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%myObliqueUhat </td>
 *   <td class="mod_symb">\f$\hat{U}\f$</td>
 *   <td class="mod_desc">The oblique paramenter \f$\hat{U}\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%myObliqueV </td>
 *   <td class="mod_symb">\f$V\f$</td>
 *   <td class="mod_desc">The oblique paramenter \f$V\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%myObliqueW </td>
 *   <td class="mod_symb">\f$W\f$</td>
 *   <td class="mod_desc">The oblique paramenter \f$W\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%myObliqueX </td>
 *   <td class="mod_symb">\f$X\f$</td>
 *   <td class="mod_desc">The oblique paramenter \f$X\f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%myObliqueY </td>
 *   <td class="mod_symb">\f$Y\f$</td>
 *   <td class="mod_desc">The oblique paramenter \f$Y\f$.</td>
 * </tr>
 * </table>
 *  
 */
class NPSTUVWXY : public NPbase {
public:
    /**
     * @brief The number of new physics parameters in the model.
     */
    static const int NSTUVWXYvars = 7;
    /**
     * @brief A string array with the names of the new physics parameters in the model.
     */
    static const std::string STUVWXYvars[NSTUVWXYvars];
    
    /**
     * @brief Constructor.
     */
    NPSTUVWXY();

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const 
    {
        return "NPSTUVWXY";
    }

    /**
     * @brief A method to initialize the model.
     * @return true is model initialization is successful
     */
    virtual bool InitializeModel();

    /**
     * @brief A method to initialize the model.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * @return true is model initialization is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars); 
    
    /**
     * @brief The update method for the model class.
     * @details This method updates all the parameters of the model every time a
     * new set of parameters is generated.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful.
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to check if all the mandatory parameters for the model have been
     * provided in the model configuration file.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A set method to fix the flags for the model.
     * @param[in] name the name of the flag
     * @param[in] value the value of the flag that can be true or false
     * @return a boolean to designate the success or failure of this procedure
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     * @brief A method to check the sanity of the set of flags.
     * @return true if the set of flags is sane.
     */
    virtual bool CheckFlags() const;
    
    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief The oblique paramenter \f$\hat{S}\f$.
     * @return the value of the oblique parameter \f$\hat{S}=\frac{\alpha}{4\sin^2{\theta_W}}S\f$
     */
    virtual double obliqueShat() const
    {
        return myObliqueShat;
    }

    /**
     * @brief The oblique paramenter \f$\hat{T}\f$.
     * @return the value of the oblique parameter \f$\hat{T}=\alpha T\f$
     */
    virtual double obliqueThat() const 
    {
        return myObliqueThat;
    }

    /**
     * @brief The oblique paramenter \f$\hat{U}\f$.
     * @return the value of the oblique parameter \f$\hat{U}=-\frac{\alpha}{4\sin^2{\theta_W}}U\f$
     */
    virtual double obliqueUhat() const 
    {
        return myObliqueUhat;
    }

    /**
     * @brief The oblique paramenter \f$V\f$.
     * @return the value of the oblique parameter \f$V\f$
     */
    virtual double obliqueV() const 
    {
        return myObliqueV;
    }

    /**
     * @brief The oblique paramenter \f$W\f$.
     * @return the value of the oblique parameter \f$W\f$
     */
    virtual double obliqueW() const 
    {
        return myObliqueW;
    }

    /**
     * @brief The oblique paramenter \f$X\f$.
     * @return the value of the oblique parameter \f$X\f$
     */
    virtual double obliqueX() const 
    {
        return myObliqueX;
    }

    /**
     * @brief The oblique paramenter \f$Y\f$.
     * @return the value of the oblique parameter \f$Y\f$
     */
    virtual double obliqueY() const 
    {
        return myObliqueY;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The parameter \f$\varepsilon_1\f$.
     * @return the value of the \f$\varepsilon_1\f$ parameter (SM plus new physics
     * corrections)
     */
    double epsilon1() const;

    /**
     * @brief The parameter \f$\varepsilon_2\f$.
     * @return the value of the \f$\varepsilon_2\f$ parameter (SM plus new physics
     * corrections)
     */
    double epsilon2() const;
    
    /**
     * @brief The parameter \f$\varepsilon_3\f$.
     * @return the value of the \f$\varepsilon_3\f$ parameter (SM plus new physics
     * corrections)
     */
    double epsilon3() const;
    
    /**
     * @brief The parameter \f$\varepsilon_b\f$.
     * @return the SM value of the \f$\varepsilon_b\f$ parameter
     */
    double epsilonb() const;

    
    ////////////////////////////////////////////////////////////////////////     
    // Combinations of the extended oblique parameters
    
    /**
     * @brief The combination of extended oblique paramenters \f$S'\f$.
     * @return the value of the oblique parameter \f$S'=\frac{4\sin^2{\theta_W}}{\alpha}\left(\hat{S}-W+\frac{X}{\sin{\theta_W} \cos{\theta_W}}-Y\right)\f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The combination of extended oblique paramenters \f$T'\f$.
     * @return the value of the oblique parameter \f$T'=\frac{1}{\alpha}\left(\hat{T}-W+2X\frac{s_W}{c_W}-Y\frac{\sin^2{\theta_W}}{\cos^2{\theta_W}}\right)\f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The combination of extended oblique paramenters \f$U'\f$.
     * @return the value of the oblique parameter \f$U'=-\frac{4\sin^2{\theta_W}}{\alpha}\left(\hat{U}-W+2X\frac{\sin{\theta_W}}{\cos{\theta_W}}-V\right)\f$
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The \f$W\f$ boson mass.
     * @return the \f$W\f$-boson mass in GeV
     */
    virtual double Mw() const;

    /**
     * @brief The (square of the) cosine of the weak angle \f$\cos^2{\theta_W}\f$.
     * @return the value of \f$\cos^2{\theta_W}\f$ in the On-mass-shell renormalization scheme,
     *  \f$\cos^2{\theta_W}=\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double cW2() const;

    /**
     * @brief The (square of the) sine of the weak angle \f$\sin^2{\theta_W}\f$.
     * @return the value of \f$\sin^2{\theta_W}\f$ in the On-mass-shell renormalization scheme,
     *  \f$\sin^2{\theta_W}=1-\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double sW2() const;

    /**
     * @brief The \f$W\f$ decay width \f$\Gamma_W\f$.
     * @return the total width of the \f$W\f$ boson in GeV
     */
    virtual double GammaW() const;
    
    
    ////////////////////////////////////////////////////////////////////////
protected:    
     /**
     * @brief A set method to fix the parameters of the model.
     * @param[in] name a string with the parameter name
     * @param[in] value the value to be asigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);
    double myObliqueShat;///< The oblique paramenter \f$\hat{S}\f$.
    double myObliqueThat;///< The oblique paramenter \f$\hat{T}\f$.
    double myObliqueUhat;///< The oblique paramenter \f$\hat{U}\f$.
    double myObliqueV;///< The oblique paramenter \f$V\f$.
    double myObliqueW;///< The oblique paramenter \f$W\f$.
    double myObliqueX;///< The oblique paramenter \f$X\f$.
    double myObliqueY;///< The oblique paramenter \f$Y\f$.

};

#endif	/* NPSTUVWXY_H */

