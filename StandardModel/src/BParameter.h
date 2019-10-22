/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BPARAMETER_H
#define	BPARAMETER_H

#include <gslpp_vector_double.h>
#include <vector>
#include <map>
#include <functional>
#include "OrderScheme.h"

/**
 * @class BParameter
 * @ingroup StandardModel
 * @brief A class for the bag parameters.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is the class for defining bag parameters, which depend on
 * a specified scale and scheme. 
 * 
 *
 * @anchor BParameterParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %BParameter are summarized below:
 * <table class="model">
 * <tr>
 *   <td class="mod_name">%BBs1 - %BBs5</td>
 *   <td class="mod_symb">@f$B^1_{B_s} - B^5_{B_s}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_5 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_s \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FBsSqrtBBs1 - %FBsSqrtBBs5</td>
 *   <td class="mod_symb">@f$f_{B_s}\sqrt{B^1_{B_s}} - f_{B_s}\sqrt{B^5_{B_s}}@f$</td>
 *   <td class="mod_desc">The decay constant times the square root of the bag parameter for \f$ O_1 - O_5 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_s \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBsscale</td>
 *   <td class="mod_symb">@f$\mu_{B_{s}}@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ B_s \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBsscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ B_s \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBsoBBd</td>
 *   <td class="mod_symb">@f$B^1_{B_s}/B^1_{B_d}@f$</td>
 *   <td class="mod_desc">The ratio \f$ B^1_{B_s}/B^1_{B_d} \f$ necessary to compute \f$ B^1_{B_d} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBd2 - %BBd5</td>
 *   <td class="mod_symb">@f$B^2_{B_d} - B^5_{B_d}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_2 - O_5 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_d \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%csi</td>
 *   <td class="mod_symb">@f$\frac{f_{B_s}\sqrt{B^1_{B_s}}}{f_{B_d}\sqrt{B^1_{B_d}}}@f$</td>
 *   <td class="mod_desc">The ratio \f$ \frac{f_{B_s}\sqrt{B^1_{B_s}}}{f_{B_d}\sqrt{B^1_{B_d}}} \f$ necessary to compute \f$ f_{B_d}\sqrt{B^1_{B_d}} \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%FBdSqrtBBd2 - %FBdSqrtBBd5</td>
 *   <td class="mod_symb">@f$f_{B_d}\sqrt{B^2_{B_d}} - f_{B_d}\sqrt{B^5_{B_d}}@f$</td>
 *   <td class="mod_desc">The decay constant times the square root of the bag parameter for \f$ O_2 - O_5 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_d \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBdscale</td>
 *   <td class="mod_symb">@f$\mu_{B_{d}}@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ B_d \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BBdscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ B_d \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BK1 - %BK5</td>
 *   <td class="mod_symb">@f$B^1_{K} - B^5_{K}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_5\f$ in \f$ \Delta s = 2 \f$ processes in \f$ K^0 \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKscale</td>
 *   <td class="mod_symb">@f$\mu_K@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ K^0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ K^0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BD1 - %BD5</td>
 *   <td class="mod_symb">@f$B^1_{D} - B^5_{D}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_5\f$ in \f$ \Delta c = 2 \f$ processes in \f$ D^0 \f$.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BDscale</td>
 *   <td class="mod_symb">@f$\mu_D@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ D_0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BDscheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ D_0 \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BK(1/2)1 - %BK(1/2)10</td>
 *   <td class="mod_symb">@f$B^1_{K(1/2)} - B^{10}_{K(1/2)}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_{10}\f$ in \f$ K_L \f$ decay in 2 pion with 0 isospin change.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BK(3/2)1 - %BK(3/2)10</td>
 *   <td class="mod_symb">@f$B^1_{K(3/2)} - B^{10}_{K(3/2)}@f$</td>
 *   <td class="mod_desc">The bag parameter for \f$ O_1 - O_{10}\f$ in \f$ K_L \f$ decay in 2 pion with double isospin change.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKd_scale</td>
 *   <td class="mod_symb">@f$\mu_{K_L}@f$</td>
 *   <td class="mod_desc">The scale at which the bag parameters are specified for the \f$ K_L \f$ system.</td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%BKd_scheme</td>
 *   <td class="mod_symb">@f$@f$</td>
 *   <td class="mod_desc">The scheme in which the bag parameters are specified for the \f$ K_L \f$ system.</td>
 * </tr>
 * </table>
 * 
 */
class BParameter {
public:
    
    /**
     * @brief Constructor.
     * @param[in] n dimension of the vector of bag parameters
     */
    BParameter(int n);
    
    /**
     * @brief Constructor.
     * @param[in] n dimension of the vector of bag parameters
     * @param[in] name_i name used to identify the B Parameter
     */
    BParameter(int n, std::string name_i);
    
    /**
     * @brief The default destructor.
     */
    virtual ~BParameter();

    /**
     * @brief A get method for the vector of the bag parameters.
     * @return the vector of the bag parameters
     */
    const gslpp::vector<double>& getBpars() const
    {
        return bpars;
    }

    /**
     * @brief A set method for a vector of the bag parameters.
     * @param[in] bpars a vector of the bag parameters
     */
    void setBpars(gslpp::vector<double> bpars) 
    {
        this->bpars = bpars;
    }

    /**
     * @brief A set method for a component of the vector of bag parameters.
     * @param[in] i the index for the component of the vector of bag parameters
     * @param[in] value the value of the bag parameters
     */
    void setBpars(int i, double value) 
    {
        this->bpars(i) = value;
    }

    /**
     * @brief A get method for the scale of the bag parameters.
     * @return the scale at which the bag parameters are defined
     */
    const double& getMu() const
    {
        return mu;
    }

    /**
     * @brief A set method for the scale of the bag parameters.
     * @param[in] mu the scale at which the bag parameters are defined
     */
    void setMu(double mu) 
    {
        this->mu = mu;
    }
    
    /**
     * @brief A get method for the scheme of the bag parameters.
     * @return the scheme in which the bag parameters are defined
     */
    schemes getScheme() const
    {
        return scheme;
    }
    
    /**
     * @brief A set method for the scheme of the bag parameters.
     * @param[in] scheme the scheme in which the bag parameters are defined
     */
    void setScheme(schemes scheme) 
    {
        this->scheme = scheme;
    }
    
    /**
     * @brief A method to get the parameters list.
     * @param[in] name the name of the bag parameter.
     */
    std::vector<std::string> parameterList(std::string name_i);
    
    /**
     * @brief A set method for setting the parameters.
     * @param[in] name the name of the parameter to be set
     * @param[in] value the value of the parameter to be set
     */
    bool setParameter(std::string name_i, double value);
    
    /**
     * @brief A method to pass the list of parameters as observables.
     * @param[in] ModelParamMap the reference to the Model Parameter Map.
     */
    void ModelParameterMapInsert(std::map< std::string, std::reference_wrapper<const double> >& ModelParamMap);
    
    /**
     * @brief A set method to set the flag for the parameter set to be used for BBs and BBd.
     * @param[in] flag the boolean value of the flag
     */
    void setFlagCsi(bool flag) {
        FlagCsi = flag;
    }
    
    /**
     * @brief A set method for BBs/BBd.
     * @param[in] BBsoBBd the parameter BBs/BBd
     */
    void setBBsoBBd(double BBsoBBd) 
    {
        this->BBsoBBd = BBsoBBd;
    }
       
    /**
     * @brief A get method for BBs/BBd.
     * @return the parameter BBs/BBd.
     */
    const double& getBBsoBBd() const
    {
        if (name.compare("BBd") == 0) return BBsoBBd;
        else throw std::runtime_error("BBsoBBd belongs to " + name);
    }
    
    /**
     * @brief A get method for csi.
     * @return the parameter csi.
     */
    const double& getcsi() const
    {
        if (name.compare("BBd") == 0) return csi;
        else throw std::runtime_error("csi belongs to " + name);
    }
    
    /**
     * @brief A get method for FBsSqrtBBs1.
     * @return the parameter FBsSqrtBBs1.
     */
    const double& getFBsSqrtBBs1() const
    {
        if (name.compare("BBs") == 0) return FBsSqrtBBs1;
        else throw std::runtime_error("FBsSqrtBBs1 belongs to " + name);
    }
    
    /**
     * @brief A get method for FBsSqrtBBs2.
     * @return the parameter FBsSqrtBBs2.
     */
    const double& getFBsSqrtBBs2() const
    {
        if (name.compare("BBs") == 0) return FBsSqrtBBs2;
        else throw std::runtime_error("FBsSqrtBBs2 belongs to " + name);
    }
    
    /**
     * @brief A get method for FBsSqrtBBs3.
     * @return the parameter FBsSqrtBBs3.
     */
    const double& getFBsSqrtBBs3() const
    {
        if (name.compare("BBs") == 0) return FBsSqrtBBs3;
        else throw std::runtime_error("FBsSqrtBBs3 belongs to " + name);
    }
    
    /**
     * @brief A get method for FBsSqrtBBs4.
     * @return the parameter FBsSqrtBBs1.
     */
    const double& getFBsSqrtBBs4() const
    {
        if (name.compare("BBs") == 0) return FBsSqrtBBs4;
        else throw std::runtime_error("FBsSqrtBBs4 belongs to " + name);
    }
    
    /**
     * @brief A get method for FBsSqrtBBs5.
     * @return the parameter FBsSqrtBBs5.
     */
    const double& getFBsSqrtBBs5() const
    {
        if (name.compare("BBs") == 0) return FBsSqrtBBs5;
        else throw std::runtime_error("FBsSqrtBBs5 belongs to " + name);
    }
    
    /**
     * @brief A get method for FBdSqrtBBd2.
     * @return the parameter FBdSqrtBBd2.
     */
    const double& getFBdSqrtBBd2() const
    {
        if (name.compare("BBd") == 0) return FBdSqrtBBd2;
        else throw std::runtime_error("FBdSqrtBBd2 belongs to " + name);
    }
    
    /**
     * @brief A get method for FBdSqrtBBd3.
     * @return the parameter FBdSqrtBBd3.
     */
    const double& getFBdSqrtBBd3() const
    {
        if (name.compare("BBd") == 0) return FBdSqrtBBd3;
        else throw std::runtime_error("FBdSqrtBBd3 belongs to " + name);
    }
    
    /**
     * @brief A get method for FBdSqrtBBd4.
     * @return the parameter FBdSqrtBBd4.
     */
    const double& getFBdSqrtBBd4() const
    {
        if (name.compare("BBd") == 0) return FBdSqrtBBd4;
        else throw std::runtime_error("FBdSqrtBBd4 belongs to " + name);
    }
    
    /**
     * @brief A get method for FBdSqrtBBd5.
     * @return the parameter FBdSqrtBBd5.
     */
    const double& getFBdSqrtBBd5() const
    {
        if (name.compare("BBd") == 0) return FBdSqrtBBd5;
        else throw std::runtime_error("FBdSqrtBBd5 belongs to " + name);
    }
    

private:
    gslpp::vector<double> bpars;///< A vector of bag parameters.
    double mu;///< The scale at which the bag parameters are defined. 
    schemes scheme;///< The scheme in which the bag parameters are defined.
    std::string name;   
    bool FlagCsi;
    double BBsoBBd;
    double csi;
    double FBsSqrtBBs1;
    double FBsSqrtBBs2;
    double FBsSqrtBBs3;
    double FBsSqrtBBs4;
    double FBsSqrtBBs5;
    double FBdSqrtBBd2;
    double FBdSqrtBBd3;
    double FBdSqrtBBd4;
    double FBdSqrtBBd5;
};

/**
 * @}
 */

#endif	/* BPARAMETER_H */
