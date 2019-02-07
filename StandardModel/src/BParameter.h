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
 * @addtogroup StandardModel
 * @brief A module for the Standard %Model.
 * @{
 */

/**
 * @class BParameter
 * @brief A class for the bag parameters.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This is the class for defining bag parameters, which depend on
 * a specified scale and scheme. 
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
