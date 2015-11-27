/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MODELPARAMETER_H
#define	MODELPARAMETER_H

#include <string>
#include <iostream>
#include <boost/tokenizer.hpp>

/**
 * @class ModelParameter
 * @ingroup Observable
 * @brief A class for model parameters. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class stores the details of the model parameters defined in the
 * SomeModel.conf file or specified by the user. The name of the model parameter has to
 * correspond to the list of names in the SomeModel.cpp file in the
 * Model src directories.\n
 * e.g. GF, which can be found in StandardModel.cpp file in the
 * StandardModel project source folder StandardModel/src.
 */
class ModelParameter {
public:
    
    /**
     * @brief Constructor.
     * @param[in] name_in the name for the model parameter
     * @param[in] ave_in the average value for the model parameter
     * @param[in] errg_in the gaussian error for the model parameter
     * @param[in] errf_in the flat error for the model parameter
     */
    ModelParameter(std::string name_in, double ave_in, double errg_in, double errf_in);
    
    /**
     * @brief Constructor.
     */
    ModelParameter();
    
    /**
     * @brief Parser for model parameters
     * @param[in] beg iterator for reading the model parameter
     */
    boost::tokenizer<boost::char_separator<char> >::iterator & ParseModelParameter(boost::tokenizer<boost::char_separator<char> >::iterator & beg);
    
    /**
     * @brief The default destructor.
     */
    virtual ~ModelParameter();
    
     /**
     * @brief A method to check if the parameters are correlated or not
     * @return a boolean which is true of the parameters are correlated
     */
    bool IsCorrelated() const
    {
        return(!cgp_name.empty());
    }
    
    /**
     * @brief A get method to get the name of each parameter
     * @return a string that contains the name
     */
    std::string getname() const
    {
        return name;
    }
    
    /**
     * @brief A get method to get the average
     * @return the average of the parameter
     */
    double getave() const
    {
        return ave;
    }
    
    /**
     * @brief A get method to get the flat error
     * @return the flat error of the parameter
     */
    double geterrf() const
    {
        return errf;
    }
    
    /**
     * @brief A get method to get the gaussian error
     * @return the Gaussian error of the parameter
     */
    double geterrg() const
    {
        return errg;
    }
    
    /**
     * @brief A get method to get the minimum value
     * @return the minimum value of the parameter
     */
    double getmin() const
    {
        return min;
    }
    
    /**
     * @brief A get method to get the maximum value
     * @return the maximum value of the parameter
     */
    double getmax() const
    {
        return max;
    }
    
    /**
     * @brief A method to check if the parameter is fixed
     * @return a boolean that is true if the parameter is fixed
     */
    bool IsFixed() const
    {
        return isFixed;
    }
    
    /**
     * @brief A get method to get the name of the set of correlated parameter
     * @return a string that contains the name
     */
    std::string getCgp_name() const
    {
        return cgp_name;
    }

    /**
     * @brief A set method to set the name of the set of correlated parameter
     */
    void setCgp_name(std::string cgp_name)
    {
        this->cgp_name = cgp_name;
    }
    
    /**
     * @brief Befriending of the std::ostream operator << to generate an
     * output stream for printing the model parameters.
     * @param[out] output the formatted output stream to print the model parameters
     * @param[in] m a reference to an object of type ModelParameter()
     */
    friend std::ostream& operator<<(std::ostream& output, const ModelParameter& m);
 
private:
    std::string name; ///< The name of the model parameter.
    double ave; ///< The average value of the model parameter.
    double errg; ///< The Gaussian error of the model parameter.
    double errf; ///< The flat error of the model parameter.
    double min; ///< The minimum value of the model parameter.
    double max; ///< The maximum value of the model parameter.
    bool isFixed; ///< A boolean flag that is true if the parameter is fixed
    std::string cgp_name;
};

#endif	/* MODELPARAMETER_H */

