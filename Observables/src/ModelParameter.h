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
     * @brief The default destructor.
     */
    virtual ~ModelParameter();
    
    
    std::string name; ///< The name of the model parameter.
    double ave; ///< The average value of the model parameter.
    double errg; ///< The Gaussian error of the model parameter.
    double errf; ///< The flat error of the model parameter.
    double min; ///< The minimum value of the model parameter.
    double max; ///< The maximum value of the model parameter.
    
    /**
     * @brief Befriending of the std::ostream operator << to generate an
     * output stream for printing the model parameters.
     * @param[out] output the formatted output stream to print the model parameters
     * @param[in] m a reference to an object of type ModelParameter()
     */
    friend std::ostream& operator<<(std::ostream& output, const ModelParameter& m);
private:
};

#endif	/* MODELPARAMETER_H */

