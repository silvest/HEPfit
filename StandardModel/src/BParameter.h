/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BPARAMETER_H
#define	BPARAMETER_H

#include <gslpp_vector_double.h>
#include "OrderScheme.h"

using namespace gslpp;

/**
 * @addtogroup StandardModel
 * @brief A module for the Standard %Model.
 * @{
 */

/**
 * @class BParameter
 * @brief A class for the bag paramters.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details The bag parameters are input values read from the configuration files
 * SomeModel.conf. They depend on a specified scale and scheme. Both the scale and the
 * scheme have to be specified in the same SomeModel.conf file. These parameters 
 * are set by the QCD class
 */
class BParameter {
public:
    /**
     * @brief BParameter constructor
     * @param[in] n dimension of the vector of bag parameters
     */
    BParameter(int n) : bpars(n,0.)
    {};
    /**
     * @brief The get method for the vector of bag parameters.
     * @return The vector of bag parameters
     */
    vector<double> getBpars() const
    {
        return bpars;
    }
    /**
     * @brief The set method for a vector of bag parameters
     * @param[in] bpars a vector of bag parameters read as input from SomeModel.conf
     */
    void setBpars(vector<double> bpars) 
    {
        this->bpars = bpars;
    }
    /**
     * @brief The set method for a component of the vector of bag parameters.
     * @param[in] i the index for the component of the vector of bag parameters
     * @param[in] value the value of the bag parameter
     */
    void setBpars(int i, double value) 
    {
        this->bpars(i) = value;
    }
    /**
     * @brief The get method for the scale of the bag parameter is specified in the SomeModel.conf file.
     * @return the scale at which the bag parameter is defined.
     */
    double getMu() const
    {
        return mu;
    }
    /**
     * @brief The set method for the scale of the bag parameter is specified in the SomeModel.conf file.
     * @param[in] mu the scale mu at which the bag parameter is defined.
     */
    void setMu(double mu) 
    {
        this->mu = mu;
    }
    /**
     * @brief The get method for the scheme in whcih the bag parameter is specified in the SomeModel.conf file.
     * @return the scheme in which the bag parameter is defined.
     */
    schemes getScheme() const
    {
        return scheme;
    }
    /**
     * @brief The set method for the scheme in which the bag parameter is specified in the SomeModel.conf file.
     * @param[in] scheme the scheme in which the bag parameter is defined.
     */
    void setScheme(schemes scheme) 
    {
        this->scheme = scheme;
    }

private:
    vector<double> bpars;
    double mu;
    schemes scheme;
};

/**
 * @}
 */

#endif	/* BPARAMETER_H */