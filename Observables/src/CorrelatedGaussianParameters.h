/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CORRELATEDGAUSSIANPARAMETERS_H
#define	CORRELATEDGAUSSIANPARAMETERS_H

#include "ModelParameter.h"
#include "gslpp.h"

/**
 * @addtogroup Observables
 * @brief A module for model parameters and observables.
 * @{
 */

/**
 * @class CorrelatedGaussianParameters
 * @ingroup Observables
 * @brief A class for correlated Gaussian parameters. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class builds the correlated Gaussian parameters that
 * are specified in the SomeModel.conf file or specified by the user.
 */
class CorrelatedGaussianParameters {
public:

    /**
     * @brief Constructor.
     * @param[in] name_i a given name for the set of correlated Gaussian parameters
     */
    CorrelatedGaussianParameters(std::string name_i);

    /**
     * @brief The copy constructor.
     */
    CorrelatedGaussianParameters(const CorrelatedGaussianParameters& orig);

    /**
     * @brief The default destructor.
     */
    virtual ~CorrelatedGaussianParameters();

    /**
     * @brief Diagonalizes the correlated Gaussian parameters set.
     * @param Corr the correlation matrix for the correlated Gaussian parameters set
     */
    void DiagonalizePars(gslpp::matrix<double> Corr);

    /**
     * @brief A method to add parameters to the list of correlated Gaussian parameters.
     * @param Par_i reference to an object of type ModelParameter
     */
    void AddPar(ModelParameter& Par_i);

    /**
     * @brief A get method to access the vector of parameters that are defined in
     * one correlated Gaussian parameters set.
     * @return a vector of type ModelParameter()
     */
    std::vector<ModelParameter> getPars() const
    {
        return Pars;
    }

    /**
     * @brief A get method to access an element of the vector of parameters that are defined in
     * one correlated Gaussian parameters set.
     * @return an element of the vector of type ModelParameter()
     */
    ModelParameter getPar(int i) const
    {
        return (Pars.at(i));
    }

    /**
     * @brief A get method to access the name of the correlated Gaussian parameters set.
     * @return the name
     */
    std::string getName() const
    {
        return name;
    }

    /**
     * @brief A set method to fix the name of the set of correlated Gaussian parameters.
     * @param name the name
     */
    void setName(std::string name)
    {
        this->name = name;
    }

    /**
     * @brief A get method to access the covariance matrix of the correlated Gaussian parameters.
     */
    gslpp::matrix<double> getCov() const
    {
        return *Cov;
    }

    std::vector<ModelParameter> getDiagPars() const
    {
        return DiagPars;
    }

    std::vector<double> getOrigParsValue(const std::vector<double>& DiagPars_i) const;

private:
    std::vector<ModelParameter> Pars; ///< A vector of parameters whose correlation will be calculated.
    gslpp::matrix<double>* Cov = NULL; ///< The covariance matrix.
    std::string name; ///< The name of the correlated Gaussian Parameters set.
    gslpp::matrix<double> * v = NULL;
    gslpp::vector<double> * e = NULL;
    std::vector<ModelParameter> DiagPars; ///< The vector of diagonal parameters
};

/** 
 * @}
 */

#endif	/* CORRELATEDGAUSSIANPARAMETERS_H */

