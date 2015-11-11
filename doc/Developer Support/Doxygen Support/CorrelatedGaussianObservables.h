/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DOXYEXAMPLE_H
#define	DOXYEXAMPLE_H

#include "Observable.h"

/**
 * @addtogroup Documentation
 * @brief A module for setting examples for documentation.
 * @{
 */

/**
 * @class DoxyExample
 * @ingroup Documentation
 * @brief A class for setting an example for Doxygen documentation.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class shows how different members of a class and the 
 * class itself can be documneted in Doxygen.
 */
class DoxyExample {
public:
    
    /**
     * @brief The default constructor.
     * @param[in] name_i a given name for the set of Doxygen observables
     */
    DoxyExample(std::string name_i);
    
    /**
     * @brief The copy constructor.
     */
    DoxyExample(const DoxyExample& orig);
    
    /**
     * @brief The default destructor.
     */
    virtual ~DoxyExample();

    /**
     * @brief Computes the covariance matrix for the Doxygen set.
     * @param Corr the correlation matrix for the Doxygen set
     */
    void ComputeCov(gslpp::matrix<double> Corr);

    /**
     * @brief A method to add observables to the list of Doxygen.
     * @param obs_i reference to an object of type Observable()
     */
    void AddObs(Observable& Obs_i);

    /**
     * @brief A get method to access the vector of observables that are defined in
     * one Doxygen set.
     * @return a vector of type Observable()
     */
    std::vector<Observable> getObs() const
    {
        return Obs;
    }

    /**
     * @brief A get method to access the name of the Doxygen set.
     * @return the name
     */
    std::string getName() const
    {
        return name;
    }

    /**
     * @brief A set method to fix the name of the set of Doxygen.
     * @param name the name
     */
    void setName(std::string name)
    {
        this->name = name;
    }

    /**
     * @brief A get method to access the covariance matrix of the Doxygen.
     * @param the covariance matrix
     */
    gslpp::matrix<double> getCov() const
    {
        return *Cov;
    }
    
private:
    std::vector<Observable> Obs;/**< A vector of observables whose correlation will be calculated. */
    gslpp::matrix<double>* Cov;/**< The covariance matrix. */
    std::string name;/**< The name of the correlated Gaussian Observables set. */
};

/** 
 * @}
 */

#endif	/* DOXYEXAMPLE_H */

