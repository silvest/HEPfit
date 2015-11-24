/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CORRELATEDGAUSSIANOBSERVABLES_H
#define	CORRELATEDGAUSSIANOBSERVABLES_H

#include "Observable.h"
#include <boost/ptr_container/ptr_vector.hpp>

class ThObsFactory;

/**
 * @addtogroup Observable
 * @brief A module for model parameters and observables.
 * @{
 */

/**
 * @class CorrelatedGaussianObservables
 * @ingroup Observable
 * @brief A class for correlated Gaussian observables. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class builds the correlated Gaussian observables that
 * are specified in the SomeModel.conf file or specified by the user.
 */
class CorrelatedGaussianObservables {
public:
    
    /**
     * @brief Constructor.
     * @param[in] name_i a given name for the set of correlated Gaussian observables
     */
    CorrelatedGaussianObservables(std::string name_i);
    
    /**
     * @brief The default Constructor.
     */
    CorrelatedGaussianObservables();
    
    /**
     * @brief The copy constructor.
     */
    CorrelatedGaussianObservables(const CorrelatedGaussianObservables& orig);
    
    /**
     * @brief The default destructor.
     */
    virtual ~CorrelatedGaussianObservables();

    /**
     * @brief Computes the covariance matrix for the correlated Gaussian observables set.
     * @param Corr the correlation matrix for the correlated Gassian observables set
     */
    void ComputeCov(gslpp::matrix<double> Corr);

    /**
     * @brief A method to compute the weight associated with the observable.
     */
    virtual double computeWeight();

    /**
     * @brief A method to add observables to the list of correlated Gaussian observables.
     * @param Obs_i reference to an object of type Observable
     */
    void AddObs(Observable& Obs_i);

    /**
     * @brief A get method to access the vector of observables that are defined in
     * one correlated Gaussian observables set.
     * @return a vector of type Observable()
     */
    std::vector<Observable> getObs() const
    {
        return Obs;
    }
    
    /**
     * @brief A get method to access an element of the vector of observables that are defined in
     * one correlated Gaussian observables set.
     * @return an element of the vector of type Observable()
     */
    Observable getObs(int i) const
    {
        return (Obs.at(i));
    }

    /**
     * @brief A get method to access the name of the correlated Gaussian observables set.
     * @return the name
     */
    std::string getName() const
    {
        return name;
    }

    /**
     * @brief A get method to access the covariance matrix of the correlated Gaussian observables.
     */
    gslpp::matrix<double> getCov() const
    {
        return *Cov;
    }
    
    /**
     * @brief The parser for CorrelatedGaussianObservables.
     * @param[in] Observables the pointer vector containing the Observables
     * @param[in] ifile the stream containing the config file to be parsed
     * @param[in] beg the iterator that parses a line in the config file
     * @param[in] infilename the name of the config file being parsed
     * @param[in] myModel a pointer to the model
     * @param[in] lineNo the current line number at which the file is being parsed
     * @param[in] rank the rank of the process that is using the parser
     * @return the line number (integer) after the parsing is done
     */
    int ParseCGO(boost::ptr_vector<Observable>& Observables, 
                 std::ifstream& ifile, 
                 boost::tokenizer<boost::char_separator<char> >::iterator& beg, 
                 std::string& infilename, 
                 ThObsFactory& myObsFactory,
                 StandardModel * myModel,
                 int lineNo,
                 int rank);
    /**
     * @brief A method to check if the end of file has been reached
     * @return a boolean which is true if the end of file has been reached
     */
    bool isEOF()
    {
        return IsEOF;
    }
    
private:
    std::vector<Observable> Obs;///< A vector of observables whose correlation will be calculated.
    gslpp::matrix<double>* Cov;///< The covariance matrix.
    std::string name;///< The name of the correlated Gaussian Observables set.
    std::string filepath;///< The path to the config file being parsed
    bool IsEOF;///< A boolean which is true if the end of file is reached.
};

/** 
 * @}
 */

#endif	/* CORRELATEDGAUSSIANOBSERVABLES_H */

