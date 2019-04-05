/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSOBSERVABLE_H
#define	HIGGSOBSERVABLE_H

#include "Observable.h"
#include <TMatrixD.h>
#include "gslpp.h"
#include <boost/tokenizer.hpp>

class ThObsFactory;
/**
 * @class HiggsObservable
 * @ingroup Observable
 * @brief A class for Higgs experimental analyses
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The class for building observables encoding Higgs experimental analyses, storing the 
 * parameters read from a file specified in the SomeModel.conf file or by the user. The names (thnames) of the observables have
 * to correspond to the allowed names of observables listed in the ThFactory class.
 */

class HiggsObservable : public Observable {
public:
    HiggsObservable(const Observable& Obs);
               
    HiggsObservable(const HiggsObservable& orig);

//    virtual ~HiggsObservable();
    
    // Read the necessary information from the config file. First row should 
    // contain the list of categories used in the analysis

    /**
     * @brief Set the parametric likelihood describing one Higgs decay channel from a config file.
     * @param filename the name of the config file
     * @param thObsV a vector of ThObservables containing the ratio of the production cross section for 
     * the specified categories in the model analyzed over the SM prediction 
     */
    virtual void setParametricLikelihood(std::string filename, std::vector<ThObservable*> thObsV);

    /**
     * @brief A method to compute the weight associated with the observable.
     */
    virtual double computeWeight();
    
    /**
     * @brief A method to set the observable to the new parametric form.
     */
    void setIsnew(bool isnew)
    {
        this->isnew = isnew;
    }
    
    /**
     * @brief A method to get the observable to the new parametric form.
     * @return a boolean which is true if the observable is in the new parametric form
     */
    bool isNew()
    {
        return isnew;
    }
    
    void getTheoryValues(std::vector<double>& theoryValues_i) 
    {
        theoryValues_i = theoryValues;
    }
    
    int getNTheoryValues() 
    {
        return thobsvsize;
    }
    
    int getNChannels() 
    {
        return channelsize;
    }
    
    bool getIsCorrelated()
    {
        return isCorrelated;
    }
    
    void setIsCorrelated(bool correlated)
    {
        isCorrelated = correlated;
    }
    
    
    /**
     * @brief the parser for HiggsObservables
     * @param[in] beg the iterator that parses a line in the config file
     * @param[in] myObsFactory a reference to the Observables Factory  
     * @param[in] myModel a pointer to the model
     * @param[in] rank the rank of the process that is using the parser
     * @return the iterator for the line being parsed
     */
    boost::tokenizer<boost::char_separator<char> >::iterator & ParseHiggsObservable(boost::tokenizer<boost::char_separator<char> >::iterator & beg, 
                                                                                    ThObsFactory& myObsFactory,
                                                                                    StandardModel * myModel,
                                                                                    int rank);

    private:
        TMatrixD channels;///< A matrix holding the information of all the channels.
        bool isnew;///< A boolean which is true if the parametrization is new.
        std::vector<ThObservable*> thObsV;///< A vector of ThObservables.
        std::vector<double> theoryValues;///< A vector to contain the theoryvalues.
        double thobsvsize;///< The size of the theory observables vector.
        double channelsize;///< The number of channels in the Higgs Observable.
        bool isCorrelated;///< A flag for correlated Higgs Observable.
        bool covarianceFromConfig;///< A flag for reading inverse covariance from config file.
        gslpp::matrix<double>* InvCov;///< The inverse covariance matrix.
        int rank;///< The rank of the process initializing this observable.
};

#endif	/* HIGGSOBSERVABLE_H */

