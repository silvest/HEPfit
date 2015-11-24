/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HIGGSOBSERVABLE_H
#define	HIGGSOBSERVABLE_H

#include "Observable.h"
#include <TMatrixD.h>
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
    HiggsObservable(const Observable& Obs): 
               Observable(Obs)
               {
                   isnew = true;
               };
               
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
};

#endif	/* HIGGSOBSERVABLE_H */

