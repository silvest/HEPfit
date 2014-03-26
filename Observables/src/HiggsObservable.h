/* 
 * File:   HiggsObservable.h
 * Author: silvest
 *
 * Created on March 18, 2014, 4:58 PM
 */

#ifndef HIGGSOBSERVABLE_H
#define	HIGGSOBSERVABLE_H

#include "Observable.h"
#include <TMatrixD.h>

/**
 * @class HiggsObservable
 * @ingroup Observable
 * @brief A class for Higgs experimental analyses
 * @author SusyFit Collaboration
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
               };
               
    HiggsObservable(const HiggsObservable& orig);

    virtual ~HiggsObservable();
        // Read the necessary information from the config file. Each row contains:
        // ggH fraction
        // VBF fraction
        // VH fraction (ttH is computed as 1-ggH-VBF-VH)

    /**
     * Set the parametric likelihood describing one Higgs decay channel from a config file.
     * @param filename the name of the config file
     * @param thObsV a vector of ThObservables containing the ratio of the production cross section for 
     * ggH, VBF, VH and ttH in the model analyzed over the SM prediction 
     */
    virtual void setParametricLikelihood(std::string filename, std::vector<ThObservable*> thObsV);

    virtual double computeWeight();
    
    private:
        TMatrixD channels;
        std::vector<ThObservable*> thObsV;
};

#endif	/* HIGGSOBSERVABLE_H */

