/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CORRELATEDGAUSSIANOBSERVABLES_H
#define	CORRELATEDGAUSSIANOBSERVABLES_H

#include "Observable.h"

/**
 * @addtogroup Observable
 * @brief A project for observable. 
 * @{
 */

/**
 * @class CorrelatedGaussianObservables
 * @ingroup Observable
 * @brief A class for correlated Gaussian observables. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class CorrelatedGaussianObservables {
public:
    CorrelatedGaussianObservables(std::string name_i);
    CorrelatedGaussianObservables(const CorrelatedGaussianObservables& orig);
    virtual ~CorrelatedGaussianObservables();

    void ComputeCov(gslpp::matrix<double> Corr);

    void AddObs(Observable& Obs_i);

    std::vector<Observable> getObs() const 
    {
        return Obs;
    }

    std::string getName() const 
    {
        return name;
    }

    void setName(std::string name)
    {
        this->name = name;
    }

    gslpp::matrix<double> getCov() const
    {
        return *Cov;
    }
    
private:
    std::vector<Observable> Obs;
    gslpp::matrix<double>* Cov;
    std::string name;
};

/** 
 * @}
 */

#endif	/* CORRELATEDGAUSSIANOBSERVABLES_H */

