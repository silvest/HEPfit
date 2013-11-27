/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAUSSIANLIKELIHOOD_H
#define	GAUSSIANLIKELIHOOD_H

#include "Likelihood.h"

/**
 * @class GaussianLikelihood
 * @ingroup Observable
 * @brief A class for likelihood function (Not used). 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for calculating the Gaussian likelihood given a 
 * certain mean and standard deviation
 */
class GaussianLikelihood : public Likelihood {
public:
    
    /**
     * @brief The default constructor.
     * @param[in] mean the mean of the Gaussian likelihood function
     * @param[in] sigma the mean of the Gaussian likelihood function
     */
    GaussianLikelihood(const double mean, const double sigma);
    
    /**
     * @brief The copy constructor.
     */
    GaussianLikelihood(const GaussianLikelihood& orig);
    
    /**
     * @brief The default destructor.
     */
    virtual ~GaussianLikelihood();
    
    /**
     * @brief A method to compute the Gaussian likelihood.
     * @param value the value at which the Gaussian likelihood will be computed
     */
    double computeLikelihood(const double value) const;
private:
    double mean; /**< The mean of the Gaussian likelihood function */
    double sigma; /**< The standard deviation of the Gaussian likelihood function */
};

#endif	/* GAUSSIANLIKELIHOOD_H */
