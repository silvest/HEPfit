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
 * @details  
 */
class GaussianLikelihood : public Likelihood {
public:
    GaussianLikelihood(const double mean, const double sigma);
    GaussianLikelihood(const GaussianLikelihood& orig);
    virtual ~GaussianLikelihood();
    double getLikelihood(const double value) const;
private:
    double mean, sigma;
};

#endif	/* GAUSSIANLIKELIHOOD_H */
