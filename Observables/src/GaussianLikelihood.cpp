/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GaussianLikelihood.h"
#include <math.h>

GaussianLikelihood::GaussianLikelihood(const double mean, const double sigma) {
    this->mean = mean;
    this->sigma = sigma;
}

GaussianLikelihood::GaussianLikelihood(const GaussianLikelihood& orig) {
}

GaussianLikelihood::~GaussianLikelihood() {
}

double GaussianLikelihood::getLikelihood(const double value) const {
    return(exp(-(value-mean)*(value-mean)/(2.*sigma*sigma)));
}
