/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LIKELIHOOD_H
#define	LIKELIHOOD_H

class Likelihood {
public:
    Likelihood();
    Likelihood(const Likelihood& orig);
    virtual ~Likelihood();
    double getLikelihood(const double) const;
    double getLikelihood(const double, const double) const;
private:

};

#endif	/* LIKELIHOOD_H */

