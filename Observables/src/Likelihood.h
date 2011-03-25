/* 
 * File:   Likelihood.h
 * Author: silvest
 *
 * Created on February 22, 2011, 12:20 PM
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

