/* 
 * File:   BParameter.h
 * Author: silvest
 *
 * Created on June 8, 2011, 5:35 PM
 */

#ifndef BPARAMETER_H
#define	BPARAMETER_H

#include <gslpp_vector_double.h>
#include "OrderScheme.h"

using namespace gslpp;

class BParameter {
public:
    BParameter(int n) : bpars(n,0.) {};

    vector<double> getBpars() const {
        return bpars;
    }

    void setBpars(vector<double> bpars) {
        this->bpars = bpars;
    }

    void setBpars(int i, double value) {
        this->bpars(i) = value;
    }

    double getMu() const {
        return mu;
    }

    void setMu(double mu) {
        this->mu = mu;
    }

    schemes getScheme() const {
        return scheme;
    }

    void setScheme(schemes scheme) {
        this->scheme = scheme;
    }

private:
    vector<double> bpars;
    double mu;
    schemes scheme;
};

#endif	/* BPARAMETER_H */

