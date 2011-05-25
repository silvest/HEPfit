/* 
 * File:   WilsonCoefficient.h
 * Author: marco
 *
 * Created on May 11, 2011, 11:01 AM
 */

#ifndef WILSONCOEFFICIENT_H
#define	WILSONCOEFFICIENT_H

#include <gslpp_vector_complex.h>
#include "WilsonTemplate.h"

using namespace gslpp;

class WilsonCoefficient : public WilsonTemplate<vector<complex> > {
public:

    WilsonCoefficient(unsigned int dim, schemes scheme, orders order) :
    WilsonTemplate(dim, scheme, order) {
    };

    vector<complex>** getCoeff() const {
        return (vector<complex>**) elem;
    }

    void setCoeff(const vector<complex>& z, orders order_i) { 
        setElem(z, order_i); 
    };

    void setCoeff(unsigned int i, complex z, orders order_i);

    vector<complex>* Coeff(orders ord) const { return Elem(ord); };

};

#endif	/* WILSONCOEFFICIENT_H */
