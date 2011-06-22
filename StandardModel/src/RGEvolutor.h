/* 
 * File:   RGEvolutor.h
 * Author: marco
 *
 * Created on May 11, 2011, 5:16 PM
 */

#ifndef RGEVOLUTOR_H
#define	RGEVOLUTOR_H

#include <gslpp_matrix_double.h>
#include "OrderScheme.h"
#include "WilsonTemplate.h"

using namespace gslpp;

class RGEvolutor : public WilsonTemplate<matrix<double> > {
public:
    RGEvolutor(unsigned int dim, schemes scheme, orders order) : 
    WilsonTemplate<matrix<double> >(dim, scheme, order) {};
    
    matrix<double>** getEvol() const {
        return (matrix<double>**) elem;
    }

    double getM() const {
        return M;
    }

    void setScales(double mu, double M) {
        this->M = M;
        this->mu = mu;
        for(int i = LO; i <= order; i++)
          *(elem[i]) = matrix<double>::Id(size);
    }

    void setM(double M) {
        this->M = M;
        for(int i = LO; i <= order; i++)
          *(elem[i]) = matrix<double>::Id(size);
    }
    
    void setMu(double mu) {
        this->mu = mu;
        for(int i = LO; i <= order; i++)
          *(elem[i]) = matrix<double>::Id(size);
    }
    
    void setEvol(unsigned int i, unsigned int j, double x, orders order_i);

    void setEvol(const matrix<double>& m, orders order_i) { 
        setElem(m, order_i); }
    
    matrix<double>* Evol(orders order) { return Elem(order); };

protected:
    double M;
};

#endif	/* RGEVOLUTOR_H */
