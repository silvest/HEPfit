/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
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
    
    RGEvolutor(unsigned int dim, schemes scheme, orders order, orders_ew order_ew) : 
    WilsonTemplate<matrix<double> >(dim, scheme, order, order_ew) {};
    
    matrix<double>** getEvol() const {
        return (matrix<double>**) elem;
    }

    double getM() const {
        return M;
    }

    void setScales(double mu, double M) {
        this->M = M;
        this->mu = mu;
        *(elem[LO]) = matrix<double>::Id(size);
        for(int i = NLO; i <= order; i++)
          *(elem[i]) = 0.;
        
        if (order_ew != NULL_ew){
            for(int i = NLO_ew; i <= order_ew; i++)
                *(elem[i]) = 0.;
        }
    }

    void setM(double M) {
        this->M = M;
        *(elem[LO]) = matrix<double>::Id(size);
        for(int i = NLO; i <= order; i++)
          *(elem[i]) = 0.;
        
        if (order_ew != NULL_ew){
            for(int i = NLO_ew; i <= order_ew; i++)
                *(elem[i]) = 0.;
        }
    }
    
    void setMu(double mu) {
        this->mu = mu;
        *(elem[LO]) = matrix<double>::Id(size);
        for(int i = NLO; i <= order; i++)
          *(elem[i]) = 0.;
        
        if (order_ew != NULL_ew){
            for(int i = NLO_ew; i <= order_ew; i++)
                *(elem[i]) = 0.;
        }
    }
    
    void setEvol(unsigned int i, unsigned int j, double x, orders order_i);
    void setEvol(unsigned int i, unsigned int j, double x, orders order_i, orders_ew order_ew) ;

    void setEvol(const matrix<double>& m, orders order_i) { 
        setElem(m, order_i); }
    
    void setEvol(const matrix<double>& m, orders_ew order_ew_i) { 
        setElem(m, order_ew_i); }
    
    matrix<double>* Evol(orders order) { 
        return Elem(order);}
        
    matrix<double>* Evol(orders_ew order_ew) { 
        return Elem(order_ew);}
    
protected:
    double M;
};

#endif	/* RGEVOLUTOR_H */

