/* 
 * File:   WilsonTemplate.h
 * Author: enrico
 *
 * Created on May 12, 2011, 10:52 AM
 */

#ifndef WILSONTEMPLATE_H
#define	WILSONTEMPLATE_H

#include "OrderScheme.h"

template <typename T> class WilsonTemplate {
public:
    WilsonTemplate(unsigned int dim, schemes scheme, orders order);
    WilsonTemplate(const WilsonTemplate<T> & orig);
    virtual ~WilsonTemplate();

    orders getOrder() const {
        return order;
    }

    double getMu() const {
        return mu;
    }

    void setMu(double mu) {
        this->mu = mu;
        for(int i = LO; i <= order; i++)
            *(elem[i]) = 0.;
    }
    
    schemes getScheme() const {
        return scheme;
    }

    void setScheme(schemes scheme) {
        this->scheme = scheme;
    }

    unsigned int  getSize() const {
        return size;
    }

protected:
    T* elem[MAXORDER+1];
    unsigned int size;
    double mu;
    schemes scheme;
    orders order;

    T * Elem(orders order) const;

    void setElem(const T &v, orders order_i);
};

#endif	/* WILSONTEMPLATE_H */

