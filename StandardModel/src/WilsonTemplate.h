/* 
 * File:   WilsonTemplate.h
 * Author: enrico
 *
 * Created on May 12, 2011, 10:52 AM
 */

#ifndef WILSONTEMPLATE_H
#define	WILSONTEMPLATE_H

#include "OrderScheme.h"
#include <sstream>

template <class T> class WilsonTemplate {
public:
    WilsonTemplate(unsigned int dim, schemes scheme_i, orders order_i){
    size = dim;
    scheme = scheme_i;
    order = order_i;
    mu = -1.;
    for (int i = LO; i <= MAXORDER; i++)
        if (i <= order)
            elem[i] = new T(size, 0.);
        else
            elem[i] = NULL;
    };

    WilsonTemplate<T>(const WilsonTemplate<T>& orig) {
        size = orig.size;
        scheme = orig.scheme;
        order = orig.order;
        mu = orig.mu;
        for (int i = LO; i <= MAXORDER; i++)
            if (orig.elem[i]!= NULL)
                elem[i] = new T(*(orig.elem[i]));
            else
                elem[i] = NULL;
    }
    
    virtual ~WilsonTemplate(){
    for (int i = LO; i <= MAXORDER; i++)
        if (elem[i] != NULL)
            delete elem[i];
    };

    orders getOrder() const {
        return order;
    }

    double getMu() const {
        return mu;
    }

    virtual void setMu(double mu) {
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
    T* elem[MAXORDER + 1];
    unsigned int size;
    double mu;
    schemes scheme;
    orders order;

    T * Elem(orders order) const {
        if (order > this->order) {
            std::stringstream out;
            out << order;
            throw "WilsonTemplate::getElem(): requested order " + out.str() +
                    "not present in the object";
        }
        return elem[order];
    };

    void setElem(const T & v, orders order_i) {
        if (order_i > order) {
            std::stringstream out;
            out << order_i;
            throw "MatchingCondition::setElem(): order " + out.str() +
                    " not implemented ";
        }
        *elem[order_i] = v;
    };
};

#endif	/* WILSONTEMPLATE_H */

