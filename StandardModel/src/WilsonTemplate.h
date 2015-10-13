/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef WILSONTEMPLATE_H
#define	WILSONTEMPLATE_H

#include "OrderScheme.h"
#include <sstream>
#include <stdexcept>

/**
 * @class WilsonTemplate
 * @ingroup StandardModel
 * @brief A template class for the Wilson coefficients. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
template <class T> class WilsonTemplate {
public:
    WilsonTemplate(unsigned int dim, schemes scheme_i, orders order_i , 
                   orders_ew order_ew_i = NULL_ew)
    {
        size = dim;
        scheme = scheme_i;
        order = order_i;
        order_ew = order_ew_i; 
        mu = -1.;
        
        for (int i = LO; i <= MAXORDER; i++)
            if (i <= order)
                elem[i] = new T(size, 0.);
            else
                elem[i] = NULL;
        elem[orders_ew(NULL_ew)] = NULL;
        //for (int i = LO_ew; i <= NLO_ew; i++){
        for (int i = LO_ew; i <= MAXORDER_EW; i++){
            if (i <= order_ew)
                elem[i] = new T(size, 0.);
            else
                elem[i] = NULL;
        }
    };
    
    WilsonTemplate<T>(const WilsonTemplate<T>& orig) 
    {
        size = orig.size;
        scheme = orig.scheme;
        order = orig.order;
        order_ew = orig.order_ew;
        mu = orig.mu;
        for (int i = LO; i <= MAXORDER_EW; i++)
            if (orig.elem[i]!= NULL)
                elem[i] = new T(*(orig.elem[i]));
            else
                elem[i] = NULL;
    }
    
    virtual ~WilsonTemplate()
    {
        for (int i = LO; i <= MAXORDER_EW; i++)
            if (elem[i] != NULL)
                delete elem[i];
    };

    orders getOrder() const 
    {
        return order;
    }
    
    orders_ew getOrder_ew() const
    {
        return order_ew;
    }

    double getMu() const
    {
        return mu;
    }

        virtual void resetCoefficient()
    {
        for(int i = LO; i <= order; i++){
            *(elem[i]) = 0.;
        }
        if (order_ew != NULL_ew){
            for(int i = LO_ew; i <= order_ew; i++){
                *(elem[i]) = 0.;
            }
        }
    }

    virtual void setMu(double mu)
    {
        this->mu = mu;
        resetCoefficient();
    }
    
    schemes getScheme() const 
    {
        return scheme;
    }

    void setScheme(schemes scheme)
    {
        this->scheme = scheme;
    }

    unsigned int  getSize() const
    {
        return size;
    }

protected:
    T* elem[MAXORDER_EW+1];
    unsigned int size;
    double mu;
    schemes scheme;
    orders order;
    orders_ew order_ew;

    T * Elem(orders order) const
    {
        if (order > this->order) {
            std::stringstream out;
            out << order;
            throw std::runtime_error("WilsonTemplate::getElem(): requested order " + out.str() +
                    " not present in the object");
        }
        return elem[order];
    };
    
    T * Elem(orders_ew order_ew) const
    {
        if ((order_ew > this->order_ew)) {
            std::stringstream out;
            out << order_ew;
            throw std::runtime_error("WilsonTemplate::getElem(): requested order_ew " + out.str() +
                    "not present in the object");
        }
        return elem[order_ew];
    };

    void setElem(const T & v, orders order_i)
    {
        if (order_i > order) {
            std::stringstream out;
            out << order_i;
            throw std::runtime_error("MatchingCondition::setElem(): order " + out.str() +
                    " not implemented ");
        }
        *elem[order_i] = v;
    };
    
    void setElem(const T & v, orders_ew order_ew_i)
    {
        if (order_ew_i > order_ew) {
            std::stringstream out;
            out << order_ew_i;
            throw std::runtime_error("MatchingCondition::setElem(): order " + out.str() +
                    " not implemented ");
        }
        *elem[order_ew_i] = v;
    };
};

#endif	/* WILSONTEMPLATE_H */

