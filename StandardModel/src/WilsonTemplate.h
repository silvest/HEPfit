/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
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
                   orders_qed order_qed_i = NO_QED)
    {
        size = dim;
        scheme = scheme_i;
        order = order_i;
        order_qed = order_qed_i; 
        mu = -1.;
        
        for (int i = LO; i <= MAXORDER; i++)
            if (i <= order)
                elem[i] = new T(size, 0.);
            else
                elem[i] = NULL;
        elem[orders_qed(NO_QED)] = NULL;
        //for (int i = LO_QED; i <= NLO_QED; i++){
        for (int i = LO_QED; i <= MAXORDER_QED; i++){
            if (i <= order_qed)
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
        order_qed = orig.order_qed;
        mu = orig.mu;
        for (int i = LO; i <= MAXORDER_QED; i++)
            if (orig.elem[i]!= NULL)
                elem[i] = new T(*(orig.elem[i]));
            else
                elem[i] = NULL;
    }
    
    virtual ~WilsonTemplate()
    {
        for (int i = LO; i <= MAXORDER_QED; i++)
            if (elem[i] != NULL)
                delete elem[i];
    };

    orders getOrder() const 
    {
        return order;
    }
    
    orders_qed getOrder_qed() const
    {
        return order_qed;
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
        if (order_qed != NO_QED){
            for(int i = LO_QED; i <= order_qed; i++){
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
    T* elem[MAXORDER_QED+1];
    unsigned int size;
    double mu;
    schemes scheme;
    orders order;
    orders_qed order_qed;

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
    
    T * Elem(orders_qed order_qed) const
    {
        if ((order_qed > this->order_qed)) {
            std::stringstream out;
            out << order_qed;
            throw std::runtime_error("WilsonTemplate::getElem(): requested order_qed " + out.str() +
                    "not present in the object");
        }
        return elem[order_qed];
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
    
    void setElem(const T & v, orders_qed order_qed_i)
    {
        if (order_qed_i > order_qed) {
            std::stringstream out;
            out << order_qed_i;
            throw std::runtime_error("MatchingCondition::setElem(): order " + out.str() +
                    " not implemented ");
        }
        *elem[order_qed_i] = v;
    };
};

#endif	/* WILSONTEMPLATE_H */

