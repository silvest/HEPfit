/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef WILSONCOEFFICIENT_H
#define	WILSONCOEFFICIENT_H

#include <gslpp_vector_complex.h>
#include "WilsonTemplate.h"

using namespace gslpp;

/**
 * @class WilsonCoefficient
 * @ingroup StandardModel
 * @brief A class for the Wilson coefficients. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class WilsonCoefficient : public WilsonTemplate<vector<complex> > {
public:

    WilsonCoefficient(unsigned int dim, schemes scheme, orders order) 
    : WilsonTemplate<vector<complex> >(dim, scheme, order) 
    {
    };
    
    WilsonCoefficient(unsigned int dim, schemes scheme, orders order, orders_ew order_ew) 
    : WilsonTemplate<vector<complex> >(dim, scheme, order, order_ew) 
    {
    };

    vector<complex>** getCoeff() const
    {
        return (vector<complex>**) elem;
    }

    void setCoeff(const vector<complex>& z, orders order_i)
    { 
        setElem(z, order_i); 
    };
    
    void setCoeff(const vector<complex>& z, orders_ew order_ew_i) 
    { 
        setElem(z, order_ew_i); 
    };

    void setCoeff(unsigned int i, complex z, orders order_i);
    
    void setCoeff(unsigned int i, complex z, orders_ew order_ew_i);

    vector<complex>* getCoeff(orders ord) const { return Elem(ord); };
    
    vector<complex>* getCoeff(orders_ew ord_ew) const { return Elem(ord_ew); };

};

#endif	/* WILSONCOEFFICIENT_H */

