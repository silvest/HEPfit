/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef WILSONCOEFFICIENT_H
#define	WILSONCOEFFICIENT_H

#include <gslpp_vector_complex.h>
#include "WilsonTemplate.h"

/**
 * @class WilsonCoefficient
 * @ingroup StandardModel
 * @brief A class for the Wilson coefficients. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class WilsonCoefficient : public WilsonTemplate<gslpp::vector<gslpp::complex> > {
public:

    WilsonCoefficient(unsigned int dim, schemes scheme, orders order) 
    : WilsonTemplate<gslpp::vector<gslpp::complex> >(dim, scheme, order) 
    {
    };
    
    WilsonCoefficient(unsigned int dim, schemes scheme, orders order, orders_ew order_ew) 
    : WilsonTemplate<gslpp::vector<gslpp::complex> >(dim, scheme, order, order_ew) 
    {
    };

    gslpp::vector<gslpp::complex>** getCoeff() const
    {
        return (gslpp::vector<gslpp::complex>**) elem;
    }

    void setCoeff(const gslpp::vector<gslpp::complex>& z, orders order_i)
    { 
        setElem(z, order_i); 
    };
    
    void setCoeff(const gslpp::vector<gslpp::complex>& z, orders_ew order_ew_i) 
    { 
        setElem(z, order_ew_i); 
    };

    void setCoeff(unsigned int i, gslpp::complex z, orders order_i);
    
    void setCoeff(unsigned int i, gslpp::complex z, orders_ew order_ew_i);

    gslpp::vector<gslpp::complex>* getCoeff(orders ord) const { return Elem(ord); };
    
    gslpp::vector<gslpp::complex>* getCoeff(orders_ew ord_ew) const { return Elem(ord_ew); };

};

#endif	/* WILSONCOEFFICIENT_H */

