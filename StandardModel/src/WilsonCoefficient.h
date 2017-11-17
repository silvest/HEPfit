/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
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
    
    WilsonCoefficient(unsigned int dim, schemes scheme, orders order, orders_qed order_qed) 
    : WilsonTemplate<gslpp::vector<gslpp::complex> >(dim, scheme, order, order_qed) 
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
    
    void setCoeff(const gslpp::vector<gslpp::complex>& z, orders_qed order_qed_i) 
    { 
        setElem(z, order_qed_i); 
    };

    void setCoeff(unsigned int i, gslpp::complex z, orders order_i);
    
    void setCoeff(unsigned int i, gslpp::complex z, orders_qed order_qed_i);

    gslpp::vector<gslpp::complex>* getCoeff(orders ord) const { return Elem(ord); };
    
    gslpp::vector<gslpp::complex>* getCoeff(orders_qed ord_qed) const { return Elem(ord_qed); };
    
};

#endif	/* WILSONCOEFFICIENT_H */
