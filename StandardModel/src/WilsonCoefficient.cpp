/* 
 * File:   WilsonCoefficient.cpp
 * Author: marco
 * 
 * Created on May 11, 2011, 11:01 AM
 */

#include "WilsonCoefficient.h"
#include <sstream>

void WilsonCoefficient::setCoeff(unsigned int i, complex z, orders order_i) {
    
    if ((unsigned int) i > size) {
        std::stringstream out;
        out << i;
        throw "WilsonCoefficient::setCoeff(): coefficient index "
        + out.str() + " out of range";
    }
    if (order_i > order) {
        std::stringstream out;
        out << order_i;
        throw "WilsonCoefficient::setCoeff(): order " + out.str() +
                " not implemented ";
    }
    elem[order_i]->assign(i, z);
}

void WilsonCoefficient::setCoeff(unsigned int i, complex z, orders_ew order_ew_i) {
    
    if ((unsigned int) i > size) {
        std::stringstream out;
        out << i;
        throw "WilsonCoefficientEW::setCoeff(): coefficient index "
        + out.str() + " out of range";
    }
    if (order_ew_i > order_ew) {
        std::stringstream out;
        out << order_ew_i;
        throw "WilsonCoefficientEW::setCoeff(): order_ew " + out.str() +
                " not implemented ";
    }
    elem[order_ew_i]->assign(i, z);
}

