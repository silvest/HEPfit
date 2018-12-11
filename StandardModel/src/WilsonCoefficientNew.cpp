/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "WilsonCoefficientNew.h"
#include <sstream>
#include <stdexcept>

WilsonCoefficientNew::WilsonCoefficientNew(unsigned int dim, schemes scheme, qcd_orders order_qcd_i, qed_orders order_qed_i)
: WilsonTemplateNew<gslpp::vector<gslpp::complex> >(dim, scheme, order_qcd_i, order_qed_i) {
};

Expanded<gslpp::complex> WilsonCoefficientNew::getCoeffElement(int i) const {
    Expanded<gslpp::complex> ret;

    if (i >= size) {
        std::stringstream out;
        out << i;
        throw std::runtime_error("WilsonTemplate::getCoeff(): requested element " + out.str() +
                " not present in the object");
    }
    std::vector<std::vector<gslpp::complex> > obj(wilson.getN1());
    for (int j = 0; j < wilson.getN1(); j++)
        for (int k = 0; k < wilson.getN2().at(j); k++)
            obj[j].push_back(wilson.getOrd(j, k)(i));
    return (Expanded<gslpp::complex>(obj));
};

void WilsonCoefficientNew::setCoeff(int i, gslpp::complex z, qcd_orders order_qcd_i, qed_orders order_qed_i) {
    if (i >= size) {
        std::stringstream out;
        out << i;
        throw std::runtime_error("WilsonTemplate::setCoeff(): coefficient index "
                + out.str() + " out of range");
    }
    if (order_qcd_i > order_qcd || order_qed_i > order_qed) {
        std::stringstream out;
        out << order_qcd_i << " and " << order_qed_i;
        throw std::runtime_error("WilsonTemplate::setCoeff(): order " + out.str() +
                " not implemented ");
    }
    gslpp::vector<gslpp::complex> tmp = wilson.getOrd(order_qcd_i, order_qed_i);
    tmp.assign(i, z);
    wilson.setOrd(order_qcd_i, order_qed_i, tmp);
}

void WilsonCoefficientNew::setCoeff(const gslpp::vector<gslpp::complex>& v, qcd_orders order_qcd_i, qed_orders order_qed_i) {
    setWilson(v, order_qcd_i, order_qed_i);
}

void WilsonCoefficientNew::resetCoeff() {
    resetWilson();
}

gslpp::vector<gslpp::complex> WilsonCoefficientNew::getCoeff(qcd_orders order_qcd_i, qed_orders order_qed_i) const {
    return getWilson(order_qcd_i, order_qed_i);
}

Expanded<gslpp::vector<gslpp::complex> > WilsonCoefficientNew::getCoeff() const {
    return getWilson();
}

void WilsonCoefficientNew::setCoeff(const Expanded<gslpp::vector<gslpp::complex> > wc) {
    wilson = wc;
}
