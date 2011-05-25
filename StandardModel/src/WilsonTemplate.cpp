/* 
 * File:   WilsonTemplate.cpp
 * Author: enrico
 * 
 * Created on May 12, 2011, 10:52 AM
 */

#include "WilsonTemplate.h"
#include <sstream>

template <typename T> WilsonTemplate<T>::WilsonTemplate(unsigned int dim, schemes scheme_i, orders order_i) {
    size = dim;
    scheme = scheme_i;
    order = order_i;
    mu = -1.;
    for (int i = LO; i <= MAXORDER; i++)
        if (i <= order)
            elem[i] = new T(size, 0.);
        else
            elem[i] = NULL;
}

template <typename T> WilsonTemplate<T>::WilsonTemplate(const WilsonTemplate<T>& orig) {
    size = orig.getSize();
    scheme = orig.getScheme();
    order = orig.getOrder();
    mu = orig.getMu();

    for (int i = LO; i <= MAXORDER; i++) {
        orders ord = orders(i);
        if (orig.Elem(ord) != NULL)
            elem[ord] = new T(*(orig.Elem(ord)));
        else
            elem[ord] = NULL;
    }
}

template <typename T> WilsonTemplate<T>::~WilsonTemplate() {
    for (int i = LO; i <= MAXORDER; i++)
        if (elem[i] != NULL)
            delete elem[i];
}

template <typename T> void WilsonTemplate<T>::setElem(const T& v, orders order_i) {
    if (order_i > order) {
        std::stringstream out;
        out << order_i;
        throw "MatchingCondition::setElem(): order " + out.str() +
                " not implemented ";
    }
    for (unsigned int i = 0; i < size; i++)
        elem[order_i]->assign(i, v(i));
}

template <typename T> T* WilsonTemplate<T>::Elem(orders order) const {
    if (order > this->order) {
        std::stringstream out;
        out << order;
        throw "WilsonTemplate::getElem(): requested order " + out.str() +
                "not present in the object";
    }
    return elem[order];
}
