/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RGEvolutor.h"

RGEvolutor::RGEvolutor(unsigned int dim, schemes scheme, orders order)
: WilsonTemplate<gslpp::matrix<double> >(dim, scheme, order)
{}

RGEvolutor::RGEvolutor(unsigned int dim, schemes scheme, orders order, orders_qed order_qed)
: WilsonTemplate<gslpp::matrix<double> >(dim, scheme, order, order_qed)
{}

RGEvolutor::~RGEvolutor()
{}

void RGEvolutor::setEvol(unsigned int i, unsigned int  j, double x, orders order_i) 
{    
    if (i > size || j > size) {
        std::stringstream out;
        out << i << " " << j;
        throw std::runtime_error("RGEvolutor::setEvol(): matrix indices " + out.str() + " out of range"); 
    }
    if (order_i > order) {
        std::stringstream out;
        out << order_i;
        throw std::runtime_error("RGEvolutor::setEvol(): order " + out.str() +" not implemented "); 
    }
    (*elem[order_i])(i,j) = x;   
}

void RGEvolutor::setEvol(unsigned int i, unsigned int  j, double x, orders order_i, orders_qed order_qed_i) 
{    
    if (i > size || j > size) {
        std::stringstream out;
        out << i << " " << j;
        throw std::runtime_error("RGEvolutor::setEvol(): matrix indices " + out.str() + " out of range"); 
    }
    if (order_i > order) {
        std::stringstream out;
        out << order_i;
        throw std::runtime_error("RGEvolutor::setEvol(): order " + out.str() +" not implemented "); 
    }
    (*elem[order_i])(i,j) = x;
    
    if (order_qed != NO_QED){
    if (i > size || j > size) {
        std::stringstream out;
        out << i << " " << j;
        throw std::runtime_error("RGEvolutor::setEvol(): matrix indices " + out.str() + " out of range"); 
    }
    if (order_qed_i > order_qed) {
        std::stringstream out;
        out << order_i;
        throw std::runtime_error("RGEvolutor::setEvol(): order " + out.str() +" not implemented "); 
    }
    (*elem[order_qed_i])(i,j) = x;
    }
}

void RGEvolutor::setEvol(const gslpp::matrix<double>& m, orders order_i)
{
    setElem(m, order_i);
}

void RGEvolutor::setEvol(const gslpp::matrix<double>& m, orders_qed order_qed_i)
{
    setElem(m, order_qed_i);
}

gslpp::matrix<double>** RGEvolutor::getEvol() const
{
    return (gslpp::matrix<double>**) elem;
}

double RGEvolutor::getM() const
{
    return M;
}

void RGEvolutor::setScales(double mu, double M)
{
    this->M = M;
    this->mu = mu;
    resetCoefficient();
    *(elem[LO]) = gslpp::matrix<double>::Id(size);
}

void RGEvolutor::setM(double M)
{
    setScales(mu, M);
}

void RGEvolutor::setMu(double mu)
{
    setScales(mu, M);
}

gslpp::matrix<double>* RGEvolutor::Evol(orders order)
{
    return Elem(order);
}

gslpp::matrix<double>* RGEvolutor::Evol(orders_qed order_qed)
{
    return Elem(order_qed);
}
