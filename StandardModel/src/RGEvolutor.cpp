/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RGEvolutor.h"

RGEvolutor::RGEvolutor(unsigned int dim, schemes scheme, orders order)
: WilsonTemplate<gslpp::matrix<double> >(dim, scheme, order)
{}

RGEvolutor::RGEvolutor(unsigned int dim, schemes scheme, orders order, orders_ew order_ew)
: WilsonTemplate<gslpp::matrix<double> >(dim, scheme, order, order_ew)
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

void RGEvolutor::setEvol(unsigned int i, unsigned int  j, double x, orders order_i, orders_ew order_ew_i) 
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
    
    if (order_ew != NULL_ew){
    if (i > size || j > size) {
        std::stringstream out;
        out << i << " " << j;
        throw std::runtime_error("RGEvolutor::setEvol(): matrix indices " + out.str() + " out of range"); 
    }
    if (order_ew_i > order_ew) {
        std::stringstream out;
        out << order_i;
        throw std::runtime_error("RGEvolutor::setEvol(): order " + out.str() +" not implemented "); 
    }
    (*elem[order_ew_i])(i,j) = x;
    }
}

void RGEvolutor::setEvol(const gslpp::matrix<double>& m, orders order_i)
{
    setElem(m, order_i);
}

void RGEvolutor::setEvol(const gslpp::matrix<double>& m, orders_ew order_ew_i)
{
    setElem(m, order_ew_i);
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
    *(elem[LO]) = gslpp::matrix<double>::Id(size);
    for(int i = NLO; i <= order; i++)
        *(elem[i]) = 0.;
    
    if (order_ew != NULL_ew){
        for(int i = NLO_ew; i <= order_ew; i++)
            *(elem[i]) = 0.;
    }
}

void RGEvolutor::setM(double M)
{
    this->M = M;
    *(elem[LO]) = gslpp::matrix<double>::Id(size);
    for(int i = NLO; i <= order; i++)
        *(elem[i]) = 0.;
    
    if (order_ew != NULL_ew){
        for(int i = NLO_ew; i <= order_ew; i++)
            *(elem[i]) = 0.;
    }
}

void RGEvolutor::setMu(double mu)
{
    this->mu = mu;
    *(elem[LO]) = gslpp::matrix<double>::Id(size);
    for(int i = NLO; i <= order; i++)
        *(elem[i]) = 0.;
    
    if (order_ew != NULL_ew){
        for(int i = NLO_ew; i <= order_ew; i++)
            *(elem[i]) = 0.;
    }
}

gslpp::matrix<double>* RGEvolutor::Evol(orders order)
{
    return Elem(order);
}

gslpp::matrix<double>* RGEvolutor::Evol(orders_ew order_ew)
{
    return Elem(order_ew);
}
