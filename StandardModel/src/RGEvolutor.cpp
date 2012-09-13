/* 
 * File:   RGEvolutor.cpp
 * Author: marco
 * 
 * Created on May 11, 2011, 5:16 PM
 */

#include "RGEvolutor.h"
#include <sstream>
#include <stdexcept>

void RGEvolutor::setEvol(unsigned int i, unsigned int  j, double x, orders order_i) {
    
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

void RGEvolutor::setEvol(unsigned int i, unsigned int  j, double x, orders order_i, orders_ew order_ew_i) {
    
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


