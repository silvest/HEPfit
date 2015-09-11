/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include <vector>
#include <sstream>
#include <math.h>
#include "CorrelatedGaussianParameters.h"

CorrelatedGaussianParameters::CorrelatedGaussianParameters(std::string name_i)
{
    name = name_i;
}

CorrelatedGaussianParameters::CorrelatedGaussianParameters(const CorrelatedGaussianParameters& orig)
{
    Pars = orig.Pars;
    name = orig.name;
    Cov = new gslpp::matrix<double>(*orig.Cov);
    v = new gslpp::matrix<double>(*orig.v);
    e = new gslpp::vector<double>(*orig.e);
}

CorrelatedGaussianParameters::~CorrelatedGaussianParameters()
{   
    //Cov = new gslpp::matrix<double>(10, 10, 0.); /** Put in to prevent seg fault during error handling in InputParser **/
    if(Cov != NULL)
        delete(Cov);
    if(v != NULL)
        delete(v);
    if(e != NULL)
        delete(e);
}

void CorrelatedGaussianParameters::AddPar(ModelParameter& Par_i)
{
    Pars.push_back(Par_i);
}

void CorrelatedGaussianParameters::DiagonalizePars(gslpp::matrix<double> Corr)
{
    unsigned int size = Pars.size();
    if (Corr.size_i() != size || Corr.size_j() != size)
        throw std::runtime_error("The size of the correlated parameters in " + name + " does not match the size of the correlation matrix!");
    Cov = new gslpp::matrix<double>(size, size, 0.);
    for (unsigned int i = 0; i < size; i++)
        for (unsigned int j = 0; j < size; j++)
            (*Cov)(i, j) = Pars.at(i).errg * Corr(i, j) * Pars.at(j).errg;
    *Cov = Cov->inverse();
    
    gslpp::matrix<gslpp::complex>vv(size, size, 0.);
    gslpp::vector<gslpp::complex>ee(size, 0.);

    Cov->eigensystem(vv, ee);

    *v = vv.real();
    *e = ee.real();
    gslpp::matrix<double> vi = v->inverse();
    
    std::cout << v << std::endl;
    std::cout << e << std::endl;
    
    gslpp::vector<double> ave_in(size, 0.);
    
    int ind = 0;
    for(std::vector<ModelParameter>::iterator it = Pars.begin(); it != Pars.end(); it++)
    {
        ave_in(ind) = it->ave;
        ind++;
    }
    
    gslpp::vector<double> ave = vi * ave_in;
    
    for (int i = 0; i < size; i++)
    {
        std::stringstream ss;
        ss << i;
        std::string namei = name + ss.str();
        DiagPars.push_back(ModelParameter(namei,ave(i),0.,1./sqrt((*e)(i))));
    }
}

