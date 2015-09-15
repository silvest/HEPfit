/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "CorrelatedGaussianObservables.h"

CorrelatedGaussianObservables::CorrelatedGaussianObservables(std::string name_i)
{
    name = name_i;
}

CorrelatedGaussianObservables::CorrelatedGaussianObservables(const CorrelatedGaussianObservables& orig)
{
    Obs = orig.Obs;
    name = orig.name;
    Cov = new gslpp::matrix<double>(*orig.Cov);
}

CorrelatedGaussianObservables::~CorrelatedGaussianObservables()
{   
    Cov = new gslpp::matrix<double>(10, 10, 0.); /** Put in to prevent seg fault during error handling in InputParser **/
    delete(Cov);
}

void CorrelatedGaussianObservables::AddObs(Observable& Obs_i)
{
    Obs.push_back(Obs_i);
}

void CorrelatedGaussianObservables::ComputeCov(gslpp::matrix<double> Corr)
{
    unsigned int size = Obs.size();
    if (Corr.size_i() != size || Corr.size_j() != size)
        throw std::runtime_error("The size of the correlated observables in " + name + " does not match the size of the correlation matrix!");
    Cov = new gslpp::matrix<double>(size, size, 0.);
    for (unsigned int i = 0; i < size; i++)
        for (unsigned int j = 0; j < size; j++)
            (*Cov)(i, j) = Obs.at(i).getErrg() * Corr(i, j) * Obs.at(j).getErrg();
    *Cov = Cov->inverse();
    
    gslpp::matrix<gslpp::complex>vv(size, size, 0.);
    gslpp::vector<gslpp::complex>ee(size, 0.);

    Cov->eigensystem(vv, ee);

    gslpp::matrix<double> v(vv.real());
    gslpp::vector<double> e(ee.real());

    std::cout << v << std::endl;
    std::cout << e << std::endl;
    for (int i=0; i< size; i++)
        std::cout << 1./sqrt(e(i)) << std::endl;
}

double CorrelatedGaussianObservables::computeWeight()
{
    gslpp::vector<double> x(Obs.size());

    for (unsigned int i = 0; i < Obs.size(); i++)
        x(i) = Obs.at(i).computeTheoryValue() - Obs.at(i).getAve();

    return (-0.5 * x * ((*Cov) * x));
}
