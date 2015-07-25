/* 
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "DoxyExample.h"

DoxyExample::DoxyExample(std::string name_i)
{
    name = name_i;
}

DoxyExample::DoxyExample(const DoxyExample& orig)
{
    Obs = orig.Obs;
    name = orig.name;
    Cov = new gslpp::matrix<double>(*orig.Cov);
}

DoxyExample::~DoxyExample()
{
    if(Cov != NULL)
    delete(Cov);
}

void DoxyExample::AddObs(Observable& Obs_i)
{
    Obs.push_back(Obs_i);
}

void DoxyExample::ComputeCov(gslpp::matrix<double> Corr) 
{
    int size = Obs.size();
    if (Corr.size_i() != size || Corr.size_j() != size)
        throw std::runtime_error("The size of the correlated observables in "+name+" does not match the size of the correlation matrix!");
    Cov = new gslpp::matrix<double>(size,size,0.);
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            (*Cov)(i,j) = Obs.at(i).getErrg()*Corr(i,j)*Obs.at(j).getErrg();
    *Cov = Cov->inverse();
}
