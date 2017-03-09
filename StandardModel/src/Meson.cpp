/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Meson.h"
#include "std_make_vector.h"

Meson::Meson(double mass, double lifetime = 5.e29, double decayconst = 0., 
        double lambdaM = 0., double gegenalpha1 = 0., double gegenalpha2 = 0.)
{
    this->mass = mass;
    this->lifetime = lifetime;
    this->decayconst = decayconst;
    this->lambdaM = lambdaM;
    gegenalpha[0] = gegenalpha1;
    gegenalpha[1] = gegenalpha2;
}

Meson::~Meson()
{}

void initializeParameters()
{}

double Meson::computeWidth() const
{
    return (HCUT / lifetime);
}

std::vector<std::string> Meson::injectParameterList(std::string mesonName_i)
{
//    lifetimeName = ("t" + mesonName_i).c_str(); 
//    gegenalphaName[0] = ("alpha1" + mesonName_i).c_str();
//    gegenalphaName[1] = ("alpha2" + mesonName_i).c_str();
//    lambdaMName; 
//    Dgamma_gammaName;
    return make_vector<std::string>() << "XXX" << "YYY" << "ZZZ";
}