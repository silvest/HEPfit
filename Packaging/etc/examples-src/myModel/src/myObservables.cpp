/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "myObservables.h"


myObservables::myObservables(const StandardModel& SM_i)
: ThObservable(SM_i), my_model(static_cast<const myModel*> (&SM_i))
{
    fact = 2. * 3000 * pow(0.95, 4.);
    kfact = 1.85;
}

myObservables::~myObservables()
{}

void myObservables::updateParameters()
{
    c1 = my_model->getc1();
    if (my_model->get_condition_flag() == true) {
        c2 = 1. - c1;
    } else {
        c2 = my_model->getc2();
    }
    c3 = my_model->getc3();
    c4 = my_model->getc4();
    sw2 = my_model->sW2();
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/

yield::yield(const StandardModel& SM_i, unsigned int bin_i)
: myObservables(SM_i)
{
    bin = bin_i;
}

double yield::computeThValue()
{
    updateParameters();
    
    double wt = kfact * sqrt(fact);
    
    if (bin == 1) return (13.0556 + 0.00433972*c4*c4 - 0.595884*c2 - 0.539262*c4*c4*c2 + 0.45963*c2*c2 - 0.844644*c1 - 0.896512*c4*c4*c1 + 1.40761*c2*c1 + 1.10521*c1*c1 - 0.0566324*c3*c3 + 0.00224709*c4*c4*c3*c3 + 0.000814862*c2*c3*c3 + 0.000885144*c1*c3*c3 + 0.252697*pow(c4,4.) + 0.000602038*pow(c3,4.)) * wt;
    
    if (bin == 2) return (2.59123 + 0.169824*c4*c4 - 0.113991*c2 - 0.637075*c4*c4*c2 + 0.52955*c2*c2 - 0.184628*c1 - 1.77388*c4*c4*c1 + 1.31494*c2*c1 + 1.45195*c1*c1 + 0.00228014*c3*c3 + 0.00281452*c4*c4*c3*c3 - 0.00121316*c2*c3*c3 + 0.00216552*c1*c3*c3 + 0.68041*pow(c4,4.) + 0.000218874*pow(c3,4.)) * wt;
    
    if (bin == 3) return (0.49923 + 0.0605644*c4*c4 - 0.0176069*c2 - 0.160722*c4*c4*c2 + 0.439888*c2*c2 - 0.0183051*c1 - 1.04927*c4*c4*c1 + 0.456769*c2*c1 + 0.695114*c1*c1 + 0.0115523*c3*c3 + 0.00023496*c4*c4*c3*c3 + 0.000200689*c2*c3*c3 + 0.00104028*c1*c3*c3 + 0.48041*pow(c4,4.) - 0.000532762*pow(c3,4.)) * wt;
    
    if (bin == 4) return (0.177548 + 0.0262649*c4*c4 - 0.00743058*c2 + 0.0546468*c4*c4*c2 + 0.524296*c2*c2 - 0.00407272*c1 - 0.748231*c4*c4*c1 + 0.14861*c2*c1 + 0.436721*c1*c1 + 0.00637587*c3*c3 + 0.000702175*c4*c4*c3*c3 + 0.0002901*c2*c3*c3 - 0.00040867*c1*c3*c3 + 0.375159*pow(c4,4.) + 5.37945e-6*pow(c3,4.)) * wt;
    
    if (bin == 5) return (0.0485093 + 0.00891872*c4*c4 - 0.00324976*c2 + 0.155788*c4*c4*c2 + 0.523622*c2*c2 + 0.00396845*c1 - 0.384512*c4*c4*c1 - 0.0424597*c2*c1 + 0.20254*c1*c1 + 0.00271808*c3*c3 + 0.000118455*c4*c4*c3*c3 + 0.000624283*c2*c3*c3 - 0.000523273*c1*c3*c3 + 0.202822*pow(c4,4.) - 5.43403e-6 *pow(c3,4.)) * wt;
    
    if (bin == 6) return (0.0165142 + 0.00240494*c4*c4 - 0.00151538*c2 + 0.443147*c4*c4*c2 + 2.01346*c2*c2 + 0.0066562*c1 - 0.285082*c4*c4*c1 - 0.329009*c2*c1 + 0.139285*c1*c1 + 0.00226039*c3*c3 + 0.000642254*c4*c4*c3*c3 - 0.000373886*c2*c3*c3 - 0.00122827*c1*c3*c3 + 0.154372*pow(c4,4.) - 0.000211814*pow(c3,4.)) * wt;
    else {
        std::cout << "Bin not defined" << std::endl;
        return (EXIT_FAILURE);
    }
}

C_3::C_3(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double C_3::computeThValue()
{
    updateParameters();
    
    return (1./2. - 4./3.*sw2)/(2. * (sqrt(sw2) * sqrt(1. - sw2))) * c3;
}

C_4::C_4(const StandardModel& SM_i)
: myObservables(SM_i)
{}

double C_4::computeThValue()
{
    updateParameters();
    
    return -1./2./(2. * (sqrt(sw2) * sqrt(1. - sw2))) * c4;
}
