/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "gg4l.h"


gg4l::gg4l(const StandardModel& SM_i)
: ThObservable(SM_i), my_model(static_cast<const myModel*> (&SM_i))
{
    fact = 2. * 3000 * pow(0.95, 4.);
    kfact = 1.85;
}

gg4l::~gg4l()
{}

void gg4l::updateParameters()
{
    ct = my_model->getct();
    if (my_model->get_onshell_flag() == true) {
        cg = 1 - ct;
    } else {
        cg = my_model->getcg();
    }
    cV = my_model->getcV();
    cA = my_model->getcA();
    sw2 = my_model->sW2();
}

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/

yield::yield(const StandardModel& SM_i, unsigned int bin_i)
: gg4l(SM_i)
{
    bin = bin_i;
}

double yield::computeThValue()
{
    updateParameters();
    
    double wt = kfact * sqrt(fact);
    
    if (bin == 1) return (13.0556 + 0.00433972*cA*cA - 0.595884*cg - 0.539262*cA*cA*cg + 0.45963*cg*cg - 0.844644*ct - 0.896512*cA*cA*ct + 1.40761*cg*ct + 1.10521*ct*ct - 0.0566324*cV*cV + 0.00224709*cA*cA*cV*cV + 0.000814862*cg*cV*cV + 0.000885144*ct*cV*cV + 0.252697*pow(cA,4.) + 0.000602038*pow(cV,4.)) * wt;
    
    if (bin == 2) return (2.59123 + 0.169824*cA*cA - 0.113991*cg - 0.637075*cA*cA*cg + 0.52955*cg*cg - 0.184628*ct - 1.77388*cA*cA*ct + 1.31494*cg*ct + 1.45195*ct*ct + 0.00228014*cV*cV + 0.00281452*cA*cA*cV*cV - 0.00121316*cg*cV*cV + 0.00216552*ct*cV*cV + 0.68041*pow(cA,4.) + 0.000218874*pow(cV,4.)) * wt;
    
    if (bin == 3) return (0.49923 + 0.0605644*cA*cA - 0.0176069*cg - 0.160722*cA*cA*cg + 0.439888*cg*cg - 0.0183051*ct - 1.04927*cA*cA*ct + 0.456769*cg*ct + 0.695114*ct*ct + 0.0115523*cV*cV + 0.00023496*cA*cA*cV*cV + 0.000200689*cg*cV*cV + 0.00104028*ct*cV*cV + 0.48041*pow(cA,4.) - 0.000532762*pow(cV,4.)) * wt;
    
    if (bin == 4) return (0.177548 + 0.0262649*cA*cA - 0.00743058*cg + 0.0546468*cA*cA*cg + 0.524296*cg*cg - 0.00407272*ct - 0.748231*cA*cA*ct + 0.14861*cg*ct + 0.436721*ct*ct + 0.00637587*cV*cV + 0.000702175*cA*cA*cV*cV + 0.0002901*cg*cV*cV - 0.00040867*ct*cV*cV + 0.375159*pow(cA,4.) + 5.37945e-6*pow(cV,4.)) * wt;
    
    if (bin == 5) return (0.0485093 + 0.00891872*cA*cA - 0.00324976*cg + 0.155788*cA*cA*cg + 0.523622*cg*cg + 0.00396845*ct - 0.384512*cA*cA*ct - 0.0424597*cg*ct + 0.20254*ct*ct + 0.00271808*cV*cV + 0.000118455*cA*cA*cV*cV + 0.000624283*cg*cV*cV - 0.000523273*ct*cV*cV + 0.202822*pow(cA,4.) - 5.43403e-6 *pow(cV,4.)) * wt;
    
    if (bin == 6) return (0.0165142 + 0.00240494*cA*cA - 0.00151538*cg + 0.443147*cA*cA*cg + 2.01346*cg*cg + 0.0066562*ct - 0.285082*cA*cA*ct - 0.329009*cg*ct + 0.139285*ct*ct + 0.00226039*cV*cV + 0.000642254*cA*cA*cV*cV - 0.000373886*cg*cV*cV - 0.00122827*ct*cV*cV + 0.154372*pow(cA,4.) - 0.000211814*pow(cV,4.)) * wt;
    else {
        std::cout << "Bin not defined" << std::endl;
        return (EXIT_FAILURE);
    }
}

C_V::C_V(const StandardModel& SM_i)
: gg4l(SM_i)
{}

double C_V::computeThValue()
{
    updateParameters();
    
    return (1/2. - 4./3.*sw2)/(2. * (sqrt(sw2) * sqrt(1. - sw2))) * cV;
}

C_A::C_A(const StandardModel& SM_i)
: gg4l(SM_i)
{}

double C_A::computeThValue()
{
    updateParameters();
    
    return -1/2./(2. * (sqrt(sw2) * sqrt(1. - sw2))) * cA;
}
