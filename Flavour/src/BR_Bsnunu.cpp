/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Bsnunu.h"

double BR_Bsnunu::getThValue(){
    double theta= asin(sqrt( (M_PI * myFlavour.getModel().getAle() )/( sqrt(2) *
                  myFlavour.getModel().getGF() * myFlavour.getModel().Mw_tree() * 
                  myFlavour.getModel().Mw_tree()) ));
    
    double z = myFlavour.getModel().getQuarks(QCD::CHARM).getMass()/
               myFlavour.getModel().getQuarks(QCD::BOTTOM).getMass();
    
    
    return(3.*myFlavour.getModel().getAle()*myFlavour.getModel().getAle()/(4.*M_PI*M_PI*pow(sin(theta),4.)) * 
           myFlavour.getModel().getBr_B_Xcenu() *
           BRBsnunu(NLO).real() * (1.+2.*myFlavour.getModel().Als(myFlavour.getModel().getMub())/3./M_PI)*(25./4.-M_PI*M_PI)
           /(1. - 8.*z*z + 8.*z*z*z*z*z*z - z*z*z*z*z*z*z*z - 24.*z*z*z*z*log(z))
           /(1. - 2.*myFlavour.getModel().Als(myFlavour.getModel().getMub())/3./M_PI*((M_PI*M_PI-31./4.)*(1.-z*z)+1.5)));
}

complex BR_Bsnunu::BRBsnunu(orders order){
    if (myFlavour.getHDB1().getCoeffsnunu().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRBsnunu::getThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    vector<complex> ** allcoeff = myFlavour.ComputeCoeffsnunu();
    
    switch(order) {
        case NLO:
                return((*(allcoeff[LO]) + *(allcoeff[NLO])) *
                       (*(allcoeff[LO]) + *(allcoeff[NLO])));
        case LO:
                return((*(allcoeff[LO])) *
                       (*(allcoeff[LO])));
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("BRBsnunu::BRBsnunu(): order " + out.str() + "not implemented");;
    }
}

