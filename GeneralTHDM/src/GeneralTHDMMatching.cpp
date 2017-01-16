/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMMatching.h"
#include "GeneralTHDM.h"
#include <gsl/gsl_integration.h>

GeneralTHDMMatching::GeneralTHDMMatching(const GeneralTHDM & GeneralTHDM_i) :

    StandardModelMatching(GeneralTHDM_i),
    myGTHDM(GeneralTHDM_i),
    mcgminus2mu(2, NDR, NLO)

{}

void GeneralTHDMMatching::updateGTHDMParameters()
{
    GF=myGTHDM.getGF();
    mMU=myGTHDM.getLeptons(StandardModel::MU).getMass();
}

double GeneralTHDMMatching::gminus2muLO() {

    updateGTHDMParameters();
    double gminus2muLO;
    gminus2muLO=0.0;
    return(gminus2muLO);
}

double GeneralTHDMMatching::gminus2muNLO() {

    updateGTHDMParameters();
    double gminus2muNLO;
    gminus2muNLO=0.0;
    return(gminus2muNLO);
}

std::vector<WilsonCoefficient>& GeneralTHDMMatching::CMgminus2mu() {

    vmcgminus2mu = StandardModelMatching::CMgminus2mu();

    double gminus2muLOvalue=gminus2muLO();
    double gminus2muNLOvalue=gminus2muNLO();
    
    switch (mcgminus2mu.getOrder()) {
        case LO:
            mcgminus2mu.setCoeff(0, gminus2muLOvalue, LO);  //g-2_muR
            mcgminus2mu.setCoeff(1, 0., LO);  //g-2_muL
            break;
        case NLO:
            mcgminus2mu.setCoeff(0, gminus2muLOvalue+gminus2muNLOvalue, NLO);  //g-2_muR
            mcgminus2mu.setCoeff(1, 0., NLO);  //g-2_muL
            break;
        case NNLO:
        default:
            std::stringstream out;
            out << mcgminus2mu.getOrder();
            throw std::runtime_error("GeneralTHDMMatching::CMgminus2mu(): order " + out.str() + " not implemented.\nOnly leading order (LO) or next-to-leading order (NLO) are allowed.");
    }

    vmcgminus2mu.push_back(mcgminus2mu);
    return(vmcgminus2mu);

}
