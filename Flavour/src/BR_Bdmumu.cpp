/* 
 * Copyright (C) 2012 SsuyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BR_Bdmumu.h"

double BR_Bdmumu::getThValue(){
    double mmu = myFlavour.getModel().getLeptons(StandardModel::MU).getMass();
    double mBd = myFlavour.getModel().getMesons(QCD::B_D).getMass();
    double FBd = myFlavour.getModel().getMesons(QCD::B_D).getDecayconst();
    double beta = sqrt(1. - pow(2. * mmu / mBd, 2.));
    double coupling = myFlavour.getModel().getGF() * myFlavour.getModel().getAle() / 4. / M_PI / myFlavour.getModel().sW2();
    double PRF = pow(coupling, 2.) / M_PI / myFlavour.getModel().getMesons(QCD::B_D).getWidth() * pow(FBd, 2.) * pow(mmu, 2.) * mBd * beta;
/*    double theta = asin(sqrt( (M_PI * myFlavour.getModel().getAle() )/( sqrt(2) * myFlavour.getModel().getGF() *
                   myFlavour.getModel().Mw_tree() * myFlavour.getModel().Mw_tree()) ));
    
   return( myFlavour.getModel().getMesons(QCD::B_D).Lifetime() * myFlavour.getModel().getGF()*myFlavour.getModel().getGF()/M_PI
            * myFlavour.getModel().getAle()*myFlavour.getModel().getAle()/(16.*M_PI*M_PI*pow(sin(theta),4.))
            * mBd * myFlavour.getModel().getMesons(QCD::B_D).getDecayconst()
            * mmu * mmu * sqrt(1.-4.*mmu*mmu/mBd/mBd) * BRBdmumu(NLO).real());*/
    double Bdmumu = PRF * BRBdmumu(FULLNLO).real();
    
    return(Bdmumu);
}

complex BR_Bdmumu::BRBdmumu(orders order){
    if (myFlavour.getHDB1().getCoeffdmumu().getOrder() < order % 3){
        std::stringstream out;
        out << order;
        throw std::runtime_error("BRBdmumu::getThValue(): required cofficient of "
                                 "order " + out.str() + " not computed");
    }
    
    vector<complex> ** allcoeff = myFlavour.ComputeCoeffdmumu();
    
    switch(order) {
        case FULLNLO:
                return((*(allcoeff[LO]) + *(allcoeff[NLO])).conjugate() *
                       (*(allcoeff[LO]) + *(allcoeff[NLO])));
            
        case LO:
                return((*(allcoeff[LO])).conjugate() *
                       (*(allcoeff[LO])));
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("BRBdmumu::BRBdmumu(): order " + out.str() + "not implemented");;
    }
}