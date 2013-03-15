/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HeffDS1.h"

HeffDS1::HeffDS1(const StandardModel & SM) :
        model(SM), u(10, NDR, NLO, NLO_ew, SM), uM(10, NDR, NLO, SM),
        coeffds1 (10, NDR, NLO, NLO_ew), 
        coeffds1p0nunu(1, NDR, NLO, NLO_ew), coeffds1ppnunu(1, NDR, NLO, NLO_ew),
        coeffds1mumu(1, NDR, NLO), coeffds1cc(10, NDR, NLO), 
        DS1cce(10, 0.), DS1cc(10, 0.){
}

HeffDS1::~HeffDS1() {
}

vector<complex>** HeffDS1::ComputeCoeffDS1PP(double mu, schemes scheme) {
    
    const std::vector<WilsonCoefficient>& mcb = model.GetMyMatching()->CMK();
    const std::vector<WilsonCoefficient>& mcbCC = model.GetMyMatching()->CMKCC();
    
    coeffds1.setMu(mu); //inizializes to zero the coefficients
    coeffds1cc.setMu(mu);
    
    orders ordDF1 =  coeffds1.getOrder();
    
    switch(scheme) {
    
        case NDR:
            for (int i = 0; i < mcb.size(); i++){
                                   
                //evolves until the chram trheshold
                if (i == 0){
                    for (int j = LO; j <= ordDF1; j++){
                        for (int k = LO; k <= j; k++){        

                            //Evolves the LO terms and the ones proportional to alpha_s 
                            coeffds1.setCoeff(*coeffds1.getCoeff(orders(j)) +
                                u.Df1Evolnlep(model.getMuc(), mcb[i].getMu(), 
                                orders(k), NULL_ew, mcb[i].getScheme())*
                                (*(mcb[i].getCoeff(orders(j - k)))), orders(j));

                            //Evolves terms proportional to alpha_e and alpha_e/aplha_s
                            coeffds1.setCoeff(*coeffds1.getCoeff(orders_ew(j+4)) +
                                u.Df1Evolnlep(model.getMuc(), mcb[i].getMu(), NNLO, 
                                orders_ew(k+4), mcb[i].getScheme()) *
                                (*(mcb[i].getCoeff(orders(j - k)))), orders_ew(j+4));

                            //Evolves the open-charm current*current part
                            coeffds1cc.setCoeff(*coeffds1cc.getCoeff(orders(j)) +
                                u.Df1Evolnlep(model.getMuc(), mcbCC[i].getMu(), 
                                orders(k), NULL_ew, mcbCC[i].getScheme()) *
                                (*(mcbCC[i].getCoeff(orders(j - k)))), orders(j));
                        }

                    coeffds1.setCoeff(*coeffds1.getCoeff(orders_ew(NLO_ew)) +
                        u.Df1Evolnlep(model.getMuc(), mcb[i].getMu(), orders(LO), 
                        NULL_ew,  mcb[i].getScheme()) *
                        (*(mcb[i].getCoeff(orders_ew(NLO_ew)))), orders_ew(NLO_ew));
                    }

                    //Matching at the charm threshold
                    CharmMatch();

                    //evolves below the chram trheshold
                    for (int j = LO; j <= ordDF1; j++){
                        for (int k = LO; k <= j; k++){  

                            coeffds1.setCoeff(*coeffds1.getCoeff(orders(j)) +
                                u.Df1Evolnlep(mu, model.getMuc(), orders(k),
                                NULL_ew, mcb[i].getScheme())*
                                (*(mcb[i].getCoeff(orders(j - k)))), orders(j));

                            coeffds1.setCoeff(*coeffds1.getCoeff(orders_ew(j+4)) +
                                u.Df1Evolnlep(mu, model.getMuc(), NNLO, 
                                orders_ew(k+4), mcb[i].getScheme()) *
                                (*(mcb[i].getCoeff(orders(j - k)))), orders_ew(j+4));
                        }
                    }

                    coeffds1.setCoeff(*coeffds1.getCoeff(orders_ew(NLO_ew)) +
                        u.Df1Evolnlep(mu, model.getMuc(), orders(LO), 
                        NULL_ew,  mcb[i].getScheme()) *
                        (*(mcb[i].getCoeff(orders_ew(NLO_ew)))), orders_ew(NLO_ew));
                }
                
                if (i>0){
                    //if (model.BasisFlag() == 0){
                    //evolves according to the Misiak Basis
                        for (int j = LO; j <= ordDF1; j++){
                            for (int k = LO; k <= j; k++){ 
                                coeffds1.setCoeff(*coeffds1.getCoeff(orders(j)) +
                                    uM.Df1Evolbsg(mu, mcb[i].getMu(), orders(k), mcb[i].getScheme())*
                                    (*(mcb[i].getCoeff(orders(j - k)))), orders(j));
                            }
                        }
                    /*}else{
                        for (int j = LO; j <= ordDF1; j++){
                            for (int k = LO; k <= j; k++){  

                                //Evolves the LO terms and the ones proportional to alpha_s 
                                coeffds1.setCoeff(*coeffds1.getCoeff(orders(j)) +
                                    u.Df1Evolnlep(mu, mcb[i].getMu(), 
                                    orders(k), NULL_ew, mcb[i].getScheme())*
                                    (*(mcb[i].getCoeff(orders(j - k)))), orders(j));

                                //Evolves terms proportional to alpha_e and alpha_e/aplha_s
                                coeffds1.setCoeff(*coeffds1.getCoeff(orders_ew(j+4)) +
                                    u.Df1Evolnlep(mu, mcb[i].getMu(), NNLO, 
                                    orders_ew(k+4), mcb[i].getScheme()) *
                                    (*(mcb[i].getCoeff(orders(j - k)))), orders_ew(j+4));
                            }
                        }
                        coeffds1.setCoeff(*coeffds1.getCoeff(orders_ew(NLO_ew)) +
                            u.Df1Evolnlep(mu, mcb[i].getMu(), orders(LO), 
                            NULL_ew,  mcb[i].getScheme()) *
                            (*(mcb[i].getCoeff(orders_ew(NLO_ew)))), orders_ew(NLO_ew));
                    }*/
                }
            }
            return coeffds1.getCoeff();

        default:
            throw "HeffDS1::ComputeCoeffDS1PP(double mu, schemes scheme): scheme not implemented";
    }
}

void HeffDS1::CharmMatch(){
    DS1cc = *coeffds1cc.getCoeff(LO);
    DS1cce = *coeffds1cc.getCoeff(LO_ew);
    
    double mc = model.Mrun(model.getMuc(), model.getQuarks(QCD::CHARM).getMass(), NNLO);

    DS1cc.assign(2, (-model.Als(model.getMuc())/24./M_PI)*(-2./3.*
                        (log(mc*mc/model.getMuc()/model.getMuc())+1.)*DS1cc(1)));
    DS1cc.assign(3, (model.Als(model.getMuc())/8./M_PI)*(-2./3.*
                        (log(mc*mc/model.getMuc()/model.getMuc())+1.)*DS1cc(1)));
    DS1cc.assign(4, (-model.Als(model.getMuc())/24./M_PI)*(-2./3.*
                        (log(mc*mc/model.getMuc()/model.getMuc())+1.)*DS1cc(1)));
    DS1cc.assign(5, (model.Als(model.getMuc())/8./M_PI)*(-2./3.*
                        (log(mc*mc/model.getMuc()/model.getMuc())+1.)*DS1cc(1)));
    DS1cc.assign(8, 0.);
    DS1cc.assign(7, 0.);
    DS1cc.assign(8, 0.);
    DS1cc.assign(9, 0.);
    
    DS1cc.assign(2, (-model.Als(model.getMuc())/24./M_PI)*(-2./3.*
                        (log(mc*mc/model.getMuc()/model.getMuc())+1.)*DS1cce(1)));
    DS1cc.assign(3, (model.Als(model.getMuc())/8./M_PI)*(-2./3.*
                        (log(mc*mc/model.getMuc()/model.getMuc())+1.)*DS1cce(1)));
    DS1cc.assign(4, (-model.Als(model.getMuc())/24./M_PI)*(-2./3.*
                        (log(mc*mc/model.getMuc()/model.getMuc())+1.)*DS1cce(1)));
    DS1cc.assign(5, (model.Als(model.getMuc())/8./M_PI)*(-2./3.*
                        (log(mc*mc/model.getMuc()/model.getMuc())+1.)*DS1cce(1)));
    DS1cce.assign(8, -model.getAle()/6./M_PI*4./9.*
                 (log(mc*mc/model.getMuc()/model.getMuc())+1.)*(3.*DS1cc(0)+DS1cc(1)));
    DS1cce.assign(7, 0.);
    DS1cce.assign(8, -model.getAle()/6./M_PI*4./9.*
                 (log(mc*mc/model.getMuc()/model.getMuc())+1.)*(3.*DS1cc(0)+DS1cc(1)));
    DS1cce.assign(9, 0.);
    DS1cce.assign(0, 0.);
    DS1cce.assign(1, 0.);

    coeffds1cc.setCoeff(DS1cc, NLO);
    coeffds1cc.setCoeff(DS1cce, NLO_ew);
    
    coeffds1.setCoeff(*coeffds1.getCoeff(NLO_ew) - *coeffds1cc.getCoeff(NLO_ew), NLO);

    coeffds1.setCoeff(*coeffds1.getCoeff(NLO) - *coeffds1cc.getCoeff(NLO), NLO);
    coeffds1.setCoeff(*coeffds1.getCoeff(LO) - *coeffds1cc.getCoeff(LO), LO);
    
    coeffds1cc.setMu(0.);
}
