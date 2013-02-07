/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HeffDF2.h"
#include <sstream>
#include <QCD.h>
#include <stdexcept>

HeffDF2::HeffDF2(const StandardModel& SM): model(SM), coeffbd(5, NDR, NLO), coeffbs(5, NDR, NLO), 
        coeffDd(5, NDR, NLO), coeffk(5, NDR, NLO), u(5, NDR, NLO, SM), drNDRLRI(5, 5, 0){
    
    double Nc = SM.getNc();
    drNDRLRI(0,0) = -(((-1. + Nc)*(-7. + log(4096.))) / Nc);                   
    drNDRLRI(1,1) = (-2.*(-1. + 6.*Nc*Nc - 8.*log(2.) + Nc*(-13. + log(1024.))))/(3.*Nc);
    drNDRLRI(1,2) = (-2.*(13. - 10.*log(2.) + Nc*(-5. + log(256.))))/(3.*Nc);   
    drNDRLRI(2,1) = (-8. + 6.*Nc*Nc + 20.*log(2.) - 8.*Nc*(1. + log(4.)))/(3.*Nc);
    drNDRLRI(2,2) = (2.*(4. + Nc - 10.*Nc*log(2.) + log(256.)))/(3.*Nc);
    drNDRLRI(3,3) = (2. - 4.*Nc*Nc + log(4.))/Nc;
    drNDRLRI(3,4) = 2. - log(4.);
    drNDRLRI(4,3) = -2.*(1. + log(2.));
    drNDRLRI(4,4) = (2. + log(4.))/Nc;
}

HeffDF2::~HeffDF2() {
}

vector<complex>** HeffDF2::ComputeCoeffBd(double mu, schemes scheme) {

    const std::vector<WilsonCoefficient>& mc = model.GetMyMatching()->CMdbd2();
    
    coeffbd.setMu(mu);
    
    orders ordDF2 = coeffbd.getOrder();
    for (int i = 0; i < mc.size(); i++){
        for (int j = LO; j <= ordDF2; j++){
            for (int k = LO; k <= j; k++){                
                coeffbd.setCoeff(*coeffbd.getCoeff(orders(j)) +
                    u.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
    ChangeScheme(scheme, mc[0].getScheme(), ordDF2, 0);

    coeffbd.setScheme(scheme);
    
    return coeffbd.getCoeff();
}

vector<complex>** HeffDF2::ComputeCoeffBs(double mu, schemes scheme) {

    const std::vector<WilsonCoefficient>& mc = model.GetMyMatching()->CMdbs2();

    coeffbs.setMu(mu);

    orders ordDF2 = coeffbs.getOrder();
    for (int i = 0; i < mc.size(); i++){
        for (int j = LO; j <= ordDF2; j++){
            for (int k = LO; k <= j; k++){
                coeffbs.setCoeff(*coeffbs.getCoeff(orders(j)) +
                    u.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }

    ChangeScheme(scheme, mc[0].getScheme(), ordDF2, 1);

    coeffbs.setScheme(scheme);

    return coeffbs.getCoeff();
}

vector<complex>** HeffDF2::ComputeCoeffdd(double mu, schemes scheme) {

    const std::vector<WilsonCoefficient>& mc = model.GetMyMatching()->CMdd2();
    
    coeffDd.setMu(mu);

    orders ordDF2 = coeffDd.getOrder();
    for (int i = 0; i < mc.size(); i++){
        for (int j = LO; j <= ordDF2; j++){
            for (int k = LO; k <= j; k++){
                coeffDd.setCoeff(*coeffDd.getCoeff(orders(j)) +
                    u.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }

    coeffDd.setScheme(scheme); 
    
    return coeffDd.getCoeff();
}

vector<complex>** HeffDF2::ComputeCoeffK(double mu, schemes scheme) {

    const std::vector<WilsonCoefficient>& mc = model.GetMyMatching()->CMdk2();
    vector<complex> zero(5,0.);
    
    coeffk.setMu(mu);

    orders ordDF2 = coeffk.getOrder();
    for (int i = 0; i < mc.size(); i++){
        if (i == 0){
            coeffk.setCoeff(0, u.etatt(mu) * model.GetMyMatching()->S0tt() +
                u.etacc(mu) * model.GetMyMatching()->S0c() + u.etact(mu) * model.GetMyMatching()->S0ct(), NLO);
            coeffk.setCoeff(zero, LO);
        }else{
            for (int j = LO; j <= ordDF2; j++){
                for (int k = LO; k <= j; k++){
                    coeffk.setCoeff(*coeffk.getCoeff(orders(j)) +
                        u.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                        (*(mc[i].getCoeff(orders(j - k)))), orders(j));
                }
            }
        }
    }
    
    coeffk.setScheme(scheme); 
    
    return coeffk.getCoeff();
}

void HeffDF2::ChangeScheme(schemes schout, schemes schin, orders order, int j) {
    if (schout == schin) return;
    switch(schin) {
        case NDR:
            switch(schout) {
                case LRI:
                    if (j == 0)
                        coeffbd.setCoeff(*coeffbd.getCoeff(NLO) -
                            model.Als(coeffbd.getMu()) / 4. / M_PI * drNDRLRI.transpose()*
                            (*coeffbd.getCoeff(LO)), NLO);
                    break;
                    if (j == 1)
                        coeffbs.setCoeff(*coeffbs.getCoeff(NLO) -
                            model.Als(coeffbs.getMu()) / 4. / M_PI * drNDRLRI.transpose()*
                            (*coeffbs.getCoeff(LO)), NLO);
                    break;
                    if (j == 2)
                        coeffk.setCoeff(*coeffk.getCoeff(NLO) -
                            model.Als(coeffk.getMu()) / 4. / M_PI * drNDRLRI.transpose()*
                            (*coeffk.getCoeff(LO)), NLO);
                    break;
                    if (j == 3)
                        coeffDd.setCoeff(*coeffDd.getCoeff(NLO) -
                            model.Als(coeffDd.getMu()) / 4. / M_PI * drNDRLRI.transpose()*
                            (*coeffDd.getCoeff(LO)), NLO);
                    break;
                default:
                    throw std::runtime_error("HeffDF2::ChangeScheme(): out scheme not implemented"); 
            }
        default:
            throw std::runtime_error("HeffDF2::ChangeScheme(): in scheme not implemented"); 
    }
}

