/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HeffDF2.h"
#include <QCD.h>

HeffDF2::HeffDF2(const StandardModel& SM)
:       model(SM),
        drNDRLRI(5, 5, 0),
        coeffbd(5, NDR, NLO),
        coeffbs(5, NDR, NLO),
        coeffDd(5, NDR, NLO),
        coeffk(5, NDR, NLO),
        coeffmk(5, NDR, NLO),
        evolDF2(5, NDR, NLO, SM)
{
    
    double Nc = SM.getNc();
    drNDRLRI(0,0) = -(((-1. + Nc) * (-7. + log(4096.))) / Nc);
    drNDRLRI(1,1) = (-2. * (-1. + 6. * Nc * Nc - 8. * log(2.) + Nc * (-13. + log(1024.)))) / (3. * Nc);
    drNDRLRI(1,2) = (-2. * (13. - 10. * log(2.) + Nc * (-5. + log(256.)))) / (3. * Nc);
    drNDRLRI(2,1) = (-8. + 6. * Nc * Nc + 20. * log(2.) - 8. * Nc * (1. + log(4.))) / (3. * Nc);
    drNDRLRI(2,2) = (2. * (4. + Nc - 10. * Nc * log(2.) + log(256.))) / (3. * Nc);
    drNDRLRI(3,3) = (2. - 4. * Nc * Nc + log(4.)) / Nc;
    drNDRLRI(3,4) = 2. - log(4.);
    drNDRLRI(4,3) = -2. * (1. + log(2.));
    drNDRLRI(4,4) = (2. + log(4.)) / Nc;
}

HeffDF2::~HeffDF2() 
{}

gslpp::vector<gslpp::complex>** HeffDF2::ComputeCoeffBd(double mu, schemes scheme) 
{

     std::vector<WilsonCoefficient>& mc = model.getMyMatching()->CMdbd2();
    
    coeffbd.setMu(mu);
    
    coeffbd.setScheme(mc[0].getScheme());
    
    orders ordDF2 = coeffbd.getOrder();
    for (unsigned int i = 0; i < mc.size(); i++){
        ChangeScheme(mc[0].getScheme(),mc[i],ordDF2);
        for (int j = LO; j <= ordDF2; j++){
            for (int k = LO; k <= j; k++){                
                coeffbd.setCoeff(*coeffbd.getCoeff(orders(j)) +
                    evolDF2.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
    ChangeScheme(scheme, coeffbd, ordDF2);

    return coeffbd.getCoeff();
}

gslpp::vector<gslpp::complex>** HeffDF2::ComputeCoeffBs(double mu, schemes scheme) 
{

     std::vector<WilsonCoefficient>& mc = model.getMyMatching()->CMdbs2();

    coeffbs.setMu(mu);
    
    coeffbs.setScheme(mc[0].getScheme());

    orders ordDF2 = coeffbs.getOrder();
    for (unsigned int i = 0; i < mc.size(); i++){
        ChangeScheme(mc[0].getScheme(),mc[i],ordDF2);
        for (int j = LO; j <= ordDF2; j++){
            for (int k = LO; k <= j; k++){
                coeffbs.setCoeff(*coeffbs.getCoeff(orders(j)) +
                    evolDF2.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }

    ChangeScheme(scheme, coeffbs, ordDF2);

    return coeffbs.getCoeff();
}

gslpp::vector<gslpp::complex>** HeffDF2::ComputeCoeffdd(double mu, schemes scheme) 
{

     std::vector<WilsonCoefficient>& mc = model.getMyMatching()->CMdd2();
    
    coeffDd.setMu(mu);

    coeffDd.setScheme(mc[0].getScheme());

    orders ordDF2 = coeffDd.getOrder();
    for (unsigned int i = 0; i < mc.size(); i++){
       ChangeScheme(mc[0].getScheme(),mc[i],ordDF2);
       for (int j = LO; j <= ordDF2; j++){
            for (int k = LO; k <= j; k++){
                coeffDd.setCoeff(*coeffDd.getCoeff(orders(j)) +
                    evolDF2.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }

    ChangeScheme(scheme, coeffDd, ordDF2);
    
    return coeffDd.getCoeff();
}

gslpp::vector<gslpp::complex>** HeffDF2::ComputeCoeffK(double mu, schemes scheme) 
{

    std::vector<WilsonCoefficient>& mc = model.getMyMatching()->CMdk2();
    gslpp::vector<gslpp::complex> zero(5,0.);
    
    coeffk.setScheme(mc[0].getScheme());

    coeffk.setMu(mu);
    coeffk.setCoeff(zero,LO);
    coeffk.setCoeff(zero,NLO);

    orders ordDF2 = coeffk.getOrder();
    for (unsigned int i = 0; i < mc.size(); i++){
        if (i == 0){
            coeffk.setCoeff(0, evolDF2.etatt(mu) * model.getMyMatching()->S0tt()
                             + evolDF2.etacc(mu) * model.getMyMatching()->S0c()
                             + evolDF2.etact(mu) * model.getMyMatching()->S0ct()
                             + evolDF2.Df2Evol(mu, model.getMuw(), LO, mc[0].getScheme())(0,0) * model.getMyMatching()->ZDPtt()
                             + model.getMyMatching()->ZDPct(),
                            NLO);
#if SUSYFIT_DEBUG & 2
    std::cout << "mu = " << mu<< ", S0tt = " << model.GetMyMatching()->S0tt() << 
            ", S0cc = " << model.GetMyMatching()->S0c() << 
            ", S0ct = " << model.GetMyMatching()->S0ct() << std::endl <<
            ", etatt = " << evolDF2.etatt(mu) << 
            ", etacc = " << evolDF2.etacc(mu) << 
            ", etact = " << evolDF2.etact(mu) << std::endl <<
            "tt = " << evolDF2.etatt(mu)*model.GetMyMatching()->S0tt() << 
            ", cc = " << evolDF2.etacc(mu)*model.GetMyMatching()->S0c() << 
            ", ct = " << evolDF2.etact(mu)*model.GetMyMatching()->S0ct() << std::endl;
#endif
            
        }
        else {
            ChangeScheme(mc[0].getScheme(),mc[i],ordDF2);
            for (int j = LO; j <= ordDF2; j++){
                for (int k = LO; k <= j; k++){
                    coeffk.setCoeff(*coeffk.getCoeff(orders(j)) +
                        evolDF2.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                        (*(mc[i].getCoeff(orders(j - k)))), orders(j));
                }
            }
        }
    }
    
    ChangeScheme(scheme, coeffk, ordDF2);
    
    return coeffk.getCoeff();
}


gslpp::vector<gslpp::complex>** HeffDF2::ComputeCoeffmK(double mu, schemes scheme) 
{

     std::vector<WilsonCoefficient>& mc = model.getMyMatching()->CMdk2();
    gslpp::vector<gslpp::complex> zero(5,0.);
    
    coeffmk.setMu(mu);

    orders ordDF2 = coeffmk.getOrder();
    for (unsigned int i = 0; i < mc.size(); i++){
        if (i == 0){
            coeffmk.setCoeff(zero, NLO);
            coeffmk.setCoeff(zero, LO);
        }
        else {
            for (int j = LO; j <= ordDF2; j++){
                for (int k = LO; k <= j; k++){
                    coeffmk.setCoeff(*coeffmk.getCoeff(orders(j)) +
                        evolDF2.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                        (*(mc[i].getCoeff(orders(j - k)))), orders(j));
                }
            }
        }
    }
    
    ChangeScheme(scheme, coeffmk, ordDF2);
    
    return coeffmk.getCoeff();
}

void HeffDF2::ChangeScheme(schemes schout, WilsonCoefficient& c_in, orders order) 
{
    schemes schin = c_in.getScheme();
    if (schout == schin || order == LO) return;
    WilsonCoefficient c_out(5, schout, order);
    switch(schin) {
        case NDR:
            switch(schout) {
                case LRI:
                        c_out.setCoeff(*c_in.getCoeff(NLO) -
                            model.Als(c_in.getMu()) / 4. / M_PI * drNDRLRI.transpose()*
                            (*c_in.getCoeff(LO)), NLO);
                    c_in.setCoeff(*c_out.getCoeff(NLO),NLO);
                    c_in.setScheme(schout);
                    break;
                 default:
                    throw std::runtime_error("HeffDF2::ChangeScheme(): out scheme not implemented"); 
            }
            break;
        default:
            throw std::runtime_error("HeffDF2::ChangeScheme(): in scheme not implemented"); 
    }
}

