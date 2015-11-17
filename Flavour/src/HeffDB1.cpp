/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDB1.h"
#include "gslpp_complex.h"

HeffDB1::HeffDB1(const StandardModel & SM) 
:       model(SM), 
        coeffnlep00qcd (10, NDR, NLO, NLO_ew), coeffnlep00 (12, NDR, NLO, NLO_ew),
        coeffnlep10qcd (10, NDR, NLO, NLO_ew), coeffnlep10 (12, NDR, NLO, NLO_ew),
        coeffnlep01 (10, NDR, NLO), coeffnlep01A(10, NDR, NLO), coeffnlep01B(4, NDR, NLO), coeffnlep00CC(10, NDR, NLO),
        coeffnlep11 (10, NDR, NLO), coeffnlep11A(10, NDR, NLO), coeffnlep11B(4, NDR, NLO), coeffnlep10CC(10, NDR, NLO),
        coeffsmumu (8, NDR, NNLO, NLO_ewt4), coeffdmumu (8, NDR, NNLO, NLO_ewt4),
        coeffbtaunu (3, NDR, LO),
        coeffsnunu (1, NDR, NLO), coeffdnunu (1, NDR, NLO),
        coeffsgamma(8,NDR, NLO),
        coeffprimesgamma(8,NDR, NLO),
        coeffBMll (13,NDR, NLO),
        coeffprimeBMll (13, NDR, NLO),
        evolDF1BMll(13, NDR, NLO, SM),
        evolDB1bsg(8, NDR, NLO, SM),
        u(10, NDR, NLO, NLO_ew, SM),
        evolbs(8, NDR, NNLO, NLO_ewt4, SM), evolbd(8, NDR, NNLO, NLO_ewt4, SM),
        nlep (12, 0.), nlep2(10, 0.),        
        nlepCC(4, 0.)
{
    
    for (unsigned int i = 0; i < 6; i++) {
        BMll_WC_cache.push_back(coeffBMll);
        BMll_Mu_cache.push_back(0.);
    }
    BMll_mu_cache = 0.;
    
    for (unsigned int i = 0; i < 6; i++) {
        BMllprime_WC_cache.push_back(coeffprimeBMll);
        BMllprime_Mu_cache.push_back(0.);
    }
    BMllprime_mu_cache = 0.;
    
    for (unsigned int i = 0; i < 6; i++) {
        Bsgamma_WC_cache.push_back(coeffsgamma);
        Bsgamma_Mu_cache.push_back(0.);
    }
    Bsgamma_mu_cache = 0.;
    
    for (unsigned int i = 0; i < 6; i++) {
        Bpsgamma_WC_cache.push_back(coeffsgamma);
        Bpsgamma_Mu_cache.push_back(0.);
    }
    
    for (unsigned int i = 0; i < 6; i++) {
        Bsmumu_WC_cache.push_back(coeffsmumu);
        Bsmumu_Mu_cache.push_back(0.);
    }
    Bsmumu_mu_cache = 0.;
    
    for (unsigned int i = 0; i < 6; i++) {
        Bdmumu_WC_cache.push_back(coeffdmumu);
        Bdmumu_Mu_cache.push_back(0.);
    }
    Bdmumu_mu_cache = 0.;
}

HeffDB1::~HeffDB1() 
{}

/*******************************************************************************
 * evolution Wilson Coefficients b-> nonlep.                                   * 
 * Buras base                                                                  *
 * deltaB=1   deltaC=0   deltaS=0                                              *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffBnlep00(double mu, schemes scheme) 
{
    
     std::vector<WilsonCoefficient>& mcb = model.getMyMatching()->CMbnlep( 0);
     std::vector<WilsonCoefficient>& mcbCC = model.getMyMatching()->CMbnlepCC( 0);
    
    coeffnlep00qcd.setMu(mu); //inizializes to zero the coefficients
    coeffnlep00CC.setMu(mu);
    coeffnlep00.setMu(mu);
    
    orders_ew ordDF1ew = coeffnlep00.getOrder_ew();
    orders ordDF1 =  coeffnlep00.getOrder();
    
    for (unsigned int i = 0; i < mcb.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++) {        
                
                //Evolves the LO terms and the ones proportional to alpha_s 
                coeffnlep00qcd.setCoeff(*coeffnlep00qcd.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), orders(k), NULL_ew, mcb[i].getScheme())*
                    (*(mcb[i].getCoeff(orders(j - k)))), orders(j));
                
                
                //Evolves terms proportional to alpha_e and alpha_e/aplha_s
                coeffnlep00qcd.setCoeff(*coeffnlep00qcd.getCoeff(orders_ew(j+4)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), NNLO, orders_ew(k+4), mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders(j - k)))), orders_ew(j+4));
                
                
                //Evolves the current*current part of the hamiltonian (the one non-proportional to lambda_t) 
                coeffnlep00CC.setCoeff(*coeffnlep00CC.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC[i].getMu(), orders(k), NULL_ew, mcbCC[i].getScheme()) *
                    (*(mcbCC[i].getCoeff(orders(j - k)))), orders(j)); 
   
            }
        }
        
                coeffnlep00qcd.setCoeff(*coeffnlep00qcd.getCoeff(orders_ew(NLO_ew)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), orders(LO), NULL_ew,  mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders_ew(NLO_ew)))), orders_ew(NLO_ew));       
        
    }
        
    coeffnlep00qcd.setScheme(scheme);
    coeffnlep00CC.setScheme(scheme);
    
    //Puts all together in a wilson coefficient object of 12 component: 
    //the first 10 terms are the ones proportional to lambda_t
    //the last two are the remainder current*current operators 
    for (int j=LO; j <= ordDF1; j++) {
        nlep2 = *coeffnlep00qcd.getCoeff(orders(j));
        for (int i = 0; i < 10; i++){
            nlep.assign(i, nlep2(i));
        }
        nlep2 = *coeffnlep00CC.getCoeff(orders(j));
        for (int i = 10; i < 12; i++){
            nlep.assign(i, nlep2(i-10));
        }        
    coeffnlep00.setCoeff(nlep, orders(j));
    }
    for (int k=LO_ew; k <= ordDF1ew; k++) {
        nlep2 = *coeffnlep00qcd.getCoeff(orders_ew(k));
        for (int l = 0; l < 10; l++){
            nlep.assign(l, nlep2(l));;
        }       
    coeffnlep00.setCoeff(nlep, orders_ew(k));
    }
    
    return coeffnlep00.getCoeff();
}

/*******************************************************************************
 * evolution Wilson Coefficienys b-> nonlep.                                   * 
 * Buras base                                                                  *
 * deltaB=1   deltaC=0   deltaS=1                                              *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffBnlep10(double mu, schemes scheme) 
{
    
     std::vector<WilsonCoefficient>& mcb = model.getMyMatching()->CMbnlep( 1);
     std::vector<WilsonCoefficient>& mcbCC = model.getMyMatching()->CMbnlepCC( 1);
    
    coeffnlep10qcd.setMu(mu);
    coeffnlep10CC.setMu(mu);
    coeffnlep10.setMu(mu);
    
    orders_ew ordDF1ew = coeffnlep10.getOrder_ew();
    orders ordDF1 =  coeffnlep10.getOrder();
    
    for (unsigned int i = 0; i < mcb.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++)   {
                
                //Evolves the LO terms and the ones proportional to alpha_s 
                coeffnlep10qcd.setCoeff(*coeffnlep10qcd.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), orders(k), NULL_ew,  mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders(j - k)))), orders(j));
                
                //Evolves terms proportional to alpha_e and alpha_e/aplha_s
                coeffnlep10qcd.setCoeff(*coeffnlep10qcd.getCoeff(orders_ew(j+4)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), NNLO, orders_ew(k+4), mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders(j - k)))), orders_ew(j+4));
                
                //Evolves the current*current part of the hamiltonian (the one non-proportional to lambda_t)
                coeffnlep10CC.setCoeff(*coeffnlep10CC.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC[i].getMu(), orders(k), NULL_ew, mcbCC[i].getScheme()) *
                    (*(mcbCC[i].getCoeff(orders(j - k)))), orders(j)); 
        
            }
        }
            
        coeffnlep10qcd.setCoeff(*coeffnlep10qcd.getCoeff(orders_ew(NLO_ew)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), orders(LO), NULL_ew, mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders(LO_ew)))), orders_ew(NLO_ew));
    }        
    
    coeffnlep10qcd.setScheme(scheme);
    coeffnlep10CC.setScheme(scheme);
    
    //Puts all together in a wilson coefficient object of 12 component: 
    //the first 10 terms are the one proportional to lambda_t
    //the last two are the remainder current*current operators 
    
    for (int j=LO; j <= ordDF1; j++) {
        nlep2 = *coeffnlep10qcd.getCoeff(orders(j));
        for (int i = 0; i < 10; i++){
            nlep.assign(i, nlep2(i));
        }
        nlep2 = *coeffnlep10CC.getCoeff(orders(j));
        for (int i = 10; i < 12; i++){
            nlep.assign(i, nlep2(i-10));
        }        
    coeffnlep10.setCoeff(nlep, orders(j));
    }
    for (int k=LO_ew; k <= ordDF1ew; k++) {
        nlep2 = *coeffnlep10qcd.getCoeff(orders_ew(k));
        for (int l = 0; l < 10; l++){
            nlep.assign(l, nlep2(l));;
        }       
    coeffnlep10.setCoeff(nlep, orders_ew(k));
    }
    
    return coeffnlep10.getCoeff();
}

/*******************************************************************************
 * evolution Wilson Coefficienys b-> nonlep.                                   * 
 * Buras base                                                                  *
 * deltaB=1   deltaC=1   deltaS=0                                              *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffBnlep01(double mu, schemes scheme) 
{
    
     std::vector<WilsonCoefficient>& mcbCC1 =model.getMyMatching()->CMbnlepCC( 2);
     std::vector<WilsonCoefficient>& mcbCC2 = model.getMyMatching()->CMbnlepCC( 3);
    
    coeffnlep01.setMu(mu);
    coeffnlep01A.setMu(mu);
    coeffnlep01B.setMu(mu);
    
    orders ordDF1 = coeffnlep01A.getOrder();
    
    //evolution of the current*current terms
    for (unsigned int i = 0; i < mcbCC1.size(); i++)
        for (int j = LO; j <= ordDF1; j++)
            for (int k = LO; k <= j; k++){
                
                coeffnlep01A.setCoeff(*coeffnlep01A.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC1[i].getMu(), orders(k), NULL_ew, mcbCC1[i].getScheme()) *
                    (*(mcbCC1[i].getCoeff(orders(j - k)))), orders(j)); 
                
                coeffnlep01B.setCoeff(*coeffnlep01B.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC2[i].getMu(), orders(k), NULL_ew, mcbCC2[i].getScheme()) *
                    (*(mcbCC2[i].getCoeff(orders(j - k)))), orders(j)); 
            }
        
    coeffnlep01A.setScheme(scheme);
    coeffnlep01B.setScheme(scheme);
    coeffnlep01.setScheme(scheme);
    
    //Puts all together in a wilson coefficient object of 4 components
    for (int j=LO; j <= ordDF1; j++) {
        nlep2 = *coeffnlep01A.getCoeff(orders(j));
        for (int i = 0; i < 2; i++){
            nlep.assign(i, nlep2(i));
        }
        nlep2 = *coeffnlep01A.getCoeff(orders(j));
        for (int i = 2; i < 4; i++){
            nlep.assign(i, nlep2(i-2));
        }        
    coeffnlep01.setCoeff(nlep, orders(j));
    }    
    return coeffnlep01.getCoeff();   
    
}

/*******************************************************************************
 * evolution Wilson Coefficienys b-> nonlep.                                   * 
 * Buras base                                                                  *
 * deltaB=1   deltaC=1   deltaS=1                                              *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffBnlep11(double mu, schemes scheme) 
{
    
     std::vector<WilsonCoefficient>& mcbCC1 = model.getMyMatching()->CMbnlepCC( 2);
     std::vector<WilsonCoefficient>& mcbCC2 = model.getMyMatching()->CMbnlepCC( 3);
    
    coeffnlep11.setMu(mu);
    coeffnlep11A.setMu(mu);
    coeffnlep11B.setMu(mu);
    
    orders ordDF1 = coeffnlep11A.getOrder();
    
    for (unsigned int i = 0; i < mcbCC1.size(); i++)
        for (int j = LO; j <= ordDF1; j++)
            for (int k = LO; k <= j; k++){
                
                coeffnlep11A.setCoeff(*coeffnlep11A.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC1[i].getMu(), orders(k), NULL_ew, mcbCC1[i].getScheme()) *
                    (*(mcbCC1[i].getCoeff(orders(j - k)))), orders(j)); 
        
                coeffnlep11B.setCoeff(*coeffnlep11B.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC2[i].getMu(), orders(k), NULL_ew, mcbCC2[i].getScheme()) *
                    (*(mcbCC2[i].getCoeff(orders(j - k)))), orders(j)); 
            }
        
    coeffnlep11A.setScheme(scheme);
    coeffnlep11B.setScheme(scheme);
    coeffnlep11.setScheme(scheme);
    
   for (int j=LO; j <= ordDF1; j++) {
        nlep2 = *coeffnlep11A.getCoeff(orders(j));
        for (int i = 0; i < 2; i++){
            nlep.assign(i, nlep2(i));
        }
        nlep2 = *coeffnlep11A.getCoeff(orders(j));
        for (int i = 2; i < 4; i++){
            nlep.assign(i, nlep2(i-2));
        }        
    coeffnlep11.setCoeff(nlep, orders(j));
    }    
    return coeffnlep11.getCoeff();   
    
}

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffsmumu(double mu, schemes scheme) 
{   
    
    
    coeffsmumu.setScheme(scheme);
    orders ordDF1 = coeffsmumu.getOrder(); 
    orders_ew ordDF1ew = coeffsmumu.getOrder_ew();
    
    const std::vector<WilsonCoefficient>& mc = model.getMyMatching() -> CMbsmm();
    
    if(mu == Bsmumu_mu_cache && scheme == Bsmumu_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mc.size(); i++){
            if (mc[i].getMu() == Bsmumu_Mu_cache[i]){
                
                
                for (int j = LO; j <= ordDF1; j++){
                    for (int k = LO; k <= j; k++){
                        for (int l = 0; l < 8; l++) {
                            check *= ((*(mc[i].getCoeff(orders(j - k))))(l) == (*(Bsmumu_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            }
        }
        check = 0;
        if (check == 1) return coeffsmumu.getCoeff();
    } 
       
    double M = 0;
    M = mc[0].getMu();
    double nf = 0;  
    nf = 5; //al the process has nf = 5, also the evolutor
    
    //int L = 6 - (int) nf;
    int j = 0;
    double alsM = evolbs.alphatilde_s(M);
    double alsmu = evolbs.alphatilde_s(mu);
    double eta = alsM / alsmu;
      
    double B00S = model.Beta0(nf), B10S = model.Beta1(nf); 

    double B00E = 80./9., B01E = 176./9.; 
    
    double logeta = log(eta);
    double fatt = (B00E * B10S /B00S /B00S - B01E /B00S);
    double app = B00E * (1. - eta)/ B00S;
    
    Bsmumu_mu_cache = mu;
    Bsmumu_scheme_cache = scheme;
    Bsmumu_WC_cache.clear();
    Bsmumu_WC_cache = mc;
   
    coeffsmumu.setMu(mu); 
           
    for (unsigned int i = 0; i < mc.size(); i++){
        Bsmumu_Mu_cache[i] = mc[i].getMu();
        for (j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){
                
                             
                
                if ((k <= NNLO) && (j <= NNLO)){
                    
                                        
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders(j)) + pow(eta,j) *
                    evolbs.Df1Evol(mu, mc[i].getMu(), orders(k), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
               
               
                }     
            }    
        }       
            
            for (j = LO_ew; j <= ordDF1ew; j++){
                
            switch(j) {   
            
            case(NLO_ewt4):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_ew(j)) +
                    (evolbs.Df1Evol(mu, mc[i].getMu(), LO, orders_ew(NULL_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders_ew(NLO_ewt4)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(NLO_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders_ew(NLO_ew)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(j), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(LO_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders_ew(NLO_ewt2)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(NLO_ewt1), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NNLO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(NLO_ewt3), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO))))) + pow((logeta) * fatt,2) *
                    (evolbs.Df1Evol(mu, mc[i].getMu(), LO, orders_ew(NULL_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO)))) ) + pow( app, 2 ) *
                    (evolbs.Df1Evol(mu, mc[i].getMu(), orders(LO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NNLO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(NULL_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NLO, orders_ew(NULL_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO)))) )+ logeta * fatt 
                    * app * 
                    (evolbs.Df1Evol(mu, mc[i].getMu(), orders(LO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NLO, orders_ew(NULL_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO))))), orders_ew(j));
                             
               break; 
            
            case(NLO_ewt3):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_ew(j)) + (1./ eta) *
                    (evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(LO_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders_ew(NLO_ew)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(NLO_ewt1), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(j), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO))))) +  ((logeta/eta) * fatt) 
                    * app * 
                    (evolbs.Df1Evol(mu, mc[i].getMu(), LO, orders_ew(NULL_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO))))) +  ( app * app/( eta) ) *
                    (evolbs.Df1Evol(mu, mc[i].getMu(), orders(LO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NLO, orders_ew(NULL_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO))))), orders_ew(j));
               
               break;
               
            case(NLO_ewt2):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_ew(j)) + eta *
                    (evolbs.Df1Evol(mu, mc[i].getMu(), orders(LO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders_ew(j)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), orders(NLO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders_ew(NLO_ew)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(NLO_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(j), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(LO_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NNLO)))))-
                    eta * logeta * fatt * 
                    (evolbs.Df1Evol(mu, mc[i].getMu(), orders(LO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), orders(NLO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO)))))
                    - eta *  app *
                    (evolbs.Df1Evol(mu, mc[i].getMu(), orders(LO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NNLO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), orders(NNLO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO)))) +
                    evolbs.Df1Evol(mu, mc[i].getMu(), orders(NLO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO)))) ), orders_ew(j));
               
               break;
            
            case(NLO_ewt1):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_ew(j)) + (1./ pow(eta,2)) *
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(j), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO)))) +  ( app 
                    * app/( eta * eta ) ) * 
                    evolbs.Df1Evol(mu, mc[i].getMu(), LO, orders_ew(NULL_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO)))), orders_ew(j));  
               
               break;
            
            case(NLO_ew):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_ew(j)) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(j), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), orders(LO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders_ew(j)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(LO_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO)))) - logeta * fatt * 
                    (evolbs.Df1Evol(mu, mc[i].getMu(), orders(LO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO))))) - app * 
                    (evolbs.Df1Evol(mu, mc[i].getMu(), LO, orders_ew(NULL_ew), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(NLO)))) + 
                    evolbs.Df1Evol(mu, mc[i].getMu(), orders(NLO), NULL_ew, mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(LO))))), orders_ew(j));
             
               
               break;
               
            case(LO_ew):
             
                coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_ew(j)) + (1./ pow(eta,1)) *
                evolbs.Df1Evol(mu, mc[i].getMu(), NNLO, orders_ew(j), mc[i].getScheme()) *
                (*(mc[i].getCoeff(orders(LO)))) - ( app /( eta ) ) *
                ( evolbs.Df1Evol(mu, mc[i].getMu(), LO, orders_ew(NULL_ew), mc[i].getScheme()) *
                (*(mc[i].getCoeff(orders(LO))))), orders_ew(j));
                
                
                
                break;
            }                
        }
    }
    return coeffsmumu.getCoeff();
}


gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffdmumu(double mu, schemes scheme) 
{   
    
    coeffdmumu.setScheme(scheme);
    orders ordDF1 = coeffdmumu.getOrder(); 
    orders_ew ordDF1ew = coeffdmumu.getOrder_ew();
    
    const std::vector<WilsonCoefficient>& mcbdm = model.getMyMatching() -> CMbdmm();
    
    if(mu == Bdmumu_mu_cache && scheme == Bdmumu_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mcbdm.size(); i++){
            if (mcbdm[i].getMu() == Bdmumu_Mu_cache[i]){
                
                
                for (int j = LO; j <= ordDF1; j++){
                    for (int k = LO; k <= j; k++){
                        for (int l = 0; l < 8; l++) {
                            check *= ((*(mcbdm[i].getCoeff(orders(j - k))))(l) == (*(Bdmumu_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            }
        }
        check = 0;
        if (check == 1) return coeffdmumu.getCoeff();
    } 
       
    double M = 0;
    M = mcbdm[0].getMu();
    double nf = 0;  
    nf = 5; //al the process has nf = 5, also the evolutor
    
    //int L = 6 - (int) nf;
    int j = 0;
    double alsM = evolbd.alphatilde_s(M);
    double alsmu = evolbd.alphatilde_s(mu);
    double eta = alsM / alsmu;
      
    double B00S = model.Beta0(nf), B10S = model.Beta1(nf); 

    double B00E = 80./9., B01E = 176./9.; 
    
    double logeta = log(eta);
    double fatt = (B00E * B10S /B00S /B00S - B01E /B00S);
    double app = B00E * (1. - eta)/ B00S;
    
    Bdmumu_mu_cache = mu;
    Bdmumu_scheme_cache = scheme;
    Bdmumu_WC_cache.clear();
    Bdmumu_WC_cache = mcbdm;
   
    coeffdmumu.setMu(mu); 
           
    for (unsigned int i = 0; i < mcbdm.size(); i++){
        Bdmumu_Mu_cache[i] = mcbdm[i].getMu();
        for (j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){
                
                             
                
                if ((k <= NNLO) && (j <= NNLO)){
                    
                                        
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders(j)) + pow(eta,j) *
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(k), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(j - k)))), orders(j));
               
               
                }     
            }    
        }       
            
            for (j = LO_ew; j <= ordDF1ew; j++){
                
            switch(j) {   
            
            case(NLO_ewt4):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_ew(j)) +
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), LO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_ew(NLO_ewt4)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(NLO_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_ew(NLO_ew)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(LO_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_ew(NLO_ewt2)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(NLO_ewt1), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NNLO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(NLO_ewt3), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO))))) + pow((logeta) * fatt,2) *
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), LO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) ) + pow( app, 2 ) *
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NNLO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NLO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) )+ logeta * fatt 
                    *  app * 
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NLO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))), orders_ew(j));
                             
               break; 
            
            case(NLO_ewt3):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_ew(j)) + (1./ eta) *
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(LO_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_ew(NLO_ew)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(NLO_ewt1), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))) +  ((logeta/eta) * fatt) 
                    * app *
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), LO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))) +  ( app * app/( eta) ) *
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NLO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))), orders_ew(j));
               
               break;
               
            case(NLO_ewt2):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_ew(j)) + eta *
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_ew(j)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(NLO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_ew(NLO_ew)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(NLO_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(LO_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NNLO)))))-
                    eta * logeta * fatt * 
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(NLO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))))
                    - eta *  app *
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NNLO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(NNLO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) +
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(NLO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) ), orders_ew(j));
               
               break;
            
            case(NLO_ewt1):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_ew(j)) + (1./ pow(eta,2)) *
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) +  ( app 
                    * app/(eta * eta ) ) * 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), LO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))), orders_ew(j));  
               
               break;
            
            case(NLO_ew):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_ew(j)) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_ew(j)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(LO_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) - logeta * fatt * 
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))) -  app * 
                    (evolbd.Df1Evol(mu, mcbdm[i].getMu(), LO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd.Df1Evol(mu, mcbdm[i].getMu(), orders(NLO), NULL_ew, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))), orders_ew(j));
             
               
               break;
               
            case(LO_ew):
             
                coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_ew(j)) + (1./ pow(eta,1)) *
                evolbd.Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_ew(j), mcbdm[i].getScheme()) *
                (*(mcbdm[i].getCoeff(orders(LO)))) - ( app /( eta ) ) *
                ( evolbd.Df1Evol(mu, mcbdm[i].getMu(), LO, orders_ew(NULL_ew), mcbdm[i].getScheme()) *
                (*(mcbdm[i].getCoeff(orders(LO))))), orders_ew(j));
                
                
                
                break;
            }                        
        }
    }   
    return coeffdmumu.getCoeff();
}

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffbtaunu() 
{
    std::vector<WilsonCoefficient>& mcb = model.getMyMatching() -> CMbtaunu();
    coeffbtaunu.resetCoefficient();
    orders ordDF1 = coeffbtaunu.getOrder();
    for (unsigned int i = 0; i < mcb.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            coeffbtaunu.setCoeff(*coeffbtaunu.getCoeff(orders(j))
                                     + *mcb[i].getCoeff(orders(j)), orders(j));
        }
    }
     return coeffbtaunu.getCoeff(); 
}

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffsnunu() 
{
    
     std::vector<WilsonCoefficient>& mcb = model.getMyMatching() -> CMBXsnn();
    
    orders ordDF1 = coeffsnunu.getOrder();
    
    for (unsigned int i = 0; i < mcb.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            coeffsnunu.setCoeff(*coeffsnunu.getCoeff(orders(j))
                                    + *mcb[i].getCoeff(orders(j)), orders(j));
        }
    }
     
    return coeffsnunu.getCoeff(); 
} 

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffdnunu() 
{
    
     std::vector<WilsonCoefficient>& mcb = model.getMyMatching() -> CMBXdnn();
    
    orders ordDF1 = coeffdnunu.getOrder();
    
    for (unsigned int i = 0; i < mcb.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            coeffdnunu.setCoeff(*coeffdnunu.getCoeff(orders(j))
                                    + *mcb[i].getCoeff(orders(j)), orders(j));
        }
    }
    
    return coeffdnunu.getCoeff(); 
} 

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffsgamma(double mu, schemes scheme) 
{
    
    coeffsgamma.setScheme(scheme);
    orders ordDF1 = coeffsgamma.getOrder();   
    
    const std::vector<WilsonCoefficient>& mcbsg = model.getMyMatching() -> CMbsg();
    
    if(mu == Bsgamma_mu_cache && scheme == Bsgamma_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mcbsg.size(); i++){
            if (mcbsg[i].getMu() == Bsgamma_Mu_cache[i]){
                for (int j = LO; j <= ordDF1; j++){
                    for (int k = LO; k <= j; k++){
                        for (int l = 0; l < 8; l++) {
                            check *= ((*(mcbsg[i].getCoeff(orders(j - k))))(l) == (*(Bsgamma_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            }
        }
        if (check == 1) return coeffsgamma.getCoeff();
    } 
    
    Bsgamma_mu_cache = mu;
    Bsgamma_scheme_cache = scheme;
    Bsgamma_WC_cache.clear();
    Bsgamma_WC_cache = mcbsg;
    
    coeffsgamma.setMu(mu); 
    
    for (unsigned int i = 0; i < mcbsg.size(); i++){
        Bsgamma_Mu_cache[i] = mcbsg[i].getMu();
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){
                coeffsgamma.setCoeff(*coeffsgamma.getCoeff(orders(j)) +
                    evolDB1bsg.Df1Evolbsg(mu, mcbsg[i].getMu(), orders(k), mcbsg[i].getScheme()) *
                    (*(mcbsg[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
    return coeffsgamma.getCoeff(); 
}

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffprimesgamma(double mu, schemes scheme) 
{
    
    coeffprimesgamma.setScheme(scheme);
    orders ordDF1 = coeffprimesgamma.getOrder();   
    
    const std::vector<WilsonCoefficient>& mcbsgp = model.getMyMatching() -> CMprimebsg();
    
    if(mu == Bsgamma_mu_cache && scheme == Bsgamma_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mcbsgp.size(); i++){
            if (mcbsgp[i].getMu() == Bpsgamma_Mu_cache[i]){
                for (int j = LO; j <= ordDF1; j++){
                    for (int k = LO; k <= j; k++){
                        for (int l = 0; l < 8; l++) {
                            check *= ((*(mcbsgp[i].getCoeff(orders(j - k))))(l) == (*(Bpsgamma_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            }
        }
        if (check == 1) return coeffprimesgamma.getCoeff();
    } 
    
    Bsgamma_mu_cache = mu;
    Bsgamma_scheme_cache = scheme;
    Bpsgamma_WC_cache.clear();
    Bpsgamma_WC_cache = mcbsgp;
    
    coeffprimesgamma.setMu(mu); 
    
    for (unsigned int i = 0; i < mcbsgp.size(); i++){
        Bpsgamma_Mu_cache[i] = mcbsgp[i].getMu();
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){
                coeffprimesgamma.setCoeff(*coeffprimesgamma.getCoeff(orders(j)) +
                    evolDB1bsg.Df1Evolbsg(mu, mcbsgp[i].getMu(), orders(k), mcbsgp[i].getScheme()) *
                    (*(mcbsgp[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
    return coeffprimesgamma.getCoeff(); 
}

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffBMll(double mu, schemes scheme) 
{
    
    coeffBMll.setScheme(scheme);
    orders ordDF1 = coeffBMll.getOrder();   
    
    const std::vector<WilsonCoefficient>& mc = model.getMyMatching() -> CMBMll();
    
    if(mu == BMll_mu_cache && scheme == BMll_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mc.size(); i++){
            if (mc[i].getMu() == BMll_Mu_cache[i]){
                for (int j = LO; j <= ordDF1; j++){
                    for (int k = LO; k <= j; k++){
                        for (int l = 0; l < 13; l++) {
                            check *= ((*(mc[i].getCoeff(orders(j - k))))(l) == (*(BMll_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            }
        }
        if (check == 1) return coeffBMll.getCoeff();
    } 
    
    BMll_mu_cache = mu;
    BMll_scheme_cache = scheme;
    BMll_WC_cache.clear();
    BMll_WC_cache = mc;
    
    coeffBMll.setMu(mu); 
    
    for (unsigned int i = 0; i < mc.size(); i++){
        BMll_Mu_cache[i] = mc[i].getMu();
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){
                coeffBMll.setCoeff(*coeffBMll.getCoeff(orders(j)) +
                    evolDF1BMll.Df1EvolMll(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
    return coeffBMll.getCoeff();
}


gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffprimeBMll(double mu, schemes scheme) 
{
    
    coeffprimeBMll.setScheme(scheme);
    orders ordDF1 = coeffprimeBMll.getOrder();  
    
    const std::vector<WilsonCoefficient>& mc = model.getMyMatching() -> CMprimeBMll();
    
    if(mu == BMllprime_mu_cache && scheme == BMllprime_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mc.size(); i++){
            if (mc[i].getMu() == BMllprime_Mu_cache[i]){
                for (int j = LO; j <= ordDF1; j++){
                    for (int k = LO; k <= j; k++){
                        for (int l = 0; l < 13; l++) {
                            check *= ((*(mc[i].getCoeff(orders(j - k))))(l) == (*(BMllprime_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            }
        }
        if (check == 1) return coeffprimeBMll.getCoeff();
    }
    
    BMllprime_mu_cache = mu;
    BMllprime_scheme_cache = scheme;
    BMllprime_WC_cache.clear();
    BMllprime_WC_cache = mc;
    
    coeffprimeBMll.setMu(mu); 
    
    for (unsigned int i = 0; i < mc.size(); i++){
        BMllprime_Mu_cache[i] = mc[i].getMu();
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){
                coeffprimeBMll.setCoeff(*coeffprimeBMll.getCoeff(orders(j)) +
                    evolDF1BMll.Df1EvolMll(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
   
    return coeffprimeBMll.getCoeff();
}
