/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDB1.h"
#include "gslpp_complex.h"
#include "StandardModel.h"
#include "EvolDB1Mll.h"
#include "EvolDB1bsg.h"
#include "EvolBsmm.h"

HeffDB1::HeffDB1(const StandardModel & SM) 
:       model(SM), 
        coeffnlep00qcd (10, NDR, NLO, NLO_QED11), coeffnlep00 (12, NDR, NLO, NLO_QED11),
        coeffnlep10qcd (10, NDR, NLO, NLO_QED11), coeffnlep10 (12, NDR, NLO, NLO_QED11),
        coeffnlep01 (10, NDR, NLO), coeffnlep01A(10, NDR, NLO), coeffnlep01B(4, NDR, NLO), coeffnlep00CC(10, NDR, NLO),
        coeffnlep11 (10, NDR, NLO), coeffnlep11A(10, NDR, NLO), coeffnlep11B(4, NDR, NLO), coeffnlep10CC(10, NDR, NLO),
        coeffsmumu (8, NDR, NNLO, NLO_QED22), coeffdmumu (8, NDR, NNLO, NLO_QED22),
        coeffbtaunu (3, NDR, LO),
        coeffsnunu (1, NDR, NLO), coeffdnunu (1, NDR, NLO),
        coeffsgamma(8,NDR, NNLO),
        coeffprimesgamma(8,NDR, NNLO),
        coeffBMll (13,NDR, NLO),
        coeffprimeBMll (13, NDR, NLO),
        evolDF1BMll(new EvolDB1Mll(13, NDR, NLO, SM)),
        evolDB1bsg(new EvolDB1bsg(8, NDR, NNLO, SM)),
        evolbs(new EvolBsmm(8, NDR, NNLO, NLO_QED22, SM)), evolbd(new EvolBsmm(8, NDR, NNLO, NLO_QED22, SM)),
        nlep (12, 0.), nlep2(10, 0.),        
        nlepCC(4, 0.)
{
    
    for (unsigned int i = 0; i < 6; i++) {
        BMll_WC_cache.push_back(coeffBMll);
        BMll_Mu_cache.push_back(0.);
    }
    BMll_mu_cache = 0.;
    BMll_order_ini_cache = 0;
    
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
    Bsgamma_order_ini_cache = 0;
    
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

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffsmumu(double mu, schemes scheme) 
{   
    
    
    coeffsmumu.setScheme(scheme);
    orders ordDF1 = coeffsmumu.getOrder(); 
    orders_qed ordDF1ew = coeffsmumu.getOrder_qed();
    
    const std::vector<WilsonCoefficient>& mcbsm = model.getMatching().CMbsmm();

    if (mu == Bsmumu_mu_cache && scheme == Bsmumu_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mcbsm.size(); i++) {
            if (mcbsm[i].getMu() == Bsmumu_Mu_cache[i]) {
                for (int j = LO; j <= ordDF1; j++) {
                    for (int k = LO; k <= j; k++) {
                        for (int l = 0; l < 8; l++) {
                            check *= ((*(mcbsm[i].getCoeff(orders(j - k))))(l) == (*(Bsmumu_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            } else check = 0;
        } 
        if (check == 1) return coeffsmumu.getCoeff();
    } 
       
    double nf = 5; //al the process has nf = 5, also the evolutor
    
    //int L = 6 - (int) nf;
    int j = 0;
    double alsM = evolbs->alphatilde_s(mcbsm[0].getMu());
    double alsmu = evolbs->alphatilde_s(mu);
    double eta = alsM / alsmu;
      
    double B00S = model.Beta0(nf), B10S = model.Beta1(nf); 

    double B00E = 80./9., B01E = 176./9.; 
    
    double logeta = log(eta);
    double fatt = (B00E * B10S /B00S /B00S - B01E /B00S);
    double app = B00E * (1. - eta)/ B00S;
    
    Bsmumu_mu_cache = mu;
    Bsmumu_scheme_cache = scheme;
    Bsmumu_WC_cache.clear();
    Bsmumu_WC_cache = mcbsm;
   
    coeffsmumu.setMu(mu); 
           
    for (unsigned int i = 0; i < mcbsm.size(); i++){
        Bsmumu_Mu_cache[i] = mcbsm[i].getMu();
        for (j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){  
                if ((k <= NNLO) && (j <= NNLO)){
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders(j)) + pow(eta,j) *
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(k), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(j - k)))), orders(j));
                }     
            }    
        }       
            
            for (j = LO_QED; j <= ordDF1ew; j++){
                
            switch(j) {   
            
            case(NLO_QED22):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) +
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), LO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders_qed(NLO_QED22)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(j), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders_qed(NLO_QED21)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED02), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NNLO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED12), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO))))) + pow((logeta) * fatt,2) *
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), LO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO)))) ) + pow( app, 2 ) *
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NNLO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NLO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO)))) )+ logeta * fatt 
                    * app * 
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NLO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO))))), orders_qed(j));
                             
               break; 
            
            case(NLO_QED12):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) + (1./ eta) *
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED02), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(j), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO))))) +  ((logeta/eta) * fatt) 
                    * app * 
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), LO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO))))) +  ( app * app/( eta) ) *
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NLO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO))))), orders_qed(j));
               
               break;
               
            case(NLO_QED21):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) + eta *
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders_qed(j)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(NLO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(j), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NNLO)))))-
                    eta * logeta * fatt * 
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(NLO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO)))))
                    - eta *  app *
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NNLO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(NNLO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO)))) +
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(NLO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO)))) ), orders_qed(j));
               
               break;
            
            case(NLO_QED02):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) + (1./ pow(eta,2)) *
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(j), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO)))) +  ( app 
                    * app/( eta * eta ) ) * 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), LO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO)))), orders_qed(j));  
               
               break;
            
            case(NLO_QED11):
                
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(j), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders_qed(j)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO)))) - logeta * fatt * 
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO))))) - app * 
                    (evolbs->Df1Evol(mu, mcbsm[i].getMu(), LO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(NLO)))) + 
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(NLO), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(LO))))), orders_qed(j));
             
               
               break;
               
            case(LO_QED):
             
                coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) + (1./ pow(eta,1)) *
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(j), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))) - ( app /( eta ) ) *
                ( evolbs->Df1Evol(mu, mcbsm[i].getMu(), LO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO))))), orders_qed(j));
                
                
                
                break;
            }                
        }
    }
    return coeffsmumu.getCoeff();
}


gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffsmumuStandardNorm(double mu, schemes scheme) 
{      
    coeffsmumu.setScheme(scheme);
    orders ordDF1 = coeffsmumu.getOrder(); 
    orders_qed ordDF1ew = coeffsmumu.getOrder_qed();
    const std::vector<WilsonCoefficient>& mcbsm = model.getMatching().CMbsmm();

    if (mu == Bsmumu_mu_cache && scheme == Bsmumu_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mcbsm.size(); i++) {
            if (mcbsm[i].getMu() == Bsmumu_Mu_cache[i]) {
                for (int j = LO; j <= ordDF1; j++) {
                    for (int k = LO; k <= j; k++) {
                        for (int l = 0; l < 8; l++) {
                            check *= ((*(mcbsm[i].getCoeff(orders(j - k))))(l) == (*(Bsmumu_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            } else check = 0;
        } 
        if (check == 1) return coeffsmumu.getCoeff();
    } 
    double nf = 5; //al the process has nf = 5, also the evolutor
    //int L = 6 - (int) nf;
    int j = 0;
    double alsM = evolbs->alphatilde_s(mcbsm[0].getMu());
    double alsmu = evolbs->alphatilde_s(mu);
    double eta = alsM / alsmu;
      
    double B00S = model.Beta0(nf), B10S = model.Beta1(nf); 

    double B00E = 80./9., B01E = 176./9.; 
    
    double logeta = log(eta);
    double fatt = (B00E * B10S /B00S /B00S - B01E /B00S);
    double app = B00E * (1. - eta)/ B00S;
    
    Bsmumu_mu_cache = mu;
    Bsmumu_scheme_cache = scheme;
    Bsmumu_WC_cache.clear();
    Bsmumu_WC_cache = mcbsm;
   
    coeffsmumu.setMu(mu); 
           
    for (unsigned int i = 0; i < mcbsm.size(); i++){
        Bsmumu_Mu_cache[i] = mcbsm[i].getMu();
        for (j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){  
                if ((k <= NNLO) && (j <= NNLO)){
               coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders(j)) + pow(eta,j) *
                    evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(k), NO_QED, mcbsm[i].getScheme()) *
                    (*(mcbsm[i].getCoeff(orders(j - k)))), orders(j));
                }     
            }    
        }       
            
            for (j = LO_QED; j <= ordDF1ew; j++){
                
            switch(j) {   
            
            case(NLO_QED22):
                coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), LO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED22)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED22), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED21)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED02), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NNLO)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED12), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NLO)))) +
                logeta*fatt * (evolbs->Df1Evol(mu, mcbsm[i].getMu(), LO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NLO))))) +
                app * (evolbs->Df1Evol(mu, mcbsm[i].getMu(), LO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED21)))) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NLO, orders_qed(NO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NLO)))) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED21), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NNLO))))), orders_qed(j));
                break; 

            case(NLO_QED21):
                coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) +
                eta * (evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED21)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(NLO), NO_QED, mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NLO)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED21), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NNLO))))), orders_qed(j));
                break;
                
            case(NLO_QED12):
                coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) +
                (1./eta) * (evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED02), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NLO)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED12), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO))))) +
                (logeta*fatt/eta) * (evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO))))) +
                (app/eta) * (evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NLO))))), orders_qed(j));
                break;
            
            case(NLO_QED02):
                coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) +
                (1./ pow(eta,2)) * evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED02), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))) +
                (app/(eta*eta)) * evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))), orders_qed(j));
                break;
            
            case(NLO_QED11):
                coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), orders(LO), NO_QED, mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders_qed(NLO_QED11)))) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))) +
                evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(LO_QED), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(NLO)))), orders_qed(j));
                break;
               
            case(LO_QED):
                coeffsmumu.setCoeff(*coeffsmumu.getCoeff(orders_qed(j)) +
                (1./ eta) * evolbs->Df1Evol(mu, mcbsm[i].getMu(), NNLO, orders_qed(j), mcbsm[i].getScheme()) *
                (*(mcbsm[i].getCoeff(orders(LO)))), orders_qed(j));
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
    orders_qed ordDF1ew = coeffdmumu.getOrder_qed();
    
    const std::vector<WilsonCoefficient>& mcbdm = model.getMatching().CMbdmm();

    if (mu == Bdmumu_mu_cache && scheme == Bdmumu_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mcbdm.size(); i++) {
            if (mcbdm[i].getMu() == Bdmumu_Mu_cache[i]) {
                for (int j = LO; j <= ordDF1; j++) {
                    for (int k = LO; k <= j; k++) {
                        for (int l = 0; l < 8; l++) {
                            check *= ((*(mcbdm[i].getCoeff(orders(j - k))))(l) == (*(Bdmumu_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            } else check = 0;
        }
        if (check == 1) return coeffdmumu.getCoeff();
    } 
       
    double nf = 5; //al the process has nf = 5, also the evolutor
    
    //int L = 6 - (int) nf;
    int j = 0;
    double alsM = evolbd->alphatilde_s(mcbdm[0].getMu());
    double alsmu = evolbd->alphatilde_s(mu);
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
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(k), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(j - k)))), orders(j));
                }     
            }    
        }       
            
            for (j = LO_QED; j <= ordDF1ew; j++){
                
            switch(j) {   
            
            case(NLO_QED22):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_qed(j)) +
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), LO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_qed(NLO_QED22)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(LO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_qed(NLO_QED21)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(NLO_QED02), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NNLO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(NLO_QED12), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO))))) + pow((logeta) * fatt,2) *
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), LO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) ) + pow( app, 2 ) *
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NNLO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NLO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) )+ logeta * fatt 
                    *  app * 
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NLO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))), orders_qed(j));
                             
               break; 
            
            case(NLO_QED12):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_qed(j)) + (1./ eta) *
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(LO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(NLO_QED02), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))) +  ((logeta/eta) * fatt) 
                    * app *
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), LO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))) +  ( app * app/( eta) ) *
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NLO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))), orders_qed(j));
               
               break;
               
            case(NLO_QED21):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_qed(j)) + eta *
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_qed(j)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(NLO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_qed(NLO_QED11)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(NLO_QED11), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(LO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NNLO)))))-
                    eta * logeta * fatt * 
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(NLO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))))
                    - eta *  app *
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NNLO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(NNLO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) +
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(NLO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) ), orders_qed(j));
               
               break;
            
            case(NLO_QED02):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_qed(j)) + (1./ pow(eta,2)) *
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) +  ( app 
                    * app/(eta * eta ) ) * 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), LO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))), orders_qed(j));  
               
               break;
            
            case(NLO_QED11):
                
               coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_qed(j)) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(j), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders_qed(j)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(LO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) - logeta * fatt * 
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(LO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))) -  app * 
                    (evolbd->Df1Evol(mu, mcbdm[i].getMu(), LO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(NLO)))) + 
                    evolbd->Df1Evol(mu, mcbdm[i].getMu(), orders(NLO), NO_QED, mcbdm[i].getScheme()) *
                    (*(mcbdm[i].getCoeff(orders(LO))))), orders_qed(j));
             
               
               break;
               
            case(LO_QED):
             
                coeffdmumu.setCoeff(*coeffdmumu.getCoeff(orders_qed(j)) + (1./ pow(eta,1)) *
                evolbd->Df1Evol(mu, mcbdm[i].getMu(), NNLO, orders_qed(j), mcbdm[i].getScheme()) *
                (*(mcbdm[i].getCoeff(orders(LO)))) - ( app /( eta ) ) *
                ( evolbd->Df1Evol(mu, mcbdm[i].getMu(), LO, orders_qed(NO_QED), mcbdm[i].getScheme()) *
                (*(mcbdm[i].getCoeff(orders(LO))))), orders_qed(j));
                
                
                
                break;
            }                        
        }
    }   
    return coeffdmumu.getCoeff();
}

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffbtaunu(QCD::meson meson_i) 
{
    const std::vector<WilsonCoefficient>& mcb = model.getMatching().CMbtaunu(meson_i);
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
    
    const std::vector<WilsonCoefficient>& mcb = model.getMatching().CMBXsnn();
    
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
    
    const std::vector<WilsonCoefficient>& mcb = model.getMatching().CMBXdnn();
    
    orders ordDF1 = coeffdnunu.getOrder();
    
    for (unsigned int i = 0; i < mcb.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            coeffdnunu.setCoeff(*coeffdnunu.getCoeff(orders(j))
                                    + *mcb[i].getCoeff(orders(j)), orders(j));
        }
    }
    
    return coeffdnunu.getCoeff(); 
}

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffsgamma(double mu, bool noSM, schemes scheme)
{
    coeffsgamma.setScheme(scheme);
    orders ordDF1 = coeffsgamma.getOrder();

    unsigned int order_ini;
    if (noSM) order_ini = 1;
    else order_ini = 0;

    const std::vector<WilsonCoefficient>& mcbsg = model.getMatching().CMbsg();

    if (mu == Bsgamma_mu_cache && scheme == Bsgamma_scheme_cache && order_ini == BMll_order_ini_cache) {
        int check = 1;
        for (unsigned int i = order_ini; i < mcbsg.size(); i++) {
            if (mcbsg[i].getMu() == Bsgamma_Mu_cache[i]) {
                for (int j = LO; j <= ordDF1; j++) {
                    for (int k = LO; k <= j; k++) {
                        for (int l = 0; l < 8; l++) {
                            check *= ((*(mcbsg[i].getCoeff(orders(j - k))))(l) == (*(Bsgamma_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            } else check = 0;
        }
        if (check == 1) return coeffsgamma.getCoeff();
    }

    Bsgamma_mu_cache = mu;
    Bsgamma_order_ini_cache = order_ini;
    Bsgamma_scheme_cache = scheme;
    Bsgamma_WC_cache.clear();
    Bsgamma_WC_cache = mcbsg;

    coeffsgamma.setMu(mu);

    for (unsigned int i = order_ini; i < mcbsg.size(); i++) {
        Bsgamma_Mu_cache[i] = mcbsg[i].getMu();
        for (int j = LO; j <= ordDF1; j++) {
            for (int k = LO; k <= j; k++) {
                coeffsgamma.setCoeff(*coeffsgamma.getCoeff(orders(j)) +
                        evolDB1bsg->Df1Evolbsg(mu, mcbsg[i].getMu(), orders(k), mcbsg[i].getScheme()) *
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
    
    const std::vector<WilsonCoefficient>& mcbsgp = model.getMatching().CMprimebsg();

    if (mu == Bsgamma_mu_cache && scheme == Bsgamma_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mcbsgp.size(); i++) {
            if (mcbsgp[i].getMu() == Bpsgamma_Mu_cache[i]) {
                for (int j = LO; j <= ordDF1; j++) {
                    for (int k = LO; k <= j; k++) {
                        for (int l = 0; l < 8; l++) {
                            check *= ((*(mcbsgp[i].getCoeff(orders(j - k))))(l) == (*(Bpsgamma_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            } else check = 0;
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
                    evolDB1bsg->Df1Evolbsg(mu, mcbsgp[i].getMu(), orders(k), mcbsgp[i].getScheme()) *
                    (*(mcbsgp[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
    return coeffprimesgamma.getCoeff(); 
}

gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffBMll(double mu, QCD::lepton lepton, bool noSM, schemes scheme) 
{
    
    coeffBMll.setScheme(scheme);
    orders ordDF1 = coeffBMll.getOrder();   
    
    const std::vector<WilsonCoefficient>& mc = model.getMatching().CMBMll(lepton);

    unsigned int order_ini;
    if (noSM) order_ini = 1;
    else order_ini = 0;
    
    if (mu == BMll_mu_cache && scheme == BMll_scheme_cache && order_ini == BMll_order_ini_cache) {
        int check = 1;
        for (unsigned int i = order_ini; i < mc.size(); i++) {
            if (mc[i].getMu() == BMll_Mu_cache[i]) {
                for (int j = LO; j <= ordDF1; j++) {
                    for (int k = LO; k <= j; k++) {
                        for (int l = 0; l < 13; l++) {
                            check *= ((*(mc[i].getCoeff(orders(j - k))))(l) == (*(BMll_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            } else check = 0;
        }
        if (check == 1) return coeffBMll.getCoeff();
    }

    BMll_mu_cache = mu;
    BMll_order_ini_cache = order_ini;
    BMll_scheme_cache = scheme;
    BMll_WC_cache.clear();
    BMll_WC_cache = mc;

    coeffBMll.setMu(mu);

    for (unsigned int i = order_ini; i < mc.size(); i++) {
        BMll_Mu_cache[i] = mc[i].getMu();
        for (int j = LO; j <= ordDF1; j++) {
            for (int k = LO; k <= j; k++) {
                coeffBMll.setCoeff(*coeffBMll.getCoeff(orders(j)) +
                        evolDF1BMll->Df1EvolMll(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                        (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
    return coeffBMll.getCoeff();
}


gslpp::vector<gslpp::complex>** HeffDB1::ComputeCoeffprimeBMll(double mu, QCD::lepton lepton, schemes scheme) 
{
    
    coeffprimeBMll.setScheme(scheme);
    orders ordDF1 = coeffprimeBMll.getOrder();  
    
    const std::vector<WilsonCoefficient>& mc = model.getMatching().CMprimeBMll(lepton);
    
    if (mu == BMllprime_mu_cache && scheme == BMllprime_scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mc.size(); i++) {
            if (mc[i].getMu() == BMllprime_Mu_cache[i]) {
                for (int j = LO; j <= ordDF1; j++) {
                    for (int k = LO; k <= j; k++) {
                        for (int l = 0; l < 13; l++) {
                            check *= ((*(mc[i].getCoeff(orders(j - k))))(l) == (*(BMllprime_WC_cache[i].getCoeff(orders(j - k))))(l));
                        }
                    }
                }
            } else check = 0;
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
                    evolDF1BMll->Df1EvolMll(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
   
    return coeffprimeBMll.getCoeff();
}
