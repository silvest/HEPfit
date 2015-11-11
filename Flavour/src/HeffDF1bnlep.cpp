/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDF1bnlep.h"
#include "gslpp.h"

HeffDF1bnlep::HeffDF1bnlep(const StandardModel & SM, StandardModelMatching & SM_Matching) 
:       model(SM), 
        coeffbnlep00qcd (10, NDR, NLO, NLO_ew), coeffbnlep00 (12, NDR, NLO, NLO_ew),
        coeffbnlep10qcd (10, NDR, NLO, NLO_ew), coeffbnlep10 (12, NDR, NLO, NLO_ew),
        coeffbnlep01 (10, NDR, NLO), coeffbnlep01A(10, NDR, NLO), coeffbnlep01B(4, NDR, NLO), coeffbnlep00CC(10, NDR, NLO),
        coeffbnlep11 (10, NDR, NLO), coeffbnlep11A(10, NDR, NLO), coeffbnlep11B(4, NDR, NLO), coeffbnlep10CC(10, NDR, NLO),
        u(10, NDR, NLO, NLO_ew, SM), 
        bnlep (12, 0.), bnlep2(10, 0.), bnlepCC(4, 0.)
{}

HeffDF1bnlep::~HeffDF1bnlep() 
{}

/*******************************************************************************
 * evolution Wilson Coefficienys b-> nonlep.                                   * 
 * Buras base                                                                  *
 * deltaB=1   deltaC=0   deltaS=0                                              *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDF1bnlep::ComputeCoeffBnlep00(double mu, schemes scheme) 
{
    
    const std::vector<WilsonCoefficient>& mcb = model.getMyMatching()->CMbnlep( 0);
    const std::vector<WilsonCoefficient>& mcbCC = model.getMyMatching()->CMbnlepCC( 0);
    
    coeffbnlep00qcd.setMu(mu); //inizializes to zero the coefficients
    coeffbnlep00CC.setMu(mu);
    coeffbnlep00.setMu(mu);
    
    orders_ew ordDF1ew = coeffbnlep00.getOrder_ew();
    orders ordDF1 =  coeffbnlep00.getOrder();
    
    for (unsigned int i = 0; i < mcb.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++) {        
                
                //Evolves the LO terms and the ones proportional to alpha_s 
                coeffbnlep00qcd.setCoeff(*coeffbnlep00qcd.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), orders(k), NULL_ew, mcb[i].getScheme())*
                    (*(mcb[i].getCoeff(orders(j - k)))), orders(j));
                
                
                //Evolves terms proportional to alpha_e and alpha_e/aplha_s
                coeffbnlep00qcd.setCoeff(*coeffbnlep00qcd.getCoeff(orders_ew(j+4)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), NNLO, orders_ew(k+4), mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders(j - k)))), orders_ew(j+4));
   
            }
        }
        
                coeffbnlep00qcd.setCoeff(*coeffbnlep00qcd.getCoeff(orders_ew(NLO_ew)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), orders(LO), NULL_ew,  mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders_ew(NLO_ew)))), orders_ew(NLO_ew));       
        
    }     
    
    //Evolves the current*current part of the hamiltonian (the one non-proportional to lambda_t)   
    for (unsigned int i = 0; i < mcbCC.size(); i++)
        for (int j = LO; j <= ordDF1; j++)
            for (int k = LO; k <= j; k++)
                coeffbnlep00CC.setCoeff(*coeffbnlep00CC.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC[i].getMu(), orders(k), NULL_ew, mcbCC[i].getScheme()) *
                    (*(mcbCC[i].getCoeff(orders(j - k)))), orders(j)); 
        
    coeffbnlep00qcd.setScheme(scheme);
    coeffbnlep00CC.setScheme(scheme);
    
    //Puts all together in a wilson coefficient object of 12 component: 
    //the first 10 terms are the ones proportional to lambda_t
    //the last two are the remainder current*current operators 
    gslpp::complex appoggio[10], appoggio1[10], appoggio2[10], appoggio3[10];
    for (int i=0; i<10; i++) {
        appoggio[i] = 0.;
        appoggio1[i] = 0.;
        appoggio2[i] = 0.;
        appoggio3[i] = 0.;
    }
    for (int j=LO; j <= ordDF1; j++) {
        bnlep2 = *coeffbnlep00qcd.getCoeff(orders(j));
        std::cout<< std::endl <<"$$$$$$$$$$$$ "<< j <<" $$$$$$$$$$$"<<std::endl;
                
        for (int i = 0; i < 10; i++){
            bnlep.assign(i, bnlep2(i));
         if(j == LO){ appoggio[i] = bnlep(i);}
         if( j == NLO){ appoggio1[i] = bnlep(i); }   
            std::cout<<"++++++++++ "<< i <<" -> "<< bnlep(i) <<" ++++++++++++"<<std::endl;
        }
        bnlep2 = *coeffbnlep00CC.getCoeff(orders(j));
        for (int i = 10; i < 12; i++){
            bnlep.assign(i, bnlep2(i-10));
            std::cout<<"++++++++++ "<< i <<" -> "<< bnlep(i) <<" ++++++++++++"<<std::endl;
            
        }        
    coeffbnlep00.setCoeff(bnlep, orders(j));
    }
    for (int k=LO_ew; k <= ordDF1ew; k++) {
        bnlep2 = *coeffbnlep00qcd.getCoeff(orders_ew(k));
        std::cout<< std::endl <<"$$$$$$$$$$$$ "<<k<<" $$$$$$$$$$$"<<std::endl;
                
        for (int l = 0; l < 10; l++){
            bnlep.assign(l, bnlep2(l));
            if(k == LO_ew){ appoggio2[l] = bnlep(l); }
            if(k == NLO_ew) { appoggio3[l] = bnlep(l); }
            std::cout<<"++++++++++ "<< l <<" -> "<< bnlep(l) <<" ++++++++++++"<<std::endl;
        }          
    
    coeffbnlep00.setCoeff(bnlep, orders_ew(k));
    }
    
    std::cout<< std::endl << " COEFFICIENTI TOTALI " << std::endl; 
    gslpp::complex aaa = 0. ;
    for (int l = 0; l < 10; l++){  
        aaa = appoggio[l] + appoggio2[l];
    std::cout<< std::endl << " C_0 + C_ALE/ALPHA_S "<< l <<" -> " << aaa  << std::endl;
       aaa = appoggio1[l] + appoggio3[l];
    std::cout<< std::endl <<" C_1 + C_ALE " << l <<" -> " << aaa  << std::endl;    
    std::cout << std::endl << " ************************************************ " << std::endl;  
    } 
    
    return coeffbnlep00.getCoeff();
}

/*******************************************************************************
 * evolution Wilson Coefficienys b-> nonlep.                                   * 
 * Buras base                                                                  *
 * deltaB=1   deltaC=0   deltaS=1                                              *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDF1bnlep::ComputeCoeffBnlep10(double mu, schemes scheme) 
{
    
    const std::vector<WilsonCoefficient>& mcb = model.getMyMatching()->CMbnlep( 1);
    const std::vector<WilsonCoefficient>& mcbCC = model.getMyMatching()->CMbnlepCC( 1);
    
    coeffbnlep10qcd.setMu(mu);
    coeffbnlep10CC.setMu(mu);
    coeffbnlep10.setMu(mu);
    
    orders_ew ordDF1ew = coeffbnlep10.getOrder_ew();
    orders ordDF1 =  coeffbnlep10.getOrder();
    
    for (unsigned int i = 0; i < mcb.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++)   {
                
                //Evolves the LO terms and the ones proportional to alpha_s 
                coeffbnlep10qcd.setCoeff(*coeffbnlep10qcd.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), orders(k), NULL_ew,  mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders(j - k)))), orders(j));
                
                //Evolves terms proportional to alpha_e and alpha_e/aplha_s
                coeffbnlep10qcd.setCoeff(*coeffbnlep10qcd.getCoeff(orders_ew(j+4)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), NNLO, orders_ew(k+4), mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders(j - k)))), orders_ew(j+4));
            }
        }
            
        coeffbnlep10qcd.setCoeff(*coeffbnlep10qcd.getCoeff(orders_ew(NLO_ew)) +
                    u.Df1Evolnlep(mu, mcb[i].getMu(), orders(LO), NULL_ew, mcb[i].getScheme()) *
                    (*(mcb[i].getCoeff(orders(LO_ew)))), orders_ew(NLO_ew));
    }        
    
    //Evolves the current*current part of the hamiltonian (the one non-proportional to lambda_t) 
    for (unsigned int i = 0; i < mcbCC.size(); i++)
        for (int j = LO; j <= ordDF1; j++)
            for (int k = LO; k <= j; k++)
                coeffbnlep10CC.setCoeff(*coeffbnlep10CC.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC[i].getMu(), orders(k), NULL_ew, mcbCC[i].getScheme()) *
                    (*(mcbCC[i].getCoeff(orders(j - k)))), orders(j)); 
        
    coeffbnlep10qcd.setScheme(scheme);
    coeffbnlep10CC.setScheme(scheme);
    
    //Puts all together in a wilson coefficient object of 12 component: 
    //the first 10 terms are the one proportional to lambda_t
    //the last two are the remainder current*current operators 
    
    for (int j=LO; j <= ordDF1; j++) {
        bnlep2 = *coeffbnlep10qcd.getCoeff(orders(j));
        for (int i = 0; i < 10; i++){
            bnlep.assign(i, bnlep2(i));
        }
        bnlep2 = *coeffbnlep10CC.getCoeff(orders(j));
        for (int i = 10; i < 12; i++){
            bnlep.assign(i, bnlep2(i-10));
        }        
    coeffbnlep10.setCoeff(bnlep, orders(j));
    }
    for (int k=LO_ew; k <= ordDF1ew; k++) {
        bnlep2 = *coeffbnlep10qcd.getCoeff(orders_ew(k));
        for (int l = 0; l < 10; l++){
            bnlep.assign(l, bnlep2(l));;
        }       
    coeffbnlep10.setCoeff(bnlep, orders_ew(k));
    }
    
    return coeffbnlep10.getCoeff();
}

/*******************************************************************************
 * evolution Wilson Coefficienys b-> nonlep.                                   * 
 * Buras base                                                                  *
 * deltaB=1   deltaC=1   deltaS=0                                              *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDF1bnlep::ComputeCoeffBnlep01(double mu, schemes scheme) 
{
    
    const std::vector<WilsonCoefficient>& mcbCC1 = model.getMyMatching()->CMbnlepCC( 2);
    const std::vector<WilsonCoefficient>& mcbCC2 = model.getMyMatching()->CMbnlepCC( 3);
    
    coeffbnlep01.setMu(mu);
    coeffbnlep01A.setMu(mu);
    coeffbnlep01B.setMu(mu);
    
    orders ordDF1 = coeffbnlep01A.getOrder();
    
    //evolution of the current*current terms
    for (unsigned int i = 0; i < mcbCC1.size(); i++)
        for (int j = LO; j <= ordDF1; j++)
            for (int k = LO; k <= j; k++)
                coeffbnlep01A.setCoeff(*coeffbnlep01A.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC1[i].getMu(), orders(k), NULL_ew, mcbCC1[i].getScheme()) *
                    (*(mcbCC1[i].getCoeff(orders(j - k)))), orders(j)); 
        
    for (unsigned int i = 0; i < mcbCC2.size(); i++)
        for (int j = LO; j <= ordDF1; j++)
            for (int k = LO; k <= j; k++)
                coeffbnlep01B.setCoeff(*coeffbnlep01B.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC2[i].getMu(), orders(k), NULL_ew, mcbCC2[i].getScheme()) *
                    (*(mcbCC2[i].getCoeff(orders(j - k)))), orders(j)); 
        
    coeffbnlep01A.setScheme(scheme);
    coeffbnlep01B.setScheme(scheme);
    coeffbnlep01.setScheme(scheme);
    
    //Puts all together in a wilson coefficient object of 4 components
    for (int j=LO; j <= ordDF1; j++) {
        bnlep2 = *coeffbnlep01A.getCoeff(orders(j));
        for (int i = 0; i < 2; i++){
            bnlep.assign(i, bnlep2(i));
        }
        bnlep2 = *coeffbnlep01A.getCoeff(orders(j));
        for (int i = 2; i < 4; i++){
            bnlep.assign(i, bnlep2(i-2));
        }        
    coeffbnlep01.setCoeff(bnlep, orders(j));
    }    
    return coeffbnlep01.getCoeff();   
    
}

/*******************************************************************************
 * evolution Wilson Coefficienys b-> nonlep.                                   * 
 * Buras base                                                                  *
 * deltaB=1   deltaC=1   deltaS=1                                              *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDF1bnlep::ComputeCoeffBnlep11(double mu, schemes scheme) 
{
    
    const std::vector<WilsonCoefficient>& mcbCC1 = model.getMyMatching()->CMbnlepCC( 2);
    const std::vector<WilsonCoefficient>& mcbCC2 = model.getMyMatching()->CMbnlepCC( 3);
    
    coeffbnlep11.setMu(mu);
    coeffbnlep11A.setMu(mu);
    coeffbnlep11B.setMu(mu);
    
    orders ordDF1 = coeffbnlep11A.getOrder();
    
    for (unsigned int i = 0; i < mcbCC1.size(); i++)
        for (int j = LO; j <= ordDF1; j++)
            for (int k = LO; k <= j; k++)
                coeffbnlep11A.setCoeff(*coeffbnlep11A.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC1[i].getMu(), orders(k), NULL_ew, mcbCC1[i].getScheme()) *
                    (*(mcbCC1[i].getCoeff(orders(j - k)))), orders(j)); 
        
    for (unsigned int i = 0; i < mcbCC2.size(); i++)
        for (int j = LO; j <= ordDF1; j++)
            for (int k = LO; k <= j; k++)
                coeffbnlep11B.setCoeff(*coeffbnlep11B.getCoeff(orders(j)) +
                    u.Df1Evolnlep(mu, mcbCC2[i].getMu(), orders(k), NULL_ew, mcbCC2[i].getScheme()) *
                    (*(mcbCC2[i].getCoeff(orders(j - k)))), orders(j)); 
        
    coeffbnlep11A.setScheme(scheme);
    coeffbnlep11B.setScheme(scheme);
    coeffbnlep11.setScheme(scheme);
    
   for (int j=LO; j <= ordDF1; j++) {
        bnlep2 = *coeffbnlep11A.getCoeff(orders(j));
        for (int i = 0; i < 2; i++){
            bnlep.assign(i, bnlep2(i));
        }
        bnlep2 = *coeffbnlep11A.getCoeff(orders(j));
        for (int i = 2; i < 4; i++){
            bnlep.assign(i, bnlep2(i-2));
        }        
    coeffbnlep11.setCoeff(bnlep, orders(j));
    }    
    return coeffbnlep11.getCoeff();   
    
}
