/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDC1.h"

HeffDC1::HeffDC1(const StandardModel & SM, StandardModelMatching & SM_Matching) 
:       model(SM), coeffdc1(10, NDR, NLO), 
        coeffdc1g(10, NDR, NLO), ug(10, NDR, NLO, SM), u(10, NDR, NLO, SM), 
        ckm(3,0.), COEFF_pi(10,0.), COEFF_K(10,0.)
{
  
    double co = - 4. * model.getGF() / sqrt(2);
    ckm = model.getVCKM();
    gslpp::complex lambdab_cu = - ckm(1,2)*ckm(0,2).conjugate(); // - V_cb * V_ub^{*}
    gslpp::complex lambdad_cu = ckm(1,0)*ckm(0,0).conjugate(); // V_cd * V_ud^{*}
    gslpp::complex lambdas_cu = ckm(1,1)*ckm(0,1).conjugate(); // V_cs * V_us^{*}
    COEFF_pi.assign(0,0, co * lambdad_cu); COEFF_pi.assign(1,1, co * lambdad_cu);
    COEFF_pi.assign(2,2, co * lambdab_cu); COEFF_pi.assign(3,3, co * lambdab_cu); 
    COEFF_pi.assign(4,4, co * lambdab_cu); COEFF_pi.assign(5,5, co * lambdab_cu);
    COEFF_pi.assign(6,6, co * lambdab_cu); COEFF_pi.assign(7,7, co * lambdab_cu);
    COEFF_pi.assign(8,8, 1.); COEFF_pi.assign(9,9, 1.);
    COEFF_K.assign(0,0, co * lambdas_cu);  COEFF_K.assign(1,1, co * lambdas_cu); 
    COEFF_K.assign(2,2, co * lambdab_cu); COEFF_K.assign(3,3, co * lambdab_cu); 
    COEFF_K.assign(4,4, co * lambdab_cu); COEFF_K.assign(5,5, co * lambdab_cu);
    COEFF_K.assign(6,6, co * lambdab_cu); COEFF_K.assign(7,7, co * lambdab_cu);
    COEFF_K.assign(8,8, 1.); COEFF_K.assign(9,9, 1.);
    
}

HeffDC1::~HeffDC1() {
}

 /******************************************************************************
 * evolution Wilson Coefficien D-> pi pi, K K                                  * 
 * Misiak basis                                                                *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDC1::ComputeCoeffDC1_pi(double mu, schemes scheme) 
{
    
    const std::vector<WilsonCoefficient>& mc = model.getMyMatching()->CMd1();
    coeffdc1.setMu(mu); 
    orders ordDF1 = coeffdc1.getOrder();
    for (unsigned int i = 0; i < mc.size(); i++){
        if(i != 1){
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){
                coeffdc1.setCoeff(*coeffdc1.getCoeff(orders(j)) +
                    ug.DC1Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
                }
            }
        }
    } 
    
    coeffdc1.setScheme(scheme);
    
    return coeffdc1.getCoeff();
}

gslpp::vector<gslpp::complex>** HeffDC1::ComputeCoeffDC1_K(double mu, schemes scheme) 
{
    
    const std::vector<WilsonCoefficient>& mc = model.getMyMatching()->CMd1();
    coeffdc1.setMu(mu); 
    
    orders ordDF1 = coeffdc1.getOrder();
    for (unsigned int i = 0; i < mc.size(); i++){
        std::cout << " SIZE i " << i << std::endl << std::endl;
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){
                coeffdc1.setCoeff(*coeffdc1.getCoeff(orders(j)) + 
                    ug.DC1Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
                }
            }
        if(i == 0){
            for (int j = LO; j <= ordDF1; j++){
            coeffdc1.setCoeff(COEFF_K * (*coeffdc1.getCoeff(orders(j))), orders(j));
            }
        }
    }
   
    coeffdc1.setScheme(scheme);
   
    return coeffdc1.getCoeff();
}

/*  
gslpp::vector<gslpp::complex>** HeffDC1::ComputeCoeffDC1g(double mu, schemes scheme) 
{
    
    const std::vector<WilsonCoefficient>& mcg = model.getMyMatching()->CMd1();
    gslpp::vector<gslpp::complex> vecdc1g(10,0.);
    coeffdc1g.setMu(mu); 
    
    orders ordDF1 = coeffdc1g.getOrder();
    for (int i = 0; i < mcg.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++)  {
                coeffdc1g.setCoeff(*coeffdc1g.getCoeff(orders(j)) +
                    ug.DC1Evol(mu, mcg[i].getMu(), orders(k), mcg[i].getScheme()) *
                    (*(mcg[i].getCoeff(orders(j - k)))), orders(j));
            }
            //STAMPA PROVE
      vecdc1g = *coeffdc1g.getCoeff(orders(j));

        std::cout<< std::endl <<" CHROMOMAGNETIC  "<< j <<" CHROMOMAGNETIC "<<std::endl; 
        std::cout<<" ++++++++++ "<< "C_8g_eff" <<" --> "<< vecdc1g(7) <<" ++++++++++++ "<<std::endl; 
        std::cout << std::endl << " \\\\\\\\\\\\\\\\\\\\  FINE STAMPE DI PROVA \\\\\\\\\\\\\\\\\\\\\\\\ " <<std::endl <<std::endl;
        }
    }
       
    coeffdc1g.setScheme(scheme);
   
    return coeffdc1g.getCoeff();
}*/

 /*   
  
  gslpp::vector<gslpp::complex> a(6,0.);  
  
  // a(0) = C_3 + C_4/3.;
  // a(1) = C_4 + C_3/3.;
  // a(2) = - 2./3./M_PI * C_8g;
  // a(3) = C_5 + C_6/3.;
  // a(4) = C_6 + C_5/3.;
  // a(5) = - 2./9./M_PI * C_8g;

  vector<double> delta1(10,0.),delta2(10,0.),delta3(10,0.),delta4(10,0.),
                 delta5(10,0.),delta6(10,0.),delta7(10,0.),delta8(10,0.),
                 delta9(10,0.),delta10(10,0.);
   
  delta1(0) = 1.; delta2(1) = 1.; delta3(2) = 1.; delta4(3) = 1.;
  delta5(4) = 1.; delta6(5) = 1.; delta7(6) = 1.; delta8(7) = 1.;
  delta9(8) = 1.; delta10(9) = 1.;
  
  gslpp::vector<gslpp::complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffDC1_K(
            mySM.getBD().getMu(), 
            mySM.getBD().getScheme());
  
  std::cout << " Wilson coefficients SM+SUSY " << *allcoeff[LO] << std::endl; 
  std::cout << " Wilson coefficients SM+SUSY " << *allcoeff[NLO] << std::endl;
  
  a(0) = (*allcoeff[LO] + *allcoeff[NLO]) * delta3 +  (*allcoeff[LO] + *allcoeff[NLO]) 
         * delta4/3.;
  
  a(1) = (*allcoeff[LO] + *allcoeff[NLO]) * delta4 +  (*allcoeff[LO] + *allcoeff[NLO]) 
         * delta3/3.; 
  
  a(2) = - 2./3./M_PI * (*allcoeff[LO] + *allcoeff[NLO]) * delta8 ; 
   
  a(3) =  (*allcoeff[LO] + *allcoeff[NLO]) * delta5 + (*allcoeff[LO] + *allcoeff[NLO]) 
          * delta6/3.;
  
  a(4) =  (*allcoeff[LO] + *allcoeff[NLO]) * delta6 + (*allcoeff[LO] + *allcoeff[NLO]) 
          * delta5/3.; 
  
  a(5) =  - 2./9./M_PI *(*allcoeff[LO] + *allcoeff[NLO]) * delta8 ;
  
  std::cout << " a_3 -->" << a(0);
  std::cout << " a_4 -->" << a(1);
  std::cout << " a_4 cromo --> " << a(2);
  std::cout << " a_5 --> " << a(3);
  std::cout << " a_6 --> " << a(4);
  std::cout << " a_6 cromo --> " << a(5); 
   
   */ 
