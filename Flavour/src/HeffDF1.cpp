/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDF1.h"
#include "gslpp_complex.h"

HeffDF1::HeffDF1(unsigned int nops, std::string blocks, const StandardModel & SM, orders order, orders_qed order_qed) 
:       model(SM),
        coeff(nops, NDR, order, order_qed),
        evolDF1(nops, blocks, NDR, SM, order, order_qed)
{

// unnecessary???
//    for (i = 0; i < 6; i++) {
//        WC_cache.push_back(coeff);
//        Vmu_cache.push_back(0.);
//    }

   this->blocks = blocks;
   this->nops = nops;
   mu_cache = 0.;

}

HeffDF1::~HeffDF1() 
{}

gslpp::vector<gslpp::complex>** HeffDF1::ComputeCoeff_s(double mu, schemes scheme) 
{
    
    coeff.setScheme(scheme);
    orders ordDF1 = coeff.getOrder();
    orders_qed ordDF1_qed = coeff.getOrder_qed();
    
    const std::vector<WilsonCoefficient> mc = model.getMatching().CMDF1(blocks, nops, NDR, ordDF1);

    if (mu == mu_cache && scheme == scheme_cache) {
        int check = 1;
        for (unsigned int i = 0; i < mc.size(); i++) {
            if (mc[i].getMu() == Vmu_cache[i]){
                for (int j = LO; j <= ordDF1; j++) {
                        for (unsigned int l = 0; l < coeff.getSize(); l++) {
                            check *= ((*(mc[i].getCoeff(orders(j))))(l) == (*(WC_cache[i].getCoeff(orders(j))))(l));
                        }
                }
            } else check = 0;
        }
        if (check == 1) return coeff.getCoeff();
    } 
    
    mu_cache = mu;
    scheme_cache = scheme;
    WC_cache.clear();
    WC_cache = mc;
    
    coeff.setMu(mu); // clear the coefficients
    
    for (unsigned int i = 0; i < mc.size(); i++){
        Vmu_cache[i] = mc[i].getMu();
        for (int j = LO; j <= ordDF1; j++){
            for (int k = LO; k <= j; k++){
                coeff.setCoeff(*coeff.getCoeff(orders(j)) +
                    evolDF1.DF1Evol(mu, mc[i].getMu(), orders(k), NO_QED, mc[i].getScheme()) * // TO BE FIXED
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }
        }
    }
    
    return coeff.getCoeff();
}
