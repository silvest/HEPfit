/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDF1.h"
#include "gslpp_complex.h"

extern std::map<std::string,unsigned int> blocks_nops;

HeffDF1::HeffDF1(std::string blocks, const StandardModel & SM, orders order, orders_qed order_qed) 
:       model(SM),
        coeff(blocks_nops.at(blocks), NDR, order, order_qed),
        evolDF1(blocks, NDR, SM, order, order_qed)
{

// unnecessary???
//    for (i = 0; i < 6; i++) {
//        WC_cache.push_back(coeff);
//        Vmu_cache.push_back(0.);
//    }

   this->blocks = blocks;
   this->nops = blocks_nops.at(blocks);
   mu_cache = 0.;

}

HeffDF1::~HeffDF1() 
{}

gslpp::vector<gslpp::complex>** HeffDF1::ComputeCoeff(double mu, schemes scheme)
{
    orders ordDF1 = coeff.getOrder();
    orders_qed ordDF1_qed = coeff.getOrder_qed();
    unsigned int i,j,k,l;

    const std::vector<WilsonCoefficient> mc = model.getMatching().CMDF1(blocks, nops);

    if (mu == mu_cache && scheme == scheme_cache)
    {
        int check = 1;
        for (i = 0; i < mc.size(); i++)
        {
            if (mc[i].getMu() == Vmu_cache[i])
            {
                for (j = LO; j <= ordDF1; j++)
                {
                    for (l = 0; l < coeff.getSize(); l++)
                    {
                        check *= ((*(mc[i].getCoeff(orders(j))))(l) == (*(WC_cache[i].getCoeff(orders(j))))(l));
                    }
                }
                for (j = LO_QED; j <= ordDF1_qed; j++)
                {
                    for (l = 0; l < coeff.getSize(); l++)
                    {
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
    
    coeff.setMu(mu); // also reset the coefficients

    for (i = 0; i < mc.size(); i++)
    {
        Vmu_cache[i] = mc[i].getMu();
        for (j = LO; j <= ordDF1; j++)
            for (k = LO; k <= j; k++)
            {
                coeff.setCoeff(*coeff.getCoeff(orders(j)) +
                        evolDF1.DF1Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                        (*(mc[i].getCoeff(orders(j - k)))), orders(j));
            }

        for (j = LO_QED; j <= ordDF1_qed; j++)
            for (k = LO; k <= (j - LO_QED) % 3; k++) // depends on the QED enum ordering !!!
            {
                coeff.setCoeff(*coeff.getCoeff(orders_qed(j)) +
                        evolDF1.DF1Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) * 
                        (*(mc[i].getCoeff(orders_qed(j - k)))) + evolDF1.DF1Evol(mu, mc[i].getMu(), orders_qed(j - k), mc[i].getScheme()) * 
                        (*(mc[i].getCoeff(orders(k)))), orders_qed(j));
            }
        
        for (j = NLO_QED02; j <= ordDF1_qed; j++)        
            for (k = LO_QED; k <= j - 3; k++)  // depends on the QED enum ordering !!!
            {
                coeff.setCoeff(*coeff.getCoeff(orders_qed(j)) +
                        evolDF1.DF1Evol(mu, mc[i].getMu(), orders_qed(k), mc[i].getScheme()) *
                        (*(mc[i].getCoeff(orders_qed(j - k + 5)))), orders_qed(j));
            }

    }
    
    return coeff.getCoeff();
}
