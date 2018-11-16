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
   this->blocks = blocks;
   this->nops = blocks_nops.at(blocks);
   mu_cache = 0.;
    // cache initialization
    for (unsigned int i = 0; i < nops; i++) {
        WC_cache.push_back(coeff);
        Vmu_cache.push_back(0.);
    }
}

HeffDF1::~HeffDF1() 
{}

gslpp::vector<gslpp::complex> HeffDF1::LowScaleCoeff(int nm)
{
    double mu = coeff.getMu(), eta, M, alsM, kM, b0, b1, b0e, b1e;
//    orders ordDF1 = coeff.getOrder();
    orders_qed ordDF1_qed = coeff.getOrder_qed();

    gslpp::vector<gslpp::complex> aux01(nops, 0.);
    gslpp::vector<gslpp::complex> aux11(nops, 0.);
    gslpp::vector<gslpp::complex> aux21(nops, 0.);

    if (mu == -1) throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): coeff not initialized.");
    if (model.Nf(mu) != 5) throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): defined for 5 flavours only.");

    M = model.getMuw();
    alsM = model.Als(M, FULLNNNLO, ordDF1_qed == NO_QED ? false : true);
    eta = alsM / model.Als(mu, FULLNNNLO, ordDF1_qed == NO_QED ? false : true);
    alsM /= 4. * M_PI; // AlsM tilde

    switch (nm)
    {
        case 00:
            return (*(coeff.getCoeff(LO)));
        case 10:
            return (*(coeff.getCoeff(NLO)) * eta / alsM);
        case 20:
            return (*(coeff.getCoeff(NNLO)) * eta / alsM * eta / alsM);
    }

    if (ordDF1_qed != NO_QED)
    {
        b0 = model.Beta_s(00, 5.);
        b0e = model.Beta_e(00, 5.);
        b1 = model.Beta_s(10, 5.);
        b1e = model.Beta_e(01, 5.);
        kM = model.Ale(M, FULLNLO) / 4. / M_PI / alsM;
  
        if(blocks == "CPL") {
            aux01.assign(6, (*(coeff.getCoeff(orders_qed(LO_QED))))(6));
            aux11.assign(6, (*(coeff.getCoeff(orders_qed(NLO_QED11))))(6));
            aux21.assign(6, (*(coeff.getCoeff(orders_qed(NLO_QED21))))(6));
            aux01.assign(7, (*(coeff.getCoeff(orders_qed(LO_QED))))(7));
            aux11.assign(7, (*(coeff.getCoeff(orders_qed(NLO_QED11))))(7));
            aux21.assign(7, (*(coeff.getCoeff(orders_qed(NLO_QED21))))(7));
        }
        else if(blocks == "CPML" || blocks == "CPMLQB") {
            aux01.assign(8, (*(coeff.getCoeff(orders_qed(LO_QED))))(8));
            aux11.assign(8, (*(coeff.getCoeff(orders_qed(NLO_QED11))))(8));
            aux21.assign(8, (*(coeff.getCoeff(orders_qed(NLO_QED21))))(8));
            aux01.assign(9, (*(coeff.getCoeff(orders_qed(LO_QED))))(9));
            aux11.assign(9, (*(coeff.getCoeff(orders_qed(NLO_QED11))))(9));
            aux21.assign(9, (*(coeff.getCoeff(orders_qed(NLO_QED21))))(9));
        }
        
        switch (nm)
        {
            case 01:
                return (*(coeff.getCoeff(LO_QED)) / kM / eta);
            case 11:
                return (*(coeff.getCoeff(NLO_QED11)) / kM / alsM);
            case 21:
                return (*(coeff.getCoeff(NLO_QED21)) * eta / kM / alsM / alsM);
            case 02:
                return (*(coeff.getCoeff(NLO_QED02)) / kM / kM / eta / eta + b0e / b0 * (1. - eta) / eta / eta * aux01 / kM);
            case 12:
                return (*(coeff.getCoeff(NLO_QED12)) / alsM / kM / kM / eta + b0e / b0 * (1. - eta) / eta * aux11 / kM / alsM
                        + log(eta) / eta * (b0e * b1 / b0 / b0 - b1e / b0) * aux01 / kM);
            case 22:
                return (*(coeff.getCoeff(NLO_QED22)) / alsM / alsM / kM / kM + b0e / b0 * (1. - eta) * aux21 / kM / alsM / alsM
                        + log(eta) * (b0e * b1 / b0 / b0 - b1e / b0) * aux11 / kM / alsM);
            default:
                throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): undefined order.");
        }
    }
    else
        throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): undefined order.");
}


gslpp::vector<gslpp::complex>** HeffDF1::ComputeCoeff(double mu, schemes scheme)
{
    orders ordDF1 = coeff.getOrder();
    orders_qed ordDF1_qed = coeff.getOrder_qed();
    unsigned int i,j,k,l;

    const std::vector<WilsonCoefficient>& mc = model.getMatching().CMDF1(blocks, nops);

    if (mu == mu_cache && scheme == scheme_cache)
    {
        int check = 1;
        for (i = 0; i < mc.size(); i++)
            if (mc[i].getMu() == Vmu_cache[i])
            {
                for (j = LO; j <= ordDF1; j++)
                    for (l = 0; l < coeff.getSize(); l++)
                        check *= ((*(mc[i].getCoeff(orders(j))))(l) == (*(WC_cache[i].getCoeff(orders(j))))(l));
                for (j = LO_QED; j <= ordDF1_qed; j++)
                    for (l = 0; l < coeff.getSize(); l++)
                        check *= ((*(mc[i].getCoeff(orders(j))))(l) == (*(WC_cache[i].getCoeff(orders(j))))(l));
            } else check = 0;
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
                        (*(mc[i].getCoeff(orders_qed(j - k + LO_QED - 3)))), orders_qed(j));
            }

    }
    
    gslpp::complex c9, c10;
    gslpp::vector<gslpp::complex> zz(nops, 0.);
    int a;
    
    for(j = NLO_QED21; j <= ordDF1_qed; j++)
    {
        a = -1;

        if(blocks == "L") break;
  
        if(blocks == "CPL") {
            c9 = (*(coeff.getCoeff(orders_qed(j))))(6);
            c10 = (*(coeff.getCoeff(orders_qed(j))))(7);
            a = 6;
        }
        else if(blocks == "CPML" || blocks == "CPMLQB") {
            c9 = (*(coeff.getCoeff(orders_qed(j))))(8);
            c10 = (*(coeff.getCoeff(orders_qed(j))))(9);
            a = 8;
        }
        
        coeff.setCoeff(zz, orders_qed(j));

        if(a != -1)
        {
          coeff.setCoeff(a, c9, orders_qed(j));
          coeff.setCoeff(a + 1, c10, orders_qed(j));
        }
        
    }
    
    return coeff.getCoeff();
}
