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
    double mu = coeff.getMu(), eta, alsmu, alsM, kM, kmu, b0, b1, b0e, b1e;
    gslpp::vector<gslpp::complex> test(nops, 0.);
    
    if(mu == -1) throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): coeff not initialized.");
    if(model.Nf(mu) != 5) throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): defined for 5 flavours only.");
    
    alsM = model.Alstilde5(model.getMuw());
    alsmu = model.Alstilde5(mu);
    eta = alsM/ alsmu; // assuming 1 matching scale at muW
    b0 = evolDF1.Beta_s(00, 5.);
    b0e = evolDF1.Beta_e(00, 5.);
    b1 = evolDF1.Beta_e(10, 5.);
    b1e = evolDF1.Beta_e(01, 5.);
    kM = model.getAle() / alsM; // WARNING: CHANGE ME!!!
    kmu = model.getAle() / alsmu; // WARNING: CHANGE ME!!!
    switch (nm)
    {
        case 10:
            return (*(coeff.getCoeff(LO)));
        case 20:
            return ( *(coeff.getCoeff(NLO)) / alsmu );
        case 30:
            return (*(coeff.getCoeff(NNLO)) / alsmu / alsmu);
        case 01:
            std::cout << "test: " << *(coeff.getCoeff(LO_QED)) / kM / eta << std::endl;
            return (*(coeff.getCoeff(LO_QED)) / kmu);
        case 11:
            test = *(coeff.getCoeff(NLO_QED11)) / kM / eta / alsmu;
            return (*(coeff.getCoeff(NLO_QED11)) / kmu / alsmu);
        case 21:
            test = *(coeff.getCoeff(NLO_QED21)) / kM / eta / alsmu / alsmu;
            return (*(coeff.getCoeff(NLO_QED21)) / kmu / alsmu / alsmu);
        case 02:
            std::cout << "test: " << (*(coeff.getCoeff(NLO_QED02)) / kM / kM / eta / eta + b0e / b0 * (1. - eta) / eta / eta * (*(coeff.getCoeff(LO_QED))) / kM) << std::endl;
            return (*(coeff.getCoeff(NLO_QED02)) / kmu / kmu);
        case 12:
            test = (*(coeff.getCoeff(NLO_QED12)) / alsmu / kM / kM / eta / eta + b0e / b0 * (1. - eta) / eta / eta * (*(coeff.getCoeff(NLO_QED11))) / kM / alsmu
                    + log(eta) / eta * (b0e * b1 / b0 / b0 - b1e / b0) * (*(coeff.getCoeff(LO_QED))) / kM);
            return (*(coeff.getCoeff(NLO_QED12)) / alsmu / kmu / kmu);
        case 22:
            std::cout << "test: " << *(coeff.getCoeff(NLO_QED22)) / alsmu / alsmu / kM / kM / eta / eta + b0e / b0 * (1. - eta) / eta / eta * (*(coeff.getCoeff(NLO_QED21))) / kM / alsmu / alsmu
                    + log(eta) / eta * (b0e * b1 / b0 / b0 - b1e / b0) * (*(coeff.getCoeff(NLO_QED11))) / kM / alsmu << std::endl;
            std::cout << "test: " << *(coeff.getCoeff(NLO_QED22)) / alsmu / alsmu / kM / kM / eta / eta  << std::endl;
            return (*(coeff.getCoeff(NLO_QED22)) / alsmu / alsmu / kmu / kmu);
        default:
            throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): undefined order.");
    }
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
    
    return coeff.getCoeff();
}
