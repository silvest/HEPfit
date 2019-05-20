/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HeffDF1.h"
#include "gslpp_complex.h"

extern std::map<std::string, unsigned int> blocks_nops;

HeffDF1::HeffDF1(std::string blocks, const StandardModel & SM, qcd_orders order_qcd, qed_orders order_qed)
: model(SM),
coeff(blocks_nops.at(blocks), NDR, order_qcd, order_qed),
evolDF1(blocks, NDR, SM, order_qcd, order_qed)
{
    this->blocks = blocks;
    this->nops = blocks_nops.at(blocks);
    mu_cache = 0.;
}

gslpp::vector<gslpp::complex> HeffDF1::LowScaleCoeff(qcd_orders order_qcd, qed_orders order_qed)
{
    double mu = coeff.getMu(), eta, M, alsM, kM, b0, b1, b0e, b1e;

    //    gslpp::vector<gslpp::complex> aux01(nops, 0.);
    //    gslpp::vector<gslpp::complex> aux11(nops, 0.);
    //    gslpp::vector<gslpp::complex> aux21(nops, 0.);

    if (mu == -1) throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): coeff not initialized.");
    if (model.Nf(mu) != 5) throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): defined for 5 flavours only.");
    if (order_qcd > coeff.getOrder_QCD() || order_qed > coeff.getOrder_QED())
        throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): order not computed at the high scale.");

    M = model.getMuw();
    alsM = model.Als(M, FULLNNNLO, order_qed == QED0 ? false : true);
    eta = alsM / model.Als(mu, FULLNNNLO, order_qed == QED0 ? false : true);
    alsM /= 4. * M_PI; // AlsM tilde
    if (order_qed != QED0)
    {
        b0 = model.Beta_s(00, 5.);
        b0e = model.Beta_e(00, 5.);
        b1 = model.Beta_s(10, 5.);
        b1e = model.Beta_e(01, 5.);
        kM = model.Ale(M, FULLNLO) / 4. / M_PI / alsM;
    }
    //        if(blocks == "CPL") {
    //            aux01.assign(6, (*(coeff.getCoeff(orders_qed(LO_QED))))(6));
    //            aux11.assign(6, (*(coeff.getCoeff(orders_qed(NLO_QED11))))(6));
    //            aux21.assign(6, (*(coeff.getCoeff(orders_qed(NLO_QED21))))(6));
    //            aux01.assign(7, (*(coeff.getCoeff(orders_qed(LO_QED))))(7));
    //            aux11.assign(7, (*(coeff.getCoeff(orders_qed(NLO_QED11))))(7));
    //            aux21.assign(7, (*(coeff.getCoeff(orders_qed(NLO_QED21))))(7));
    //        }
    //        else if(blocks == "CPML" || blocks == "CPMLQB") {
    //            aux01.assign(8, (*(coeff.getCoeff(orders_qed(LO_QED))))(8));
    //            aux11.assign(8, (*(coeff.getCoeff(orders_qed(NLO_QED11))))(8));
    //            aux21.assign(8, (*(coeff.getCoeff(orders_qed(NLO_QED21))))(8));
    //            aux01.assign(9, (*(coeff.getCoeff(orders_qed(LO_QED))))(9));
    //            aux11.assign(9, (*(coeff.getCoeff(orders_qed(NLO_QED11))))(9));
    //            aux21.assign(9, (*(coeff.getCoeff(orders_qed(NLO_QED21))))(9));

    switch (order_qed)
    {
        case QED0:
            switch (order_qcd)
            {
                case QCD0:
                    return coeff.getCoeff(QCD0, QED0);
                case QCD1:
                    return (coeff.getCoeff(QCD1, QED0) * eta / alsM);
                case QCD2:
                    return (coeff.getCoeff(QCD2, QED0) * eta / alsM * eta / alsM);
                default:
                    throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): undefined order.");

            }
        case QED1:
            switch (order_qcd)
            {
                case QCD0:
                    return coeff.getCoeff(QCD0, QED1) / kM / eta;
                case QCD1:
                    return (coeff.getCoeff(QCD1, QED1) / kM / alsM);
                case QCD2:
                    return (coeff.getCoeff(QCD2, QED1) * eta / kM / alsM / alsM);
                default:
                    throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): undefined order.");
            }
        case QED2:
            switch (order_qcd)
            {
                case QCD0:
                    return (coeff.getCoeff(QCD0, QED2) / kM / kM / eta / eta + b0e / b0 * (1. - eta) / eta / eta * coeff.getCoeff(QCD0, QED1) / kM);
                case QCD1:
                    return (coeff.getCoeff(QCD1, QED2) / alsM / kM / kM / eta + b0e / b0 * (1. - eta) / eta * coeff.getCoeff(QCD1, QED1) / kM / alsM
                            + log(eta) / eta * (b0e * b1 / b0 / b0 - b1e / b0) * coeff.getCoeff(QCD0, QED1) / kM);
                case QCD2:
                    return (coeff.getCoeff(QCD2, QED2) / alsM / alsM / kM / kM + b0e / b0 * (1. - eta) * coeff.getCoeff(QCD2, QED1) / kM / alsM / alsM
                            + log(eta) * (b0e * b1 / b0 / b0 - b1e / b0) * coeff.getCoeff(QCD1, QED1) / kM / alsM);
                default:
                    throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): undefined order.");
            }
        default:
            throw std::runtime_error("Error in HeffDF1::LowScaleCoeff(): undefined order.");
    }
}

Expanded<gslpp::vector<gslpp::complex> > HeffDF1::ComputeCoeff(double mu, schemes scheme)
{
    const std::vector<WilsonCoefficientNew>& mc = model.getMatching().CMDF1(blocks, nops);
    uint i;

    if (mu == mu_cache && scheme == scheme_cache)
    {
        bool check = true;
        for (i = 0; i < mc.size(); i++)
            if (mc[i].getMu() == Vmu_cache[i])
                check = check && (mc[i].getCoeff() == WC_cache[i].getCoeff());
            else check = false;
        if (check) return coeff.getCoeff();
    }

    mu_cache = mu;
    scheme_cache = scheme;
    WC_cache = mc;
    Vmu_cache.clear();

    coeff.setMu(mu); // also reset the coefficients

    for (i = 0; i < mc.size(); i++)
    {
        Vmu_cache.push_back(mc[i].getMu());
        coeff.setCoeff(coeff.getCoeff() + evolDF1.DF1Evol(mu, mc[i].getMu(), mc[i].getScheme()) * mc[i].getCoeff()); // multiple matching scales wrong for EW corrections *** TO BE FIXED 
    }

    return coeff.getCoeff();
}
