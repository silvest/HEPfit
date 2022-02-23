/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HeffDS1.h"
#include "StandardModel.h"
#include "EvolDF1nlep.h"
#include "EvolDB1Mll.h"

HeffDS1::HeffDS1(const StandardModel & SM)
: model(SM),
coeffds1(10, NDR, NLO, NLO_QED11), coeffds1cc(10, NDR, NLO, NLO_QED11),
coeffds1pnunu(1, NDR, NLO, NLO_QED11), coeffds1mumu(1, NDR, NLO),
u(new EvolDF1nlep(10, NDR, NLO, NLO_QED11, SM)), uM(new EvolDB1Mll(13, NDR, NLO, SM)),
DS1ccLO(10, 0.),DS1ccLO_QED(10, 0.),DS1ccNLO(10, 0.),DS1ccNLO_QED(10, 0.)
{
}

HeffDS1::~HeffDS1()
{
}

gslpp::vector<gslpp::complex>** HeffDS1::ComputeCoeffDS1PPv(double mu, schemes scheme)
{

    const std::vector<WilsonCoefficient>& mcb = model.getMatching().CMK();
    
    coeffds1.setMu(mu); //inizializes to zero the coefficients

    orders ordDF1 = coeffds1.getOrder();
    orders_qed ordqedDF1 = coeffds1.getOrder_qed();

    switch (scheme) {
        case NDR:
            for (unsigned int i = 0; i < mcb.size(); i++) {
                if (i == 0) {
                    for (int j = LO; j <= ordDF1; j++) {
                        for (int k = LO; k <= j; k++) {
                            //Evolves the LO terms and the ones proportional to alpha_s 
                            coeffds1.setCoeff(*coeffds1.getCoeff(orders(j)) +
                                    u->Df1Evolnlep(mu, mcb[i].getMu(),
                                    orders(k), NO_QED, mcb[i].getScheme())*
                                    (*(mcb[i].getCoeff(orders(j - k)))), orders(j));
                        }
                    }
                    switch (ordqedDF1) {
                        case NLO_QED11:
                            coeffds1.setCoeff(u->Df1Evolnlep(mu, mcb[i].getMu(),
                                    LO, LO_QED, mcb[i].getScheme()) * (*(mcb[i].getCoeff(NLO))) +
                                    u->Df1Evolnlep(mu, mcb[i].getMu(),
                                    LO, NO_QED, mcb[i].getScheme()) * (*(mcb[i].getCoeff(NLO_QED11))) +
                                    u->Df1Evolnlep(mu, mcb[i].getMu(),
                                    LO, NLO_QED11, mcb[i].getScheme()) * (*(mcb[i].getCoeff(LO))) +
                                    u->Df1Evolnlep(mu, mcb[i].getMu(),
                                    NLO, NO_QED, mcb[i].getScheme()) * (*(mcb[i].getCoeff(LO_QED))), NLO_QED11);
                        case LO_QED:
                            coeffds1.setCoeff(u->Df1Evolnlep(mu, mcb[i].getMu(),
                                    LO, LO_QED, mcb[i].getScheme()) * (*(mcb[i].getCoeff(LO))) +
                                    u->Df1Evolnlep(mu, mcb[i].getMu(),
                                    LO, NO_QED, mcb[i].getScheme()) * (*(mcb[i].getCoeff(LO_QED))), LO_QED);
                            break;
                        default:
                            throw std::runtime_error("Error in EvolDF1nlep::Df1Evolnlep()");
                    }
                }
            }
            return coeffds1.getCoeff();
        default:
            throw "HeffDS1::ComputeCoeffDS1PPv(double mu, schemes scheme): scheme not implemented";
    }
}

gslpp::vector<gslpp::complex>** HeffDS1::ComputeCoeffDS1PPz(double muc, schemes scheme)
{

    const std::vector<WilsonCoefficient>& mcbCC = model.getMatching().CMKCC();
    coeffds1cc.setMu(muc);
    double mu = 0.;
    /*
    if(muc < model.getMuc()){
        mu = model.getMuc();
    }
    else{
        mu = muc;
    }*/
    // here below assuming a theory of 3 flavours even above muc in light of current lattice results
    mu = model.getMuc();
    
    orders ordDF1 = coeffds1cc.getOrder();
    orders_qed ordqedDF1 = coeffds1cc.getOrder_qed();

    switch (scheme) {
        case NDR:
            for (unsigned int i = 0; i < mcbCC.size(); i++) {
                if (i == 0) {
                    for (int j = LO; j <= ordDF1; j++) {
                        for (int k = LO; k <= j; k++) {
                            //Evolves the LO terms and the ones proportional to alpha_s 
                            coeffds1cc.setCoeff(*coeffds1cc.getCoeff(orders(j)) +
                                    u->Df1Evolnlep(mu, mcbCC[i].getMu(),
                                    orders(k), NO_QED, mcbCC[i].getScheme())*
                                    (*(mcbCC[i].getCoeff(orders(j - k)))), orders(j));
                        }
                    }
                    switch (ordqedDF1) {
                        case NLO_QED11:
                            coeffds1cc.setCoeff(u->Df1Evolnlep(mu, mcbCC[i].getMu(),
                                    LO, LO_QED, mcbCC[i].getScheme()) * (*(mcbCC[i].getCoeff(NLO))) +
                                    u->Df1Evolnlep(mu, mcbCC[i].getMu(),
                                    LO, NO_QED, mcbCC[i].getScheme()) * (*(mcbCC[i].getCoeff(NLO_QED11))) +
                                    u->Df1Evolnlep(mu, mcbCC[i].getMu(),
                                    LO, NLO_QED11, mcbCC[i].getScheme()) * (*(mcbCC[i].getCoeff(LO))) +
                                    u->Df1Evolnlep(mu, mcbCC[i].getMu(),
                                    NLO, NO_QED, mcbCC[i].getScheme()) * (*(mcbCC[i].getCoeff(LO_QED))), NLO_QED11);
                        case LO_QED:
                            coeffds1cc.setCoeff(u->Df1Evolnlep(mu, mcbCC[i].getMu(),
                                    LO, LO_QED, mcbCC[i].getScheme()) * (*(mcbCC[i].getCoeff(LO))) +
                                    u->Df1Evolnlep(mu, mcbCC[i].getMu(),
                                    LO, NO_QED, mcbCC[i].getScheme()) * (*(mcbCC[i].getCoeff(LO_QED))), LO_QED);
                            break;
                        default:
                            throw std::runtime_error("Error in HeffDS1::ComputeCoeffDS1PPz()");
                    }
                    DS1ccLO = *coeffds1cc.getCoeff(LO);
                    DS1ccLO_QED = *coeffds1cc.getCoeff(LO_QED);
                    DS1ccNLO = *coeffds1cc.getCoeff(NLO);
                    DS1ccNLO_QED = *coeffds1cc.getCoeff(NLO_QED11);
                    for (int l = 2; l < 10; l++){ 
                        DS1ccLO.assign(l, 0.);
                        DS1ccLO_QED.assign(l, 0.);
                        DS1ccNLO.assign(l, 0.);
                        DS1ccNLO_QED.assign(l, 0.);
                    }
                    if(mu == model.getMuc()) CharmMatch();
                    // if(muc < mu){
                    if(muc != mu){
                        switch (ordDF1) {
                            case NLO:
                                coeffds1cc.setCoeff(u->Df1Evolnlep3flav(muc,mu,LO, NO_QED, mcbCC[i].getScheme()) * DS1ccNLO +
                                        u->Df1Evolnlep(muc,mu,NLO, NO_QED, mcbCC[i].getScheme()) * DS1ccLO,NLO);
                            case LO:
                                coeffds1cc.setCoeff(u->Df1Evolnlep3flav(muc,mu,LO, NO_QED, mcbCC[i].getScheme()) * DS1ccLO,LO);                            
                            break;
                            default:
                                throw std::runtime_error("Error in HeffDS1::ComputeCoeffDS1PPz()");
                        }
                        switch (ordqedDF1) {
                            case NLO_QED11:
                                coeffds1cc.setCoeff(u->Df1Evolnlep3flav(muc,mu,LO, LO_QED, mcbCC[i].getScheme()) * DS1ccNLO +
                                        u->Df1Evolnlep3flav(muc,mu,LO, NO_QED, mcbCC[i].getScheme()) * DS1ccNLO_QED +
                                        u->Df1Evolnlep3flav(muc,mu,LO, NLO_QED11, mcbCC[i].getScheme()) * DS1ccLO +
                                        u->Df1Evolnlep3flav(muc,mu,NLO, NO_QED, mcbCC[i].getScheme()) * DS1ccLO_QED, NLO_QED11);
                            case LO_QED:
                                coeffds1cc.setCoeff(u->Df1Evolnlep3flav(muc,mu,LO, LO_QED, mcbCC[i].getScheme()) * DS1ccLO +
                                        u->Df1Evolnlep3flav(muc,mu,LO, NO_QED, mcbCC[i].getScheme()) * DS1ccLO_QED, LO_QED);
                                break;
                            default:
                                throw std::runtime_error("Error in HeffDS1::ComputeCoeffDS1PPz()");
                        }
                    }
                    else{
                        switch (ordDF1) {
                            case NLO:
                                coeffds1cc.setCoeff(DS1ccNLO,NLO);
                            case LO:
                                coeffds1cc.setCoeff(DS1ccLO,LO);                            
                            break;
                            default:
                                throw std::runtime_error("Error in HeffDS1::ComputeCoeffDS1PPz()");
                        }
                        switch (ordqedDF1) {
                            case NLO_QED11:
                                coeffds1cc.setCoeff(DS1ccNLO_QED, NLO_QED11);
                            case LO_QED:
                                coeffds1cc.setCoeff(DS1ccLO_QED, LO_QED);
                                break;
                            default:
                                throw std::runtime_error("Error in HeffDS1::ComputeCoeffDS1PPz()");
                        }
                    }
                }     
            }
            return coeffds1cc.getCoeff();
        default:
            throw "HeffDS1::ComputeCoeffDS1PPz(double mu, schemes scheme): scheme not implemented";
    }
}


gslpp::vector<gslpp::complex>** HeffDS1::ComputeCoeffDS1pnunu()
{

    const std::vector<WilsonCoefficient>& mcb = model.getMatching().CMkpnn();

    orders ordDF1 = coeffds1pnunu.getOrder();
    orders_qed ordDF1_ew = coeffds1pnunu.getOrder_qed();

    for (unsigned int i = 0; i < mcb.size(); i++) {
        for (int j = LO; j <= ordDF1; j++) {
            coeffds1pnunu.setCoeff(*coeffds1pnunu.getCoeff(orders(j))
                    + *mcb[i].getCoeff(orders(j)), orders(j));
        }
        for (int j = LO_QED; j <= ordDF1_ew; j++) {
            coeffds1pnunu.setCoeff(*coeffds1pnunu.getCoeff(orders(j))
                    + *mcb[i].getCoeff(orders(j)), orders(j));
        }
    }

    return coeffds1pnunu.getCoeff();
}

gslpp::vector<gslpp::complex>** HeffDS1::ComputeCoeffDS1mumu()
{

    const std::vector<WilsonCoefficient>& mcb = model.getMatching().CMkmm();

    orders ordDF1 = coeffds1mumu.getOrder();

    for (unsigned int i = 0; i < mcb.size(); i++) {
        for (int j = LO; j <= ordDF1; j++) {
            coeffds1mumu.setCoeff(*coeffds1mumu.getCoeff(orders(j))
                    + *mcb[i].getCoeff(orders(j)), orders(j));
        }
    }

    return coeffds1mumu.getCoeff();
}

void HeffDS1::CharmMatch()
{
    double mc = model.Mrun(model.getMuc(), model.getQuarks(QCD::CHARM).getMass(), FULLNNLO);
    double alphaSmuC = model.Als(model.getMuc());
    double logmc2OmuC2 = log(mc * mc / model.getMuc() / model.getMuc());

    DS1ccNLO.assign(2, (-alphaSmuC / 24. / M_PI)*(-2. / 3. * (logmc2OmuC2 + 1.) * DS1ccLO(1)));
    DS1ccNLO.assign(3, (alphaSmuC / 8. / M_PI)*(-2. / 3. * (logmc2OmuC2 + 1.) * DS1ccLO(1)));
    DS1ccNLO.assign(4, (-alphaSmuC / 24. / M_PI)*(-2. / 3. * (logmc2OmuC2 + 1.) * DS1ccLO(1)));
    DS1ccNLO.assign(5, (alphaSmuC / 8. / M_PI)*(-2. / 3. * (logmc2OmuC2 + 1.) * DS1ccLO(1)));

    DS1ccNLO_QED.assign(2, (-alphaSmuC / 24. / M_PI)*(-2. / 3. * (logmc2OmuC2 + 1.) * DS1ccLO_QED(1)));
    DS1ccNLO_QED.assign(3, (alphaSmuC / 8. / M_PI)*(-2. / 3. * (logmc2OmuC2 + 1.) * DS1ccLO_QED(1)));
    DS1ccNLO_QED.assign(4, (-alphaSmuC / 24. / M_PI)*(-2. / 3. * (logmc2OmuC2 + 1.) * DS1ccLO_QED(1)));
    DS1ccNLO_QED.assign(5, (alphaSmuC / 8. / M_PI)*(-2. / 3. * (logmc2OmuC2 + 1.) * DS1ccLO_QED(1)));
    DS1ccNLO_QED.assign(6, -model.getAle() / 6. / M_PI * 4. / 9. * (logmc2OmuC2 + 1.)*(3. * DS1ccLO(0) + DS1ccLO(1)));
    DS1ccNLO_QED.assign(8, -model.getAle() / 6. / M_PI * 4. / 9. * (logmc2OmuC2 + 1.)*(3. * DS1ccLO(0) + DS1ccLO(1)));
}
