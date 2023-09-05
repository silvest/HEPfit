/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Charm_Kpnunu.h"
#include "StandardModel.h"

Charm_Kpnunu::Charm_Kpnunu(const StandardModel& model_i)
: model(model_i),
cp(3, 0.), dcp(3, 0.), c_p(3, 0.), cpmuW0(3, 0.), cpmuW1(3, 0.), cpmuW2(3, 0.),
cb(2, 0.), dcb(2, 0.), c_b(2, 0.), cbmuW0(2, 0.), cbmuW1(2, 0.), cbmuW2(2, 0.),
U4p(3, 0.), U5p(3, 0.), J5p1(3, 0.), J4p1(3, 0.), J5p2(3, 0.), J4p2(3, 0.), dc_p(3, 0.),
U4b(2, 0.), U5b(2, 0.), J5b1(2, 0.), J4b1(2, 0.), J5b2(2, 0.), J4b2(2, 0.), dc_b(2, 0.)
{
   
    CP = 0.;
    CBe = 0.;
    CBt = 0.;
    CP_qed=0.;
    CBe_qed=0.;
    CBt_qed=0.;
    
}

Charm_Kpnunu::~Charm_Kpnunu()
{}

double Charm_Kpnunu::getEtab(){
    return model.Als(model.getMuw()) / model.Als(model.getMub());
}

double Charm_Kpnunu::getEtacb(){
    return model.Als(model.getMub()) / model.Als(model.getMuc());
}
    
double Charm_Kpnunu::getmc_mc(){
    return model.getQuarks(QCD::CHARM).getMass();
}

double Charm_Kpnunu::getEtac(){
    return model.Als(model.getMuc()) / model.Als(model.Mrun(model.getMuc(), model.getQuarks(QCD::CHARM).getMass_scale(),
        model.getQuarks(QCD::CHARM).getMass(), FULLNNLO));
}

double Charm_Kpnunu::getKc(){
    return pow(getEtac(), 24. / 25.);
}

double Charm_Kpnunu::getXi1c(){
    return 15212. / 1875. * (getEtac() - 1.) / getEtac();
}

double Charm_Kpnunu::getXi2c(){
    return 966966391. / 10546875. - 231404944. / 3515625. * (1. / getEtac()) - 272751559. / 10546875. *
        (1. / getEtac())*(1. / getEtac()) - 128. / 5. * (1. - (1. / getEtac())*(1. / getEtac())) * gsl_sf_zeta_int(3);
}

double Charm_Kpnunu::getxc_mc(){
    return getmc_mc() * getmc_mc() / model.Mw() / model.Mw();
}

double Charm_Kpnunu::getxc_mc_qed(){
    return sqrt(2.) * model.sW2_MSbar_Approx() * model.getGF() / M_PI / model.Ale(model.getMz(),FULLNLO) * getmc_mc() * getmc_mc() ;
}

double Charm_Kpnunu::getxc(){
    return getKc() * (1. + model.Als(model.getMuc()) / 4. / M_PI * getXi1c() + (model.Als(model.getMuc()) / 4. / M_PI)*(model.Als(model.getMuc()) / 4. / M_PI) * getXi2c()) * getxc_mc();
}



gslpp::vector<double> Charm_Kpnunu::Cp(orders order) // C_inqcd_penguin
{

    double x = model.getMatching().x_t(model.getMuw());
    double L = log(model.getMuw() * model.getMuw() / model.Mw() / model.Mw());

    switch (order) {
        case(LO):
            cp(0) = 4.;
            cp(1) = 4.;
            cp(2) = 0.;

            return (cp);

        case(NLO):
            cp(0) = 4. / 3. * (11. + 6. * L);
            cp(1) = -8. / 3. * (11. + 6. * L);
            cp(2) = 8. * (2. + L);

            return (cp);

        case(NNLO):
            cp(0) = 4. * (-(135677. - 124095.) / 3600. + 58. / 18. * M_PI * M_PI - 0.5 * (2. / 3.)*
                    (112. / 9. + 32. * x + (20. / 3. + 16. * x) * log(x) - (8. + 16. * x) *
                    sqrt(4. * x - 1.) * gsl_sf_clausen(2. * asin(1. / (2. * sqrt(x)))))
                    +(5. / 36. * 238. * L) + 58. / 6. * L * L);
            cp(1) = 4. * (-(135677. + 124095.) / 3600. - 44. / 18. * M_PI * M_PI + 0.5 * (4. / 3.)*
                    (112. / 9. + 32. * x + (20. / 3. + 16. * x) * log(x) - (8. + 16. * x) *
                    sqrt(4. * x - 1.) * gsl_sf_clausen(2. * asin(1. / (2. * sqrt(x)))))
                    - (5. / 36. * 260. * L) - 44. / 6. * L * L);
            cp(2) = 4. * model.getCF()*(33. + 4. * M_PI * M_PI + 34. * L + 12. * L * L);

            return (cp);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::Cp() order" + out.str() + " not implemented");
    }
}

gslpp::matrix<double> Charm_Kpnunu::RGevolP(int nf, orders order)
{

    gslpp::matrix<double> evo(3, 3, 0.);
    double etab = getEtab();
    double etacb = getEtacb();

    switch (order) {
        case (LO):
            switch (nf) {
                case(5): //U5P(0)
                    evo(0, 0) = pow(etab, 6. / 23.);
                    evo(1, 1) = pow(etab, -12. / 23.);
                    evo(2, 0) = 12. / 5. * (pow(etab, 6. / 23.) - pow(etab, 1. / 23.));
                    evo(2, 1) = 6. / 13. * (pow(etab, -12. / 23.) - pow(etab, 1. / 23.));
                    evo(2, 2) = pow(etab, 1. / 23.);

                    return (evo);

                case(4): //U4P(0)
                    evo(0, 0) = pow(etacb, 6. / 25.);
                    evo(1, 1) = pow(etacb, -12. / 25.);
                    evo(2, 0) = 12. / 7. * (pow(etacb, 6. / 25.) - pow(etacb, -1. / 25.));
                    evo(2, 1) = 6. / 11. * (pow(etacb, -12. / 25.) - pow(etacb, -1. / 25.));
                    evo(2, 2) = pow(etacb, -1. / 25.);

                    return (evo);

                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolP() " + out.str() + " wrong number of falvours");
            }
        case(NLO):
            switch (nf) {
                case(5): //J5P(1)
                    evo(0, 0) = 5165. / 3174.;
                    evo(1, 1) = -2267. / 1587.;
                    evo(2, 0) = -15857. / 1587.;
                    evo(2, 1) = 15305. / 3174.;
                    evo(2, 2) = -14924. / 1587.;

                    return (evo);

                case(4): //J4P(1)
                    evo(0, 0) = 6719. / 3750.;
                    evo(1, 1) = -3569. / 1875.;
                    evo(2, 0) = -15931. / 1875.;
                    evo(2, 1) = 5427. / 1250.;
                    evo(2, 2) = -15212. / 1875.;

                    return (evo);

                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolP() " + out.str() + " wrong number of flavours");
            }
        case(NNLO):
            switch (nf) {
                case(5): //J5P(2)
                    evo(0, 0) = -7.35665;
                    evo(1, 1) = -54.9107;
                    evo(2, 0) = 17.7699;
                    evo(2, 1) = -1.7514;
                    evo(2, 2) = 18.3025;

                    return (evo);

                case(4): //J4P(2)
                    evo(0, 0) = -10.2451;
                    evo(1, 1) = -50.3422;
                    evo(2, 0) = 8.0325;
                    evo(2, 1) = -0.3657;
                    evo(2, 2) = 4.91177;

                    return (evo);

                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolP() " + out.str() + " wrong number of flavours");
            }
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::RgevolP() " + out.str() + " wrong order assignment ");
    }
}

gslpp::vector<double> Charm_Kpnunu::ThresholdCp(orders order)
{

    double mub = model.getMub();
    double Mb = model.Mrun(model.getMub(), model.getQuarks(QCD::BOTTOM).getMass_scale(),
            model.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double L = log(mub * mub / Mb / Mb);
    double etab = getEtab();

    switch (order) {
        case(LO):
            dcp = 0.;

            return (dcp);

        case(NLO):
            dcp = 0.;

            return (dcp);

        case(NNLO): // maybe dcp(0) and dcp(1) need to be multiplied by 4 (not excplitly said in the article but the vector Cp follow this pattern)
            dcp(0) = -pow(etab, 6. / 23.)*(2. / 3. * L * ((631. + 9699.) / 6348. * (1 - etab) + etab * Cp(NLO)(0))
                    - 2. / 3. * (59. / 36. + 1. / 3. * L + L * L));
            dcp(1) = -pow(etab, -12. / 23.)*(2. / 3. * L * ((631. - 9699.) / 6348. * (1 - etab) + etab * Cp(NLO)(1))
                    + 4. / 3. * (59. / 36. + 1. / 3. * L + L * L));
            dcp(2) = (-2. / 3.) * L * ((284704. / 2645. + 694522. / 20631. * etab) * pow(etab, 1. / 23.)
                    -(1033492. / 7935. + 8264. / 529. * etab) * pow(etab, 6. / 23.) + (3058. / 1587. + 18136. / 6877. * etab)
                    * pow(etab, -12. / 23.) + etab * (pow(etab, 1. / 23.) * Cp(NLO)(2) + 48. / 5. * (pow(etab, 6. / 23.)
                    - pow(etab, 1. / 23.)) * Cp(NLO)(0) + 24. / 13. * (pow(etab, -12. / 23.) - pow(etab, 1. / 23.)) * Cp(NLO)(1)));

            return (dcp);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::MatchingCp() " + out.str() + " wrong order assignment ");
    }
}

gslpp::vector<double> Charm_Kpnunu::C_p(orders order)
{

    cpmuW0 = Cp(LO);
    cpmuW1 = Cp(NLO);
    cpmuW2 = Cp(NNLO);

    U4p = RGevolP(4, LO);
    U5p = RGevolP(5, LO);
    J4p1 = RGevolP(4, NLO);
    J5p1 = RGevolP(5, NLO);
    J4p2 = RGevolP(4, NNLO);
    J5p2 = RGevolP(5, NNLO);
    
    double etab= getEtab();
    double etacb= getEtacb();

    for (int i = 0; i < 3; i++) {
        dc_p(i, i) = ThresholdCp(NNLO)(i);
    }

    switch (order) {
        case(LO):
            c_p = U4p * U5p*cpmuW0;

            return (c_p);

        case(NLO):
            c_p = J4p1 * U4p * U5p * cpmuW0 + etacb * U4p * (J5p1 - J4p1) * U5p * cpmuW0 +
                    etab * etacb * U4p * U5p * (cpmuW1 - J5p1 * cpmuW0);

            return (c_p);

        case(NNLO):
            //threshold term need to be checked (not really clear in (86) of hep-ph/0603079v2!)
            c_p = J4p2 * U4p * U5p * cpmuW0 + etacb * J4p1 * U4p * (J5p1 - J4p1) * U5p * cpmuW0 +
                    etab * etacb * J4p1 * U4p * U5p * (cpmuW1 - J5p1 * cpmuW0) + etacb * etacb * U4p *
                    (J5p2 - J4p2 - J4p1 * (J5p1 - J4p1) - dc_p) * U5p * cpmuW0 + etab * etacb * etacb *
                    U4p * (J5p1 - J4p1) * U5p * (cpmuW1 - J5p1 * cpmuW0) + etab * etab * etacb * etacb * U4p * U5p *
                    (cpmuW2 - J5p1 * cpmuW1 - (J5p2 - J5p1 * J5p1) * cpmuW0);

            return (c_p);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_p() order" + out.str() + " not implemented");
    }

}

double Charm_Kpnunu::C_P(orders order)
{
    
    double mc_mc = getmc_mc();
    double kc = getKc();
    double xi1c= getXi1c();
    double xi2c= getXi2c();
    double xc_mc=getxc_mc();
    
    double L = log(model.getMuc() * model.getMuc() / mc_mc / mc_mc); // L(mu_c) with mc=mc(mc)
    
    double rhoP1p = 4 * (1. - L) + 4. * log(kc);
    double rhoP1m = -2. * (1. - L) - 2 * log(kc);
    double rhoP2p = 11. + 20. * L - 12. * L * L - 20. * log(kc) - 12. * log(kc) * log(kc) +
            24. * log(kc) * L + 4. * xi1c;
    double rhoP2m = -7. + 12. * L + 12. * L * L - 12. * log(kc) + 12. * log(kc) * log(kc)
            - 24. * log(kc) * L - 2. * xi1c;    
    
    double C_P0 = C(LO,NO_QED,1)(2);
    double C_P1 = C(NLO,NO_QED,1)(2) + 4. * (C(LO,NO_QED,1)(0) * rhoP1p + C(LO,NO_QED,1)(1) * rhoP1m) + xi1c * C(LO,NO_QED,1)(2);
    double C_P2 = C(NNLO,NO_QED,1)(2) + 4. * (C(NLO,NO_QED,1)(0) * rhoP1p + C(NLO,NO_QED,1)(1) * rhoP1m + C(LO,NO_QED,1)(0) * rhoP2p
            + C(LO,NO_QED,1)(1) * rhoP2m) + xi1c * (C(NLO,NO_QED,1)(2) + 4. * (C(LO,NO_QED,1)(0) * rhoP1p + C(LO,NO_QED,1)(1)
            * rhoP1m)) + xi2c * C(LO,NO_QED,1)(2);

    switch (order) {
        case(LO):
            CP = kc * xc_mc / 32. * (4. * M_PI / model.Als(model.getMuc()) * C_P0);
            
            return (CP);
            break;

        case(NLO):
            CP = kc * xc_mc / 32. * (4. * M_PI / model.Als(model.getMuc()) * C_P0 + C_P1);
            
            return (CP);
            break;

        case(NNLO):
            CP = kc * xc_mc / 32. * (4. * M_PI / model.Als(model.getMuc()) * C_P0 + C_P1 + model.Als(model.getMuc())
                    / 4. / M_PI * C_P2);
            
            return (CP);
            break;    
            
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_P() order" + out.str() + " not implemented");
    }
}

gslpp::vector<double> Charm_Kpnunu::Cb(orders order) // C_inqcd_box
{

    double L = log(model.getMuw() * model.getMuw() / model.Mw() / model.Mw());

    switch (order) {
        case(LO):
            cb(0) = 4.; // not sure (see eq 98)
            cb(1) = 0.;

            return (cb);

        case(NLO):
            cb(0) = 0.; // not sure (see eq 98)
            cb(1) = -4. * (9. + 4. * L);

            return (cb);

        case(NNLO):
            cb(0) = 0.; // not sure (see eq 98)
            cb(1) = -8. * model.getCF()*(20. + 2. * M_PI * M_PI + 25. * L + 6. * L * L);

            return (cb);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::Cb() order" + out.str() + " not implemented");
    }
}

gslpp::matrix<double> Charm_Kpnunu::RGevolB(int nf, orders order)
{

    gslpp::matrix<double> evo(2, 2, 0.);
    
    double etab= getEtab();
    double etacb= getEtacb();

    switch (order) {
        case (LO):
            switch (nf) {
                case(5):
                    evo(0, 0) = 1.;
                    evo(0, 1) = 0.;
                    evo(1, 0) = 12. * (1. - pow(etab, 1. / 23.));
                    evo(1, 1) = pow(etab, 1. / 23.);

                    return (evo);

                case(4):
                    evo(0, 0) = 1.;
                    evo(0, 1) = 0.;
                    evo(1, 0) = -12. * (1. - pow(etacb, -1. / 25.));
                    evo(1, 1) = pow(etacb, -1. / 25.);

                    return (evo);

                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolB() " + out.str() + " wrong number of falvours");
            }
        case(NLO):
            switch (nf) {
                case(5):
                    evo(0, 0) = 0.;
                    evo(0, 1) = 0.;
                    evo(1, 0) = 2402. / 1587.;
                    evo(1, 1) = -14924. / 1587.;

                    return (evo);

                case(4):
                    evo(0, 0) = 0.;
                    evo(0, 1) = 0.;
                    evo(1, 0) = 581. / 1875.;
                    evo(1, 1) = -15212. / 1875.;

                    return (evo);

                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolB() " + out.str() + " wrong number of flavours");
            }
        case(NNLO):
            switch (nf) {
                case(5):
                    evo(0, 0) = 0.;
                    evo(0, 1) = 0.;
                    evo(1, 0) = 1296371522. / 39457581. - 34624. / 1081. * gsl_sf_zeta_int(3);
                    evo(1, 1) = -177621017. / 7555707. + 800. / 23. * gsl_sf_zeta_int(3);

                    return (evo);

                case(4):
                    evo(0, 0) = 0.;
                    evo(0, 1) = 0.;
                    evo(1, 0) = 684990354. / 19140625. - 6976. / 245. * gsl_sf_zeta_int(3);
                    evo(1, 1) = -272751559. / 10546875. + 128. / 5. * gsl_sf_zeta_int(3);

                    return (evo);

                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolB() " + out.str() + " wrong number of flavours");
            }
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::RgevolB() " + out.str() + " wrong order assignment ");
    }
}

gslpp::vector<double> Charm_Kpnunu::ThresholdCb(orders order)
{

    double mub = model.getMub();
    double Mb = model.Mrun(model.getMub(), model.getQuarks(QCD::BOTTOM).getMass_scale(),
            model.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double L = log(mub * mub / Mb / Mb);
    
    double etab= getEtab();

    switch (order) {
        case(LO):
            dcb = 0.;

            return (dcb);

        case(NLO):
            dcb = 0.;

            return (dcb);

        case(NNLO):
            dcb(0) = 0.;
            dcb(1) = -2. / 3. * L * ((238784. / 529. - 9608. / 1587 * etab) * pow(etab, 1. / 23.) - 1336. / 3. +
                    pow(etab, 24. / 23.) * Cb(NLO)(1));

            return (dcb);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::MatchingCp() " + out.str() + " wrong order assignment ");
    }
}

gslpp::vector<double> Charm_Kpnunu::C_b(orders order)
{

    cbmuW0 = Cb(LO);
    cbmuW1 = Cb(NLO);
    cbmuW2 = Cb(NNLO);

    U4b = RGevolB(4, LO);
    U5b = RGevolB(5, LO);
    J4b1 = RGevolB(4, NLO);
    J5b1 = RGevolB(5, NLO);
    J4b2 = RGevolB(4, NNLO);
    J5b2 = RGevolB(5, NNLO); 
    
    double etab= getEtab();
    double etacb= getEtacb();

    for (int i = 0; i < 2; i++) {    // Not clear, check on the article.
        dc_b(i, i) = ThresholdCb(NNLO)(i);
    }

    switch (order) {
        case(LO):
            c_b = U4b * U5b*cbmuW0;

            return (c_b);

        case(NLO):
            c_b = J4b1 * U4b * U5b * cbmuW0 + etacb * U4b * (J5b1 - J4b1) * U5b * cbmuW0 +
                    etab * etacb * U4b * U5b * (cbmuW1 - J5b1 * cbmuW0);

            return (c_b);

        case(NNLO):
            //threshold term need to be checked (not really clear in (86) of hep-ph/0603079v2!)
            c_b = J4b2 * U4b * U5b * cbmuW0 + etacb * J4b1 * U4b * (J5b1 - J4b1) * U5b * cbmuW0 +
                    etab * etacb * J4b1 * U4b * U5b * (cbmuW1 - J5b1 * cbmuW0) + etacb * etacb * U4b *
                    (J5b2 - J4b2 - J4b1 * (J5b1 - J4b1) - dc_b) * U5b * cbmuW0 + etab * etacb * etacb *
                    U4b * (J5b1 - J4b1) * U5b * (cbmuW1 - J5b1 * cbmuW0) + etab * etab * etacb * etacb * U4b * U5b *
                    (cbmuW2 - J5b1 * cbmuW1 - (J5b2 - J5b1 * J5b1) * cbmuW0);

            return (c_b);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_p() order" + out.str() + " not implemented");
    }

}

double Charm_Kpnunu::C_Be(orders order)
{
    
    double mc_mc = getmc_mc();
    double kc = getKc();
    double xi1c= getXi1c();
    double xi2c= getXi2c();
    double xc_mc= getxc_mc();
    
    double L = log(model.getMuc() * model.getMuc() / mc_mc / mc_mc); // L(mu_c) with mc=mc(mu_c)
    double rhoB1_e = 5. + 4. * L - 4. * log(kc);
    double rhoB2_e = -2. * model.getCF()*(9. - L - 6. * L * L) - 8. / 3. * log(kc) + 16. * log(kc) * log(kc)
            - 32. * log(kc) * L - 4. * xi1c;
    double C_B0e = C(LO,NO_QED,0)(1);
    double C_B1e = C(NLO,NO_QED,0)(1) + 4. * rhoB1_e + xi1c * C(LO,NO_QED,0)(1);
    double C_B2e = C(NNLO,NO_QED,0)(1) + 4. * rhoB2_e + xi1c * C(NLO,NO_QED,0)(1) + 4. * xi1c * rhoB1_e  + xi2c * C(LO,NO_QED,0)(1);

    switch (order) {
        case(LO): 
            CBe = kc * xc_mc / 16. * (4. * M_PI / model.Als(model.getMuc()) * C_B0e);
            
            return (CBe);

        case(NLO): 
            CBe = kc * xc_mc / 16. * (4. * M_PI / model.Als(model.getMuc()) * C_B0e + C_B1e);
            
            return (CBe);

        case(NNLO): 
            CBe = kc * xc_mc / 16. * (4. * M_PI / model.Als(model.getMuc()) * C_B0e + C_B1e + model.Als(model.getMuc()) / 4. / M_PI * C_B2e);
            
            return (CBe);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_Be() order" + out.str() + " not implemented");
    }
}

double Charm_Kpnunu::C_Bt(orders order)
{
    
    
    double mc_mc = getmc_mc();
    double kc = getKc();
    double xi1c= getXi1c();
    double xi2c= getXi2c();
    double xc_mc= getxc_mc();
    
    double L = log(model.getMuc() * model.getMuc() / mc_mc / mc_mc); // L(mu_c) with mc=mc(mc)
    double x_t = model.getLeptons(model.TAU).getMass() * model.getLeptons(model.TAU).getMass() / mc_mc / mc_mc; // 

    double rhoB1_t = 5. + 4. * L + 4. * x_t * log(x_t) / (1. - x_t) + 4. / (x_t - kc)*(kc * log(kc) -
            x_t * (1. - kc) / (1. - x_t) * log(x_t)); 
    double rhoB2_t = -2. * model.getCF()*((9. + 7. * x_t) / (1. - x_t) + x_t * (3. + 13. * x_t) /
            (1. - x_t) / (1. - x_t) * log(x_t) - 12. * x_t / (1. - x_t) * gsl_sf_dilog(1. - x_t) 
            - ((1. - 13. * x_t) / (1. - x_t) - 12. * x_t * x_t / (1. - x_t) / (1. - x_t) * log(x_t))
            * L - 6. * L * L) + 32. / (x_t - kc)*(4. * x_t * (1 - kc) / 3. / (1. - x_t) - x_t
            * ((x_t * (13. - 29. * x_t)) + kc * (3. + 29. * x_t * x_t) - kc * kc * (3. + 13. * x_t))
            / 12. / (x_t - kc) / (1. - x_t) / (1. - x_t) * log(x_t) + kc * (17. * x_t - kc) / 12.
            / (x_t - kc) * log(kc) + x_t * x_t / (x_t - kc) * log(x_t) * log(kc) - (x_t * x_t +
            2. * x_t * kc - kc * kc) / 2. / (x_t - kc) * log(kc) * log(kc) - x_t * (x_t - kc) / 
              (1. - x_t) * gsl_sf_dilog(1. - x_t) - x_t * gsl_sf_dilog(1. - x_t / kc))
            + 32. / (x_t - kc) / (1. - x_t)*(x_t * (1. - kc) + kc * (1. - x_t)*(2. * x_t - kc)
            / (x_t - kc) * log(kc) - x_t * x_t * (1. - kc)*(1. - 2. * x_t + kc) / (x_t - kc) 
            / (1. - x_t) * log(x_t)) * L + 4. * kc / (x_t - kc)*(1. - x_t / (x_t - kc) * log(x_t) 
            + x_t / (x_t - kc) * log(kc)) * xi1c;
    
    
    double C_B0t = C(LO,NO_QED,0)(1);
    double C_B1t = C(NLO,NO_QED,0)(1) + 4. * rhoB1_t + xi1c * C(LO,NO_QED,0)(1);
    double C_B2t = C(NNLO,NO_QED,0)(1) + 4. * rhoB2_t + xi1c * C(NLO,NO_QED,0)(1) + 4. * xi1c * rhoB1_t + xi2c * C(LO,NO_QED,0)(1);

    switch (order) {
        case(LO):
            CBt = kc * xc_mc / 16. * (4. * M_PI / model.Als(model.getMuc()) * C_B0t);

            return (CBt);

        case(NLO):
            CBt = kc * xc_mc / 16. * (4. * M_PI / model.Als(model.getMuc()) * C_B0t + C_B1t);

            return (CBt);

        case(NNLO):
            CBt = kc * xc_mc / 16. * (4. * M_PI / model.Als(model.getMuc()) * C_B0t + C_B1t
                    + model.Als(model.getMuc()) / 4. / M_PI * C_B2t);

            return (CBt);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_Bt() order" + out.str() + " not implemented");
    }
}





double Charm_Kpnunu::C_P_qed(orders_qed order_qed=LO_QED)
{
    double Cnu_0, Cnu_e, Cnu_es, Cnu_1, Cp_1, Cm_1 ;
    double Cp_CA_0 , Cm_CA_0 ,Cp_CA_e , Cm_CA_e ;
    //double Cp_CA_es , Cm_CA_es; //not useful in calculations

    double nf=3.;
    double nfu=1,nfd=2; 
    double etac = getEtac();
    double mc_mc = getmc_mc();
    double kc = getKc();
    double xc_mc_qed=getxc_mc_qed();
    

    Cnu_0  = C(LO,NO_QED,1)(2);
    Cnu_e  = C(LO,NLO_QED11,1)(2);
    Cnu_es = C(LO,NLO_QED21,1)(2);   

    Cm_CA_0  = C(LO,NO_QED,1)(1) / 4.;
    Cm_CA_e  = C(LO,NLO_QED11,1)(1) / 4.;
    //Cm_CA_es = C(LO,NLO_QED21,1)(1) / 4.;  //not useful in calculations
    
    Cp_CA_0  = C(LO,NO_QED,1)(0) / 4.;
    Cp_CA_e  = C(LO,NLO_QED11,1)(0) / 4.;
    //Cp_CA_es = C(LO,NLO_QED21,1)(0) / 4.;  //not useful in calculations  
    
    Cp_1 = C(NLO,NO_QED,1)(0) / 4. ;
    Cm_1 = C(NLO,NO_QED,1)(1) / 4. ;
    Cnu_1 = C(NLO,NO_QED,1)(2);

    double beta0= 11. - 2. / 3. * nf ;
    double beta1= 102. - 38. / 3. * nf ;//taken from hep-ph/0603079v3
    double betaes  = -8. / 9. * (nfu + nfd/4. ) ;

    double gamma_m0= 8. ;
    double gamma_m1= 404. / 3. - 40. / 9. * nf; //taken from hep-ph/0603079v3
    double gamma_me= 8. / 3.;
    double gamma_mes= 32. / 9.;

    double xi_c1= (gamma_m1 / beta0 - gamma_m0 * beta1 / beta0 / beta0) * (1. - 1. / etac) ;
    double xi_ce= gamma_me / beta0 * (etac - 1.);
    double xi_ces= (gamma_mes / beta0 - betaes * gamma_m0 / beta0 / beta0 - beta1 * gamma_me / beta0 / beta0) * 
                   log(etac) + gamma_me / beta0 * (gamma_m0 * beta1 / beta0 / beta0 - gamma_m1 / beta0) * (1. - 1. / etac) * (1. - etac);

    
    double L = log(model.getMuc() * model.getMuc() / mc_mc / mc_mc);
    
    double rhop_e  = 0.;
    double rhom_e  = 0.;
    double rhop_es = 4. * xi_ce;
    double rhom_es = -2. * xi_ce;
    double rhop_1 =  4. * (1. - L) + 4. * log(kc);
    double rhom_1 = -2. * (1. - L) - 2. * log(kc); 


    double C_Pe= Cnu_e + Cnu_0 * xi_ce + 4. * Cp_CA_0 * rhop_e + 4. * Cm_CA_0 * rhom_e;
    double C_Pes= Cnu_es + Cnu_e * xi_c1 + Cnu_1 * xi_ce + Cnu_0 * xi_ces +
                  4. * (rhop_es + rhop_e * xi_c1 + rhop_1 * xi_ce) * Cp_CA_0 + 
                  4. * (rhom_es + rhom_e * xi_c1 + rhom_1 * xi_ce) * Cm_CA_0 +
                  4. * rhop_1 * Cp_CA_e +
                  4. * rhom_1 * Cm_CA_e +
                  4. * rhop_e * Cp_1 + 4. * rhom_e * Cm_1 ;
    switch (order_qed)
    {
    case(LO_QED):
        return (0.);
    
    case(NLO_QED11):
        CP_qed = kc * xc_mc_qed / 32. * (4. * M_PI * model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc()) / model.Als(model.getMuc()) * C_Pe);
        
        return (CP_qed);
    
    case(NLO_QED21):
        CP_qed = kc * xc_mc_qed / 32. * (4. * M_PI * model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc()) / model.Als(model.getMuc()) * C_Pe +  model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc()) * C_Pes);
        
        return (CP_qed);

    default:
        std::stringstream out;
            out << order_qed;
            throw std::runtime_error("Charm_Kpnunu::C_P_qed() order" + out.str() + " not implemented");
    }

}




double Charm_Kpnunu::C_Be_qed(orders_qed order_qed=LO_QED)
{
    double Cnu_0, Cnu_e, CW_0_sq, Cnu_es, Cnu_1, CW_e_sq;

    Cnu_0  = C(LO,NO_QED,0)(1);
    Cnu_e  = C(LO,NLO_QED11,0)(1);
    Cnu_es = C(LO,NLO_QED21,0)(1);   

    CW_0_sq = C(LO,NO_QED,0)(0) / 4.;   
    CW_e_sq = C(LO,NLO_QED11,0)(0) / 4.; 

    Cnu_1 = C(NLO,NO_QED,0)(1); 
    
    

    double nf=3.;
    double nfu=1,nfd=2; 
    
    double etac = getEtac();
    double mc_mc = getmc_mc();
    double kc = getKc();
    double xc_mc_qed=getxc_mc_qed();

    double L = log(model.getMuc() * model.getMuc() / mc_mc / mc_mc);

    double beta0= 11. - 2. / 3. * nf ;
    double beta1= 102. - 38. / 3. * nf ;//taken from hep-ph/0603079v3
    double betaes  = -8. / 9. * (nfu + nfd/4. ) ;

    double gamma_m0= 8. ;
    double gamma_m1= 404. / 3. - 40. / 9. * nf; //taken from hep-ph/0603079v3
    double gamma_me= 8. / 3.;
    double gamma_mes= 32. / 9.;

    double xi_c1= (gamma_m1 / beta0 - gamma_m0 * beta1 / beta0 / beta0) * (1. - 1. / etac) ;
    double xi_ce= gamma_me / beta0 * (etac - 1.);
    double xi_ces= (gamma_mes / beta0 - betaes * gamma_m0 / beta0 / beta0 - beta1 * gamma_me / beta0 / beta0) * 
                   log(etac) + gamma_me / beta0 * (gamma_m0 * beta1 / beta0 / beta0 - gamma_m1 / beta0) * (1. - 1. / etac) * (1. - etac);

    double rho_e = 0.;
    double rho_es = -4. * xi_ce;
    double rho_1 = 5. + 4. * L - 4. * log(kc); 

    double C_Be_e  = Cnu_e + Cnu_0 * xi_ce + 4. * CW_0_sq * rho_e;
    double C_Be_es = Cnu_es + Cnu_e * xi_c1 + Cnu_1 * xi_ce + Cnu_0 * xi_ces +
                     4. * CW_0_sq * rho_es + 4. * CW_0_sq * rho_e * xi_c1 +  
                     8. * CW_e_sq * rho_1 + 4. * CW_0_sq * rho_1 * xi_ce;
                    
    
    switch (order_qed)
    {
    case(LO_QED):
        return (0.);
    
    case(NLO_QED11):
        CBe_qed = kc * xc_mc_qed / 16. * (4. * M_PI * model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc()) / model.Als(model.getMuc()) * C_Be_e);
        return (CBe_qed);
    
    case(NLO_QED21):
        CBe_qed = kc * xc_mc_qed / 16. * (4. * M_PI * model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc()) / model.Als(model.getMuc()) * C_Be_e +  model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc()) * C_Be_es);
        return (CBe_qed);
    

    default:
        std::stringstream out;
            out << order_qed;
            throw std::runtime_error("Charm_Kpnunu::C_Be_qed() order" + out.str() + " not implemented");
    }

}



double Charm_Kpnunu::C_Bt_qed(orders_qed order_qed=LO_QED)
{
    double Cnu_0, Cnu_e, CW_0_sq, Cnu_es, Cnu_1, CW_e_sq;

    Cnu_0  = C(LO,NO_QED,0)(1);
    Cnu_e  = C(LO,NLO_QED11,0)(1);
    Cnu_es = C(LO,NLO_QED21,0)(1);   

    CW_0_sq = C(LO,NO_QED,0)(0) / 4.;   
    CW_e_sq = C(LO,NLO_QED11,0)(0) / 4.; 

    Cnu_1 = C(NLO,NO_QED,0)(1); 

    double nf=3.;
    double nfu=1,nfd=2; 
    
    double mc_mc = getmc_mc();
    double kc = getKc();
    double xc_mc_qed=getxc_mc_qed();
    double etac=getEtac();

    double L = log(model.getMuc() * model.getMuc() / mc_mc / mc_mc);
    double xt = model.getLeptons(model.TAU).getMass() * model.getLeptons(model.TAU).getMass() / mc_mc / mc_mc;

    double beta0= 11. - 2. / 3. * nf ;
    double beta1= 102. - 38. / 3. * nf ;//taken from hep-ph/0603079v3
    double betaes  = -8. / 9. * (nfu + nfd/4. ) ;

    double gamma_m0= 8. ;
    double gamma_m1= 404. / 3. - 40. / 9. * nf; //taken from hep-ph/0603079v3
    double gamma_me= 8. / 3.;
    double gamma_mes= 32. / 9.;

    double xi_c1= (gamma_m1 / beta0 - gamma_m0 * beta1 / beta0 / beta0) * (1. - 1. / etac) ;
    double xi_ce= gamma_me / beta0 * (etac - 1.);
    double xi_ces= (gamma_mes / beta0 - betaes * gamma_m0 / beta0 / beta0 - beta1 * gamma_me / beta0 / beta0) * 
                   log(etac) + gamma_me / beta0 * (gamma_m0 * beta1 / beta0 / beta0 - gamma_m1 / beta0) * (1. - 1. / etac) * (1. - etac);

    double rho_e = 0.;
    double rho_es = -4. * kc * xi_ce * (kc - xt * (1 - log(xt / kc) ) ) / (kc - xt) / (kc - xt) ;
    double rho_1 = 5. + 4. * L + 4. * xt / (1. - xt) + 4. / (xt - kc)*(kc * log(kc)
            - xt * (1. - kc) / (1. - xt) * log(xt));
    

    double C_Bt_e  = Cnu_e + Cnu_0 * xi_ce + 4. * CW_0_sq * rho_e; ;
    double C_Bt_es = Cnu_es + Cnu_e * xi_c1 + Cnu_1 * xi_ce + Cnu_0 * xi_ces +
                     4. * CW_0_sq * rho_es + 4. * CW_0_sq * rho_e * xi_c1 +  
                     8. * CW_e_sq * rho_1 + 4. * CW_0_sq * rho_1 * xi_ce;
                     


    switch (order_qed)
    {
    case(LO_QED):
        return (0.);
    
    case(NLO_QED11):
        CBt_qed = kc * xc_mc_qed / 16. * (4. * M_PI * model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc()) / model.Als(model.getMuc()) * C_Bt_e);
        
        return (CBt_qed);
    
    case(NLO_QED21):
        CBt_qed = kc * xc_mc_qed / 16. * (4. * M_PI * model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc()) / model.Als(model.getMuc()) * C_Bt_e +  model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc()) * C_Bt_es);
        
        return (CBt_qed);
    

    default:
        std::stringstream out;
            out << order_qed;
            throw std::runtime_error("Charm_Kpnunu::C_Bt_qed() order" + out.str() + " not implemented");
    }

}

double Charm_Kpnunu::P_C(orders order, orders_qed order_qed=LO_QED)
{  
    double Xe = C_P(order) + C_Be(order) + C_Be_qed(order_qed) + C_P_qed(order_qed);
    double Xt = C_P(order) + C_Bt(order) + C_Bt_qed(order_qed) + C_P_qed(order_qed);      
    double lambda4 = model.getCKM().getLambda() * model.getCKM().getLambda() * model.getCKM().getLambda() * model.getCKM().getLambda();
    double pc = 1. / lambda4 * (2. / 3. * Xe + 1. / 3. * Xt);
    
    return (pc);
}

double Charm_Kpnunu::C_TOT(orders order, orders_qed order_qed)
{

    /*
    * double xt = model.getMatching().x_t(model.getMut());
    * double Ale = model.getAle(); 
    * double a = 1. / model.getMatching().mt2omh2(model.getMut()); // same thing for CMkpnn in StandardModelMatching.cpp
    */
    gslpp::complex lambdat = model.getCKM().computelamt();
    gslpp::complex lambdac = model.getCKM().computelamc();
    double lambda = model.getCKM().getLambda();
    double lambda5 = model.getCKM().getLambda() * model.getCKM().getLambda() * model.getCKM().getLambda()
            * model.getCKM().getLambda() * model.getCKM().getLambda();
    double IBT = model.getOptionalParameter("DeltaP_cu"); 
    double X = 0.;
    
    
    if ((order == NNLO) && (order_qed == NLO_QED21)) {
        X = lambdat.real() / lambda5 * model.getMatching().Xt(NLO,NLO_QED11); 
        //std::cout << "top (charm): " << X << std::endl;
        //std::cout << "charm (charm): " << (lambdac.real() / lambda)*(P_C(NNLO,NLO_QED21) + IBT) << std::endl;
        return (X + (lambdac.real() / lambda)*(P_C(NNLO,NLO_QED21) + IBT));
    }

    if ((order == NNLO) && (order_qed == NLO_QED11)) {
        X = lambdat.real() / lambda5 * model.getMatching().Xt(NLO,NLO_QED11);
        return (X + (lambdac.real() / lambda)*(P_C(NNLO,NLO_QED11) + IBT));
    }

    if ((order == NLO) && (order_qed == NLO_QED11)) {
        X = lambdat.real() / lambda5 * model.getMatching().Xt(NLO,NLO_QED11);
        return (X + (lambdac.real() / lambda)*(P_C(NLO,NLO_QED11) + IBT));
    }

    if ((order == NLO) && (order_qed == LO_QED)) {
        X = lambdat.real() / lambda5 * model.getMatching().Xt(NLO,LO_QED);
        return (X + (lambdac.real() / lambda)*(P_C(NLO,LO_QED) + IBT));
    }

    if ((order == LO) && (order_qed == LO_QED)) {
        X = lambdat.real() / lambda5 * model.getMatching().Xt(LO,LO_QED) ;
        return (X + (lambdac.real() / lambda)*(P_C(LO,LO_QED) + IBT));
    }
    else {
        std::stringstream out;
        out << order_qed;
        throw std::runtime_error("BRKppnunu::BRKppnunu(): order_qed " + out.str() + "not implemented\n");
        out << order;
        throw std::runtime_error("BRKppnunu::BRKppnunu(): order " + out.str() + "not implemented");
    }
}

gslpp::matrix<double> Charm_Kpnunu::ADM(orders order,orders_qed order_qed, double nf,int contribution){

    double nfu,nfd; 

    if(nf==1){
        nfu=1.;
        nfd=0;
    }else{
        nfu = floor(nf/2.) ;
        nfd = nf-nfu ; 
    }
   

    // QED ADMs are taken from:  0805.4119v1
    // QCD ADMs are taken from:  hep-ph/0603079v3
    // There is a "transpose difference" between the penguin ADM of the 2 articles. For now the QCD article is used for the definition.
    switch (contribution) 
    {
    case (0): // BOX
        {
        switch (order_qed)
        {
        case (NO_QED): // no qed => go to qcd
            {
                switch (order)
                {
                    case (LO): // 0 (same of LO_QED)
                    {
                    double gamma_nuB = -8. ;
                    double gamma_M = 8. ;
                    double beta  = 11. - 2. / 3. * nf ;
                    double gamma_nu = 2. * ( gamma_M - beta ) ;
            
                    gslpp::matrix<double> gamma_T(2, 0.);

                    gamma_T(0,0) = 0.;
                    gamma_T(1,0) = gamma_nuB;
                    gamma_T(1,1) = gamma_nu;

                    return (gamma_T);

                    break;
                    
                    }
                    case (NLO): // 1
                    {
                    double gamma_nuB = 8. * model.getCF() ;
                    double gamma_M = 404. / 3. - 40. / 9. * nf ;
                    double beta  = 102. - 38. / 3. * nf ;
                    double gamma_nu = 2. * ( gamma_M - beta ) ;
            
                    gslpp::matrix<double> gamma_T(2, 0.);

                    gamma_T(0,0) = 0.;
                    gamma_T(1,0) = gamma_nuB;
                    gamma_T(1,1) = gamma_nu;

                    return (gamma_T);

                    break;
                    
                    }
                    case(NNLO): // 2
                    {
                    double gamma_nuB = 2. * model.getCF() * (69. / 3. - 458./3.*3. - (48. / 3. - 96. * 3.) * gsl_sf_zeta_int(3) + 38. / 3. * nf) ;
                    double gamma_M = 2498. - (4432./27. + 320.*gsl_sf_zeta_int(3)/3.) * nf - (280. * nf * nf )/81. ;
                    double beta  = 2857./2. - 5033.*nf/18. + 325.*nf*nf/54. ;
                    double gamma_nu = 2. * ( gamma_M - beta ) ;
            
                    gslpp::matrix<double> gamma_T(2, 0.);

                    gamma_T(0,0) = 0.;
                    gamma_T(1,0) = gamma_nuB;
                    gamma_T(1,1) = gamma_nu;

                    return (gamma_T);

                    break;
                    
                    }
                    default:
                        std::stringstream out;
                        out << order;
                        throw std::runtime_error("Charm_Kpnunu::ADM: order " + out.str() + " not implemented\n");
                    
                }
            }
        case (LO_QED): // 0
            {
            double gamma_W = 0. ;
            double gamma_nuB = -8. ;
            double gamma_M = 8. ;
            double beta  = 11. - 2. / 3. * nf ;
            double gamma_nu = 2. * ( gamma_M - beta ) ;
            
            gslpp::matrix<double> gamma_T(2, 0.);

            gamma_T(0,0) = 2. * gamma_W;
            gamma_T(1,0) = gamma_nuB;
            gamma_T(1,1) = gamma_nu;

            return (gamma_T);

            break;
            }
        case (NLO_QED11): // e
            {
            double gamma_W = -4. ;
            double gamma_nuB = 0. ;
            double gamma_M = 8. / 3. ;
            double beta  = 0. ;
            double gamma_nu = 2. * ( gamma_M - beta )  ;
            
            gslpp::matrix<double> gamma_T(2, 0.);

            gamma_T(0,0) = 2. * gamma_W;
            gamma_T(1,0) = gamma_nuB;
            gamma_T(1,1) = gamma_nu;

            return (gamma_T);
            
            break;
            }
        case (NLO_QED21): // es
            {
            double gamma_W = 4. ;
            double gamma_nuB = -316. / 9. ;
            double gamma_M = 32. / 9. ;
            double beta  = -8. / 9. * (nfu + nfd/4. ) ;
            double gamma_nu = 2. * ( gamma_M - beta )  ;
            
            gslpp::matrix<double> gamma_T(2, 0.);

            gamma_T(0,0) = 2. * gamma_W;
            gamma_T(1,0) = gamma_nuB;
            gamma_T(1,1) = gamma_nu;

            return (gamma_T);

            break;
            }
        default:
            std::stringstream out, out2;
            out << order_qed;
            out2 << contribution;
            throw std::runtime_error("Charm_Kpnunu::ADM: order_qed " + out.str() + " of contribution " + out2.str() + " not implemented\n");
            }
        }
    case (1): // PENGUIN
        {
        switch (order_qed)
        {
        case (NO_QED): // no qed => go to qcd
            {
                switch (order)
                {
                    case (LO): // 0 (is the same of LO_QED)
                    {
                    gslpp::matrix<double> gamma_pm(2,0.);  // from ref 29 of BG08
                    gamma_pm(0,0)= +6.*(1.-1./3.);
                    gamma_pm(0,1)= 0.; 
                    gamma_pm(1,0)= 0.;
                    gamma_pm(1,1)= -6.*(1.+1./3.);

                    double gamma_pnu = -0.5 * ( -4. * ( 1. + 3. ) ) ;
                    double gamma_mnu = -0.5 * ( -4. * ( 1. - 3. ) ) ;
                    double gamma_M = 8. ;
                    double beta  = 11. - 2. / 3. * nf ;
                    double gamma_nu = 2. * ( gamma_M - beta )  ;
            
                    gslpp::matrix<double> gamma_T(3, 0.);

                    gamma_T(0,0) = gamma_pm(0,0);
                    gamma_T(1,0) = gamma_pm(0,1); // inverted beacuse is transposed
                    gamma_T(0,1) = gamma_pm(1,0);
                    gamma_T(1,1) = gamma_pm(1,1);
            
                    gamma_T(2,0) = gamma_pnu;
                    gamma_T(2,1) = gamma_mnu;
            
                    gamma_T(2,2) = gamma_nu;
            
                    return (gamma_T);
                    break;
                    
                    }
                    case (NLO): // 1
                    {
                    gslpp::matrix<double> gamma_pm(2,0.);  // from ref 29 of BG08
                    gamma_pm(0,0)=  (-21./2. + 2.*nf/3.) * (1. - 1./3.);
                    gamma_pm(0,1)= 0.; 
                    gamma_pm(1,0)= 0.;
                    gamma_pm(1,1)=  (-21./2. - 2.*nf/3.) * (1. + 1./3.);

                    double gamma_pnu = -0.5 * ( 16. * (2. - 11.) ) ;
                    double gamma_mnu = -0.5 * ( 16. * (2. + 11.) ) ;
                    double gamma_M = 404. / 3. - 40. / 9. * nf ;
                    double beta  = 102. - 38. / 3. * nf ;
                    double gamma_nu = 2. * ( gamma_M - beta ) ;
            
                    gslpp::matrix<double> gamma_T(3, 0.);

                    gamma_T(0,0) = gamma_pm(0,0);
                    gamma_T(1,0) = gamma_pm(0,1); // inverted beacuse is transposed
                    gamma_T(0,1) = gamma_pm(1,0);
                    gamma_T(1,1) = gamma_pm(1,1);
            
                    gamma_T(2,0) = gamma_pnu;
                    gamma_T(2,1) = gamma_mnu;
            
                    gamma_T(2,2) = gamma_nu;
            
                    return (gamma_T);
                    break;
                    
                    }
                    case(NNLO): // 2
                    {
                    gslpp::matrix<double> gamma_pm(2,0.);  // from ref 29 of BG08
                    gamma_pm(0,0)= 1./300. * (349049. + 201485.) - 1./1350. * (115577. - 9795.)*nf - 130./27. * (1. - 1./3.)*nf*nf - (672. + 80. * (1. - 1./3.)*nf ) * gsl_sf_zeta_int(3);
                    gamma_pm(0,1)= 0.; 
                    gamma_pm(1,0)= 0.;
                    gamma_pm(1,1)= 1./300. * (349049. - 201485.) - 1./1350. * (115577. + 9795.)*nf + 130./27. * (1. + 1./3.)*nf*nf + (672. + 80. * (1. + 1./3.)*nf ) * gsl_sf_zeta_int(3);

                    double gamma_pnu = -0.5 * ( -2./225. * (45124. + 484917.) + 32.*(13. + 15.)*gsl_sf_zeta_int(3) + 144.*nf ) ;
                    double gamma_mnu = -0.5 * ( -2./225. * (45124. - 484917.) + 32.*(13. - 15.)*gsl_sf_zeta_int(3) - 144.*nf );
                    double gamma_M = 2498. - (4432./27. + 320.*gsl_sf_zeta_int(3)/3.) * nf - (280. * nf * nf )/81. ;
                    double beta  = 2857./2. - 5033.*nf/18. + 325.*nf*nf/54. ;
                    double gamma_nu = 2. * ( gamma_M - beta ) ;  ;
            
                    gslpp::matrix<double> gamma_T(3, 0.);

                    gamma_T(0,0) = gamma_pm(0,0);
                    gamma_T(1,0) = gamma_pm(0,1); // inverted beacuse is transposed
                    gamma_T(0,1) = gamma_pm(1,0);
                    gamma_T(1,1) = gamma_pm(1,1);
            
                    gamma_T(2,0) = gamma_pnu;
                    gamma_T(2,1) = gamma_mnu;
            
                    gamma_T(2,2) = gamma_nu;
            
                    return (gamma_T);
                    break;
                    
                    }
                    default:
                        std::stringstream out;
                        out << order;
                        throw std::runtime_error("Charm_Kpnunu::ADM: order " + out.str() + " not implemented\n");
                    
                }
            }
        case (LO_QED): // 0
            {
            gslpp::matrix<double> gamma_pm(2,0.);  // from ref 29 of BG08
            gamma_pm(0,0)= +6.*(1.-1./3.);
            gamma_pm(0,1)= 0.; 
            gamma_pm(1,0)= 0.;
            gamma_pm(1,1)= -6.*(1.+1./3.);

            double gamma_pnu = 2. * 4. ;
            double gamma_mnu = 2. * (-2.) ;
            double gamma_M = 8. ;
            double beta  = 11. - 2. / 3. * nf ;
            double gamma_nu = 2. * ( gamma_M - beta )  ;
            
            gslpp::matrix<double> gamma_T(3, 0.);

            gamma_T(0,0) = gamma_pm(0,0);
            gamma_T(1,0) = gamma_pm(0,1);
            gamma_T(0,1) = gamma_pm(1,0);
            gamma_T(1,1) = gamma_pm(1,1);
            
            gamma_T(2,0) = gamma_pnu;
            gamma_T(2,1) = gamma_mnu;
            
            gamma_T(2,2) = gamma_nu;
            
            return (gamma_T);
            break;
            }
        case (NLO_QED11): // e
            {
            gslpp::matrix<double> gamma_pm(2,0.); // from ref 29 of BG08
            gamma_pm(0,0)= -8./3.;
            gamma_pm(0,1)= 0.; 
            gamma_pm(1,0)= 0.;
            gamma_pm(1,1)= -8./3.;

            double gamma_pnu = 0. ;
            double gamma_mnu = 0. ;
            double gamma_M = 8. / 3. ;
            double beta  = 0. ;
            double gamma_nu = 2. * ( gamma_M - beta )  ;
            
            gslpp::matrix<double> gamma_T(3, 0.);

            gamma_T(0,0) = gamma_pm(0,0);
            gamma_T(1,0) = gamma_pm(0,1);
            gamma_T(0,1) = gamma_pm(1,0);
            gamma_T(1,1) = gamma_pm(1,1);
            
            gamma_T(2,0) = gamma_pnu;
            gamma_T(2,1) = gamma_mnu;
            
            gamma_T(2,2) = gamma_nu;
            
            return (gamma_T);

            break;
            }
        case (NLO_QED21): // es
            {
            gslpp::matrix<double> gamma_12(2,0.); // from ref 29 of BG08 (look out for the different basis)
            gslpp::matrix<complex> trash(2,0.); 
            gslpp::vector<complex> ev_gamma12(2,0.);

            gamma_12.assign(0,0,194./9.);
            gamma_12.assign(0,1,-2./3.); 
            gamma_12.assign(1,0,25./3.);
            gamma_12.assign(1,1,-49./9.) ;

            gamma_12.eigensystem(trash,ev_gamma12);

            double gamma_pnu = 0. ;
            double gamma_mnu = 0. ;
            double gamma_M = 32. / 9. ;
            double beta  = -8. / 9. * (nfu + nfd/4. ) ;
            double gamma_nu = 2. * ( gamma_M - beta )  ;
            
            gslpp::matrix<double> gamma_T(3, 0.);

            gamma_T(0,0) = ev_gamma12(0).real();
            gamma_T(1,1) = ev_gamma12(1).real();
            
            gamma_T(2,0) = gamma_pnu;
            gamma_T(2,1) = gamma_mnu;
            
            gamma_T(2,2) = gamma_nu;
            
            return (gamma_T);

            break;
            }
        default:
            std::stringstream out, out2;
            out << order_qed;
            out2 << contribution;
            throw std::runtime_error("Charm_Kpnunu::ADM: order_qed " + out.str() + " of contribution " + out2.str() + " not implemented\n");
            break;
        }
        break;
        }
    default:
        std::stringstream out ;
        out << contribution;
        throw std::runtime_error("Charm_Kpnunu::ADM: contribution " + out.str() + " not implemented\n");
        break;
    }


}

gslpp::vector<double> Charm_Kpnunu::CWin_muw(orders order,orders_qed order_qed, int contribution)
{

    double mt = model.getMatching().getMt_mut();
    double MH = model.getMHl(); 
    double MW = model.Mw();  
    double MZ = model.getMz(); 
    double L = log(model.getMuw() * model.getMuw() / MZ / MZ); 
    double LW = log(model.getMuw() * model.getMuw() / MW / MW); 
    double x = model.getMatching().x_t(model.getMut());
    
    double mt2 = mt * mt;
    double MW2 = MW * MW;
    double MZ2 = MZ * MZ;
    double MH2 = MH * MH;
    double MH4 = MH2 * MH2;
    double sw2 = model.sW2_MSbar_Approx();
    double sw4 = sw2 * sw2;
    double cw2 = 1. - sw2;


    switch (contribution)
    {
    case(0): // box
        {
        switch (order_qed)
        {
        case(NO_QED): // no qed => go to qcd
        {
            switch (order)
            {
            
                case(LO): // 0 (same of LO_QED)
                {
                gslpp::vector<double> CWin(2, 0.);
                double Cnu= 0. ; 

                CWin(0)= 4.  ;
                CWin(1)= Cnu ;

                return (CWin);
                break;
                
                }
                case (NLO):
                {
                gslpp::vector<double> CWin(2, 0.);
                double Cnu= -4. * (9. + 4. * LW) ; 
                
                CWin(0)= 0.  ;
                CWin(1)= Cnu ;

                return (CWin);
                break;
            
                }
                case (NNLO):
                {
                gslpp::vector<double> CWin(2, 0.);
                double Cnu= -8. * model.getCF()*(20. + 2. * M_PI * M_PI + 25. * LW + 6. * LW * LW) ; 

                CWin(0)= 0.  ;
                CWin(1)= Cnu ;

                return (CWin);
                break;
                
                }
            
                default:
                    std::stringstream out, out2;
                    out << order;
                    out2 << contribution;
                    throw std::runtime_error("Charm_Kpnunu::CWin_muw: order " + out.str() + " of contribution " + out2.str() +" not implemented\n");
                    break;
                    
            }
        }
        case(LO_QED): // 0
            {
            gslpp::vector<double> CWin(2, 0.);
            double CW= 1.;
            double Cnu= 0. ; 

            CWin(0)= 4. * CW * CW ;
            CWin(1)= Cnu ;

            return (CWin);
            break;
            }
        case(NLO_QED11):// e
            {
            gslpp::vector<double> CWin(2, 0.);
            double CW= 0.;
            double Cnu= 0. ; 

            CWin(0)= 4. * CW * CW ;
            CWin(1)= Cnu ;

            return (CWin);
            break;
            }
        case(NLO_QED21):// es
            {
            gslpp::vector<double> CWin(2, 0.);
            double CW= -11. / 3. - 2. * L ;
            double Cnu= 0. ; 

            CWin(0)= 4. * CW * CW ;
            CWin(1)= Cnu ;

            return (CWin);
            break;
            }
        default:
            std::stringstream out, out2;
            out << order_qed;
            out2 << contribution;
            throw std::runtime_error("Charm_Kpnunu::Cwin_qed: order_qed " + out.str() + " of contribution " + out2.str() +" not implemented\n");
            break;
        }
        break;
        }
    
    
    case(1): // penguin
        {
        switch (order_qed)
        {
        case(NO_QED): // no qed => go to qcd
        {
            switch (order)
            {
            
                case(LO):
                {
                gslpp::vector<double> CWin(3, 0.);
                double Cp= 1. ;
                double Cm= 1. ;
                double Cnu= 0. ; 

                CWin(0)= 4. * Cp  ;
                CWin(1)= 4. * Cm ;
                CWin(2)= Cnu ;

                return (CWin);
        
                break;
            
                }
                case (NLO):
                {
                gslpp::vector<double> CWin(3, 0.);
            
                double Cp=  0.5*(1. - 1./3.) * (11. + 6. * LW) ;
                double Cm= -0.5*(1. + 1./3.) * (11. + 6. * LW) ;
                double Cnu= 8. * (2. + LW) ; 

                CWin(0)= 4. * Cp  ;
                CWin(1)= 4. * Cm  ;
                CWin(2)= Cnu ;

                return (CWin);
        
                break;
            
                }
                case (NNLO):
                {
                gslpp::vector<double> CWin(3, 0.);
            
                double Cp= (-(135677. - 124095.) / 3600. + 58. / 18. * M_PI * M_PI - 0.5 * (2. / 3.)*
                    (112. / 9. + 32. * x + (20. / 3. + 16. * x) * log(x) - (8. + 16. * x) *
                    sqrt(4. * x - 1.) * gsl_sf_clausen(2. * asin(1. / (2. * sqrt(x)))))
                    +(5. / 36. * 238. * LW) + 58. / 6. * LW * LW) ;
                double Cm= (-(135677. + 124095.) / 3600. - 44. / 18. * M_PI * M_PI + 0.5 * (4. / 3.)*
                    (112. / 9. + 32. * x + (20. / 3. + 16. * x) * log(x) - (8. + 16. * x) *
                    sqrt(4. * x - 1.) * gsl_sf_clausen(2. * asin(1. / (2. * sqrt(x)))))
                    - (5. / 36. * 260. * LW) - 44. / 6. * LW * LW);
                double Cnu= 4. * model.getCF() * (33. + 4. * M_PI * M_PI + 34. * LW + 12. * LW * LW) ;

                CWin(0)= 4. * Cp ;
                CWin(1)= 4. * Cm ;
                CWin(2)= Cnu ;

                return (CWin);
        
                break;
                
                }
            
                default:
                    std::stringstream out, out2;
                    out << order;
                    out2 << contribution;
                    throw std::runtime_error("Charm_Kpnunu::CWin_muw: order " + out.str() + " of contribution " + out2.str() +" not implemented\n");
                    break;
            
            }
        }

        case(LO_QED): // 0
            {
            gslpp::vector<double> CWin(3, 0.);
            double Cp= 1. ;
            double Cm= 1. ;
            double CA= 1. ;
            double Cnu= 0. ; 

            CWin(0)= 4. * Cp * CA ;
            CWin(1)= 4. * Cm * CA ;
            CWin(2)= Cnu ;

            return (CWin);
        
            break;
            }
        case(NLO_QED11): // e
            {
            gslpp::vector<double> CWin(3, 0.);
            double Cp= 0. ;
            double Cm= 0. ;
            double CA= 0. ;
            double Cnu= 0. ; 

            CWin(0)= 4. * Cp * CA ;
            CWin(1)= 4. * Cm * CA ;
            CWin(2)= Cnu ;

            return (CWin);
        
            break;
            }
        case(NLO_QED21): // es
            {
            gslpp::vector<double> CWin(3, 0.);
            double Cp= -22. / 9. - 4. / 3. * L ;
            double Cm= -22. / 9. - 4. / 3. * L ;
            double CA= 3. * mt2 / 4. / sw2 / MW2 + (11. * sw2 - 6.) / 4. / sw2  / cw2 - 
                       3. / 4. * (MW2 - cw2 * MH2) / (MH2 - MW2) / sw4 * log(MW2 / MZ2) + 
                    3. * MH4 / 4. / (MH2 - MW2) / (MW2 - cw2 * MH2) * log(MH2 / MZ2);
            double Cnu= 0. ; 

            CWin(0)= 4. * Cp * CA ;
            CWin(1)= 4. * Cm * CA ;
            CWin(2)= Cnu ;

            return (CWin);
            break;
            }
        default:
            std::stringstream out, out2;
            out << order_qed;
            out2 << contribution;
            throw std::runtime_error("Charm_Kpnunu::Cwin_qed: order_qed " + out.str() + " of contribution " + out2.str() +" not implemented\n");
            break;
        }
        break;
        }
    default:
        std::stringstream out ;
        out << contribution;
        throw std::runtime_error("Charm_Kpnunu::Cwin_qed: contribution " + out.str() +" not implemented\n");
        break;
        
    }


}


gslpp::matrix<double> Charm_Kpnunu::RGevol_J(orders order, int nf, int contribution)
{
    int dimension;
    
    switch (contribution)
    {
    case(0): // box
        dimension=2;
        break;
    
    case(1): // penguin
        dimension=3;
        break;
    
    default:
        std::stringstream out ;
        out << contribution;
        throw std::runtime_error("Charm_Kpnunu::RGevol_J : contribution " + out.str() +" not implemented\n");
        break;

    }
    
    gslpp::matrix<gslpp::complex> V(dimension,0.)  , V_inv(dimension,0.) ; 
    gslpp::matrix<gslpp::complex> G1(dimension,0.) ;
    gslpp::matrix<gslpp::complex> S1(dimension,0.) ;
    gslpp::matrix<gslpp::complex> J1(dimension,0.) ;
    gslpp::matrix<double> ADM_0(dimension,0.) , ADM_1(dimension,0.) , ADM_2(dimension,0.);
    
    
    gslpp::vector<gslpp::complex> e(dimension,0.) ;
    
    double beta_0  = 11. - 2. / 3. * nf ; 
    double beta_1  = 102. - 38. / 3. * nf ;
    double beta_2 = 2857. / 2. - 5033. / 18. * nf + 325. / 54. * nf * nf  ;
    
    
    ADM_0  = ADM(LO,NO_QED,nf,contribution); // remember that these matrices are all transposed as in the reference
    ADM_1  = ADM(NLO,NO_QED,nf,contribution);
    ADM_2  = ADM(NNLO,NO_QED,nf,contribution);
    
    ADM_0.eigensystem(V,e);
    
    V_inv=V.inverse();

    G1  = V_inv  *  ADM_1  * V ;
    
    
    
    //S1
    for(unsigned int i = 0; i < G1.size_i(); i++){
        for (unsigned int j = 0; j < G1.size_j(); j++) {
            S1.assign(i , j, beta_1 / 2. / beta_0 / beta_0 * e(i).real() * (i==j) - G1(i,j).real() / (2. * beta_0 + e(i).real() - e(j).real() ));
        }
    }

    J1 = V * S1 * V_inv ;
    
    if (order==LO){
        switch (contribution){
            case(0): //box
                return RGevolB(nf,LO);
            case(1): // penguin
                return RGevolP(nf,LO);
            default:
                std::stringstream out ;
                out << contribution;
                throw std::runtime_error("Charm_Kpnunu::RGevol_J : contribution " + out.str() +" not implemented\n");                   
        }
    } else if (order==NLO){
        return J1.real();
    } else if (order==NNLO){
        gslpp::matrix<gslpp::complex> S2(dimension,0.) , G2(dimension,0.) , J2(dimension,0.);
        G2  = V_inv  *  ADM_2  * V ;
        for(unsigned int i = 0; i < G1.size_i(); i++){
            for (unsigned int j = 0; j < G1.size_j(); j++) {
                double term=0. ;
                for (unsigned int k = 0; k < G1.size_j(); k++){
                    term += ( 1. + e(i).real() / 2. / beta_0 - e(k).real() / 2. / beta_0 ) / ( 2. + e(i).real() / 2. / beta_0 - e(j).real() / 2. / beta_0 ) * ( S1(i,k).real() * S1(k,j).real() - beta_1 / beta_0 * S1(i,j).real() * (j==k) );
                }
                S2.assign(i , j , beta_2 / 4. / beta_0 / beta_0 * e(i).real() * (i==j) + term - G2(i,j).real() / (4. * beta_0 + e(i).real() - e(j).real() )); 
            }
        }
        J2 = V * S2 * V_inv ;
        return J2.real();
    
    } else{
        std::stringstream out ;
        out << order;
        throw std::runtime_error("Charm_Kpnunu::RGevol_J : order " + out.str() +" not implemented\n");   
    }
}
                              

gslpp::matrix<double> Charm_Kpnunu::RGevol_R(orders_qed order_qed, int nf, int contribution)
{
    int dimension;
    double etab=getEtab();
    double etacb=getEtacb();
    
    
    
    switch (contribution)
    {
    case(0): // box
        dimension=2;
        break;
    
    case(1): // penguin
        dimension=3;
        break;
    
    default:
        std::stringstream out ;
        out << contribution;
        throw std::runtime_error("Charm_Kpnunu::RGevol_K : contribution " + out.str() +" not implemented\n");
        break;
    }
    
    
    gslpp::matrix<gslpp::complex> V(dimension,0.)  , V_inv(dimension,0.) ; 
    gslpp::matrix<gslpp::complex> M_0(dimension,0.) , K_0(dimension,0.) ; 
    gslpp::matrix<double> ADM_0(dimension,0.) , ADM_1(dimension,0.) , ADM_2(dimension,0.), ADM_e(dimension,0.), ADM_es(dimension,0.);
    gslpp::matrix<gslpp::complex> M_1(dimension,0.), S_1(dimension,0.) , K1_1(dimension,0.) , K2_1(dimension,0.) , K3_1(dimension,0.) ;
                
    
    gslpp::vector<gslpp::complex> e(dimension,0.) ;
    
    double beta_0  = 11. - 2. / 3. * nf ; 
    double beta_1  = 102. - 38. / 3. * nf ;
    //double beta_2 = 2857. / 2. - 5033. / 18. * nf + 325. / 54. * nf * nf  ;
    
    
    ADM_0  = ADM(LO,NO_QED,nf,contribution); // remember that these matrices are all transposed as in the reference
    ADM_e  = ADM(LO,NLO_QED11,nf,contribution); // remember that these matrices are all transposed as in the reference
    ADM_es = ADM(LO,NLO_QED21,nf,contribution); // remember that these matrices are all transposed as in the reference
    ADM_1  = ADM(NLO,NO_QED,nf,contribution);
    ADM_2  = ADM(NNLO,NO_QED,nf,contribution);
    
    ADM_0.eigensystem(V,e);
    V_inv=V.inverse();
    
    M_0 = V_inv * ADM_e * V;
    M_1 = V_inv * ( (ADM_es - beta_1 / beta_0 * ADM_e) + (ADM_e * RGevol_J(NLO,nf,contribution) - RGevol_J(NLO,nf,contribution) * ADM_e)  ) * V;
    S_1 =  V_inv * RGevol_J(NLO,nf,contribution) * V;
    
    switch(nf){
        case(4):
            for(unsigned int i = 0; i < M_0.size_i(); i++){
                for (unsigned int j = 0; j < M_0.size_j(); j++) {
                    if (e(i).real() != e(j).real() + 2. * beta_0){
                        K_0.assign(i , j, M_0(i,j) / (e(i).real() / 2. / beta_0 - e(j).real() / 2. / beta_0 - 1.) * ( pow(etacb,e(j).real() /2. / beta_0) - pow(etacb,e(i).real() /2. / beta_0)  / etacb)  );
                    } else {
                        K_0.assign(i , j, M_0(i,j) * pow(etacb, e(j).real() / 2. / beta_0) * log( 1. / etacb) );
                    }
                }
            }
            break;
        case(5):
            for(unsigned int i = 0; i < M_0.size_i(); i++){
                for (unsigned int j = 0; j < M_0.size_j(); j++) {
                    if (e(i).real() != e(j).real() + 2. * beta_0){
                        K_0.assign(i , j, M_0(i,j) / (e(i).real() / 2. / beta_0 - e(j).real() / 2. / beta_0 - 1.) * ( pow(etab,e(j).real() /2. / beta_0) / etacb - pow(etab,e(i).real() /2. / beta_0)  / etacb / etab) );
                    } else {
                        K_0.assign(i , j, M_0(i,j) * pow(etab, e(j).real() / 2. / beta_0) * log( 1. / etab) / etacb );
                    }
                }
            }
            break;
        default:
            std::stringstream out ;
            out << nf;
            throw std::runtime_error("Charm_Kpnunu::RGevol_K : nf " + out.str() +" not implemented\n");
            break;       
    }
    
    
    switch (nf)
    {
    case(4): //b -> c
        switch(order_qed){
            case(LO_QED):
                return ((-2. * M_PI / beta_0) * (V.real() * K_0.real() * V_inv.real()) );
                break;
            case(NLO_QED11):
                for(unsigned int i = 0; i < M_1.size_i(); i++){
                    for (unsigned int j = 0; j < M_1.size_j(); j++) {
                        if (i!=j){
                            K1_1.assign(i , j, M_1(i,j) / (e(i).real() / 2. / beta_0 - e(j).real() / 2. / beta_0) * ( pow(etacb,e(j).real() /2. / beta_0) - pow(etacb,e(i).real() /2. / beta_0) ) );
                        } else {
                            K1_1.assign(i , j, M_1(i,j) * pow(etacb, e(i).real() / 2. / beta_0) * log( 1. / etacb) );
                        }
                    }
                }
                K2_1 = - etacb * K_0 * S_1;
                K3_1 = S_1 * K_0;
                
                return ( ( -0.5 / beta_0) * V.real() * (K1_1.real() + K2_1.real() + K3_1.real()) * V_inv.real()); 
                break;
            default:
                std::stringstream out ;
                out << order_qed;
                throw std::runtime_error("Charm_Kpnunu::RGevol_K : order " + out.str() +" not implemented\n");
                break;
        }
    case(5): //w -> b
        switch(order_qed){
            case(LO_QED):
                return ((-2. * M_PI / beta_0) * (V.real() * K_0.real() * V_inv.real()) );
                break;
            case(NLO_QED11):
                for(unsigned int i = 0; i < M_1.size_i(); i++){
                    for (unsigned int j = 0; j < M_1.size_j(); j++) {
                        if (i!=j){
                            K1_1.assign(i , j, M_1(i,j) / (e(i).real() / 2. / beta_0 - e(j).real() / 2. / beta_0) * ( pow(etab,e(j).real() /2. / beta_0) - pow(etab,e(i).real() /2. / beta_0) ) );
                        } else {
                            K1_1.assign(i , j, M_1(i,j) * pow(etab, e(i).real() / 2. / beta_0) * log( 1. / etab) );
                        }
                    }
                }
                K2_1 = - etacb * etab * K_0 * S_1;
                K3_1 = etacb * S_1 * K_0;
                
                return ( ( -0.5 / beta_0) * V.real() * (K1_1.real() + K2_1.real() + K3_1.real()) * V_inv.real()); 
                break;
            default:
                std::stringstream out ;
                out << order_qed;
                throw std::runtime_error("Charm_Kpnunu::RGevol_K : order " + out.str() +" not implemented\n");
                break;
        }
    
    default:
        std::stringstream out ;
        out << nf;
        throw std::runtime_error("Charm_Kpnunu::RGevol_K : nf " + out.str() +" not implemented\n");
        break;
    }
    
    
}




gslpp::vector<double> Charm_Kpnunu::C(orders order,orders_qed order_qed, int contribution){
    
    gslpp::vector<double> CW0 =  CWin_muw(LO,NO_QED ,contribution);
    gslpp::vector<double> CW1 =  CWin_muw(NLO,NO_QED ,contribution);
    gslpp::vector<double> CW2 =  CWin_muw(NNLO,NO_QED ,contribution);
    gslpp::vector<double> CWe =  CWin_muw(LO,NLO_QED11 ,contribution);
    gslpp::vector<double> CWes =  CWin_muw(LO,NLO_QED21 ,contribution);
    
    gslpp::matrix<double> U0_4 = RGevol_J(LO,4,contribution);
    gslpp::matrix<double> U0_5 = RGevol_J(LO,5,contribution);
    gslpp::matrix<double> J1_4 = RGevol_J(NLO,4,contribution);
    gslpp::matrix<double> J2_4 = RGevol_J(NNLO,4,contribution);
    gslpp::matrix<double> J1_5 = RGevol_J(NLO,5,contribution);
    gslpp::matrix<double> J2_5 = RGevol_J(NNLO,5,contribution);
    
    gslpp::vector<double> threshold(CW0.size());
    
    double etab = getEtab();
    double etacb = getEtacb();
    
    if (contribution==0){ // box
        threshold = ThresholdCb(NNLO);
    }else if (contribution==1){ // penguin
        threshold = ThresholdCp(NNLO);
    };
    
    switch(order_qed){
        case (NO_QED):{
            switch (order)
            {
                case (LO):{
                    return U0_4 * U0_5 * CW0 ;
                    break;
                }
                case (NLO):{
                    return J1_4 * U0_4 * U0_5 * CW0 + etacb * U0_4 * (J1_5 - J1_4) * U0_5 * CW0 +
                            etab * etacb * U0_4 * U0_5 * (CW1 - J1_5 * CW0) ;

                    break;
                }
                case(NNLO):{ 
                    return J2_4 * U0_4 * U0_5 * CW0 + etacb * J1_4 * U0_4 * (J1_5 - J1_4) * U0_5 * CW0 +
                            etab * etacb * J1_4 * U0_4 * U0_5 * (CW1 - J1_5 * CW0) + etacb * etacb * U0_4 *
                            (J2_5 - J2_4 - J1_4 * (J1_5 - J1_4)) * U0_5 * CW0 + etab * etacb * etacb *
                            U0_4 * (J1_5 - J1_4) * U0_5 * (CW1 - J1_5 * CW0) + etab * etab * etacb * etacb * U0_4 * U0_5 *
                            (CW2 - J1_5 * CW1 - (J2_5 - J1_5 * J1_5) * CW0) - threshold;

                    break;
                }
                default:{
                    std::stringstream out ;
                    out << order;
                    throw std::runtime_error("Charm_Kpnunu::C : order " + out.str() +" not implemented\n");
                    break;
                }
            
            }
        }
        case(NLO_QED11):{
            gslpp::matrix<double> R0_4 = RGevol_R(LO_QED,4,contribution);
            gslpp::matrix<double> R0_5 = RGevol_R(LO_QED,5,contribution);
            
            return ( (U0_4*R0_5/(4.*M_PI) + R0_4*U0_5/(4.*M_PI))*CW0 + U0_4*U0_5*CWe  );
            break;
        }
        case(NLO_QED21):{
            gslpp::matrix<double> R0_4 = RGevol_R(LO_QED,4,contribution);
            gslpp::matrix<double> R0_5 = RGevol_R(LO_QED,5,contribution);
            gslpp::matrix<double> R1_4 = RGevol_R(NLO_QED11,4,contribution);
            gslpp::matrix<double> R1_5 = RGevol_R(NLO_QED11,5,contribution);
            
            return ( (U0_4*R1_5 + (J1_4*U0_4*R0_5)/(4*M_PI) - (U0_4*J1_4*R0_5)/(4*M_PI) + etacb/(4*M_PI)*(R0_4*J1_5*U0_5) -
                    (etacb*etab)/(4*M_PI)*(R0_4*U0_5*J1_5) + R1_4*U0_5)*CW0 + ( (etacb*etab)/(4*M_PI)*U0_4*R0_5 + (etacb*etab)/(4*M_PI)*R0_4*U0_5)*CW1 +
                    (U0_4*U0_5)*CWes + (J1_4*U0_4*U0_5 + etacb*U0_4*J1_5*U0_5 - etacb*U0_4*J1_4*U0_5 - etacb*etab*U0_4*U0_5*J1_5)*CWe  );
            break;   
        }
        default:{
            std::stringstream out ;
            out << order_qed;
            throw std::runtime_error("Charm_Kpnunu::C : order_qed " + out.str() +" not implemented\n");
            break;
        }
            
            
    }
       
}