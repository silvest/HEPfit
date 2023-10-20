/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Charm_Kpnunu.h"
#include "StandardModel.h"

Charm_Kpnunu::Charm_Kpnunu(const StandardModel& model_i)
:   model(model_i),evoCkpnn(1, NDR, NNLO, NLO_QED11), dcp(3, 0.) ,dcb(2, 0.), 
    CW0b(2,0.),CW1b(2,0.),CW2b(2,0.),CWeb(2,0.),CWesb(2,0.),
    U0_4b(2,0.),U0_5b(2,0.),J1_4b(2,0.),J2_4b(2,0.),J1_5b(2,0.),J2_5b(2,0.),
    R0_4b(2,0.),R0_5b(2,0.),R1_4b(2,0.),R1_5b(2,0.),  
    CW0p(3,0.),CW1p(3,0.),CW2p(3,0.),CWep(3,0.),CWesp(3,0.),
    U0_4p(3,0.),U0_5p(3,0.),J1_4p(3,0.),J2_4p(3,0.),J1_5p(3,0.),J2_5p(3,0.),
    R0_4p(3,0.),R0_5p(3,0.),R1_4p(3,0.),R1_5p(3,0.)      
{
    mc_mc=model.getQuarks(QCD::CHARM).getMass();
    etab=model.Als(model.getMuw()) / model.Als(model.getMub());
    etacb=model.Als(model.getMub()) / model.Als(model.getMuc());
    etac=model.Als(model.getMuc()) / model.Als(model.Mrun(model.getMuc(), model.getQuarks(QCD::CHARM).getMass_scale(),
        model.getQuarks(QCD::CHARM).getMass(), FULLNNLO));
    kc= pow(etac, 24. / 25.);
    xc_mc_qed=sqrt(2.) * model.sW2_ND() * model.getGF() / M_PI / model.alphaMz() * mc_mc * mc_mc ;
    L = log(model.getMuc() * model.getMuc() / mc_mc / mc_mc);
    
    xi1c=15212. / 1875. * (etac - 1.) / etac;
    xi2c=966966391. / 10546875. - 231404944. / 3515625. * (1. / etac) - 272751559. / 10546875. *
        (1. / etac)*(1. / etac) - 128. / 5. * (1. - (1. / etac)*(1. / etac)) * gsl_sf_zeta_int(3);;
    xice= 8. / 3. / (11. - 2. / 3. * 3. ) * (etac - 1.);
    xices=(32. / 9. / (11. - 2. / 3. * 3. ) - (-8. / 9. * (1. + 2./4. ) ) * 8. / (11. - 2. / 3. * 3. ) / (11. - 2. / 3. * 3. ) - (102. - 38. / 3. * 3. ) * 8. / 3. / (11. - 2. / 3. * 3. ) / (11. - 2. / 3. * 3. )) * 
                   log(etac) + 8. / 3. / (11. - 2. / 3. * 3. ) * (8. * (102. - 38. / 3. * 3. ) / (11. - 2. / 3. * 3. ) / (11. - 2. / 3. * 3. ) - (404. / 3. - 40. / 9. * 3.) / (11. - 2. / 3. * 3. )) * (1. - 1. / etac) * (1. - etac);

    //box
    CW0b =  CWin_muw(LO,NO_QED ,0);
    CW1b =  CWin_muw(NLO,NO_QED ,0);
    CW2b =  CWin_muw(NNLO,NO_QED ,0);
    CWeb =  CWin_muw(LO,LO_QED ,0);
    CWesb =  CWin_muw(LO,NLO_QED11 ,0);

    U0_4b = RGevol_J(LO,4,0);
    U0_5b = RGevol_J(LO,5,0);
    J1_4b = RGevol_J(NLO,4,0);
    J2_4b = RGevol_J(NNLO,4,0);
    J1_5b = RGevol_J(NLO,5,0);
    J2_5b = RGevol_J(NNLO,5,0);
    
    R0_4b = RGevol_R(LO_QED,4,0);
    R0_5b = RGevol_R(LO_QED,5,0);
    R1_4b = RGevol_R(NLO_QED11,4,0);
    R1_5b = RGevol_R(NLO_QED11,5,0);
     
    //penguin
    CW0p =  CWin_muw(LO,NO_QED ,1);
    CW1p =  CWin_muw(NLO,NO_QED ,1);
    CW2p =  CWin_muw(NNLO,NO_QED ,1);
    CWep =  CWin_muw(LO,LO_QED ,1);
    CWesp =  CWin_muw(LO,NLO_QED11 ,1);

    U0_4p = RGevol_J(LO,4,1);
    U0_5p = RGevol_J(LO,5,1);
    J1_4p = RGevol_J(NLO,4,1);
    J2_4p = RGevol_J(NNLO,4,1);
    J1_5p = RGevol_J(NLO,5,1);
    J2_5p = RGevol_J(NNLO,5,1);
    
    R0_4p = RGevol_R(LO_QED,4,1);
    R0_5p = RGevol_R(LO_QED,5,1);
    R1_4p = RGevol_R(NLO_QED11,4,1);
    R1_5p = RGevol_R(NLO_QED11,5,1); 
}       

Charm_Kpnunu::~Charm_Kpnunu()
{}

gslpp::matrix<double> Charm_Kpnunu::RGevolP(int nf, orders order)
{

    gslpp::matrix<double> evo(3, 3, 0.);    
    
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
    double Lb = log(mub * mub / Mb / Mb);
    
    switch (order) {
        case(LO):
            dcp = 0.;

            return (dcp);

        case(NLO):
            dcp = 0.;

            return (dcp);

        case(NNLO): 
            dcp(0) = -pow(etab, 6. / 23.)*(2. / 3. * Lb * ((631. + 9699.) / 6348. * (1 - etab) + etab * CW1p(0)) 
                    - 2. / 3. * (59. / 36. + 1. / 3. * Lb + Lb * Lb));
            dcp(1) = -pow(etab, -12. / 23.)*(2. / 3. * Lb * ((631. - 9699.) / 6348. * (1 - etab) + etab * CW1p(1))
                    + 4. / 3. * (59. / 36. + 1. / 3. * Lb + Lb * Lb));
            dcp(2) = (-2. / 3.) * Lb * ((284704. / 2645. + 694522. / 20631. * etab) * pow(etab, 1. / 23.)
                    -(1033492. / 7935. + 8264. / 529. * etab) * pow(etab, 6. / 23.) + (3058. / 1587. + 18136. / 6877. * etab)
                    * pow(etab, -12. / 23.) + etab * (pow(etab, 1. / 23.) * CW1p(2) + 48. / 5. * (pow(etab, 6. / 23.)
                    - pow(etab, 1. / 23.)) * CW1p(0) + 24. / 13. * (pow(etab, -12. / 23.) - pow(etab, 1. / 23.)) * CW1p(1)));

            return (dcp);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::MatchingCp() " + out.str() + " wrong order assignment ");
    }
}

double Charm_Kpnunu::C_P(orders order)
{
    switch (order) {
        case(LO):{
            return (C(LO,NO_QED,1)(2));
        }
        case(NLO):{
            double rhoP1p = 4 * (1. - L) + 4. * log(kc);
            double rhoP1m = -2. * (1. - L) - 2 * log(kc);
            gslpp::vector<double> C_LO = C(LO,NO_QED,1);
    
            return (C(NLO,NO_QED,1)(2) + 4. * (C_LO(0) * rhoP1p +C_LO(1) * rhoP1m) + xi1c * C_LO(2));
        }
        case(NNLO):{
            double rhoP1p = 4 * (1. - L) + 4. * log(kc);
            double rhoP1m = -2. * (1. - L) - 2 * log(kc);
            double rhoP2p = 11. + 20. * L - 12. * L * L - 20. * log(kc) - 12. * log(kc) * log(kc) +
                            24. * log(kc) * L + 4. * xi1c;
            double rhoP2m = -7. + 12. * L + 12. * L * L - 12. * log(kc) + 12. * log(kc) * log(kc)
                            - 24. * log(kc) * L - 2. * xi1c;  
            gslpp::vector<double> C_NLO = C(NLO,NO_QED,1);
            gslpp::vector<double> C_LO = C(LO,NO_QED,1);
            
            return (C(NNLO,NO_QED,1)(2) + 4. * (C_NLO(0) * rhoP1p + C_NLO(1) * rhoP1m + C_LO(0) * rhoP2p
                    + C_LO(1) * rhoP2m) + xi1c * (C_NLO(2) + 4. * (C_LO(0) * rhoP1p + C_LO(1)
                    * rhoP1m)) + xi2c * C_LO(2));    
        }
        default:{
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_P() order" + out.str() + " not implemented");
        }
    }
}   

gslpp::matrix<double> Charm_Kpnunu::RGevolB(int nf, orders order)
{

    gslpp::matrix<double> evo(2, 2, 0.);
    
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
    double Lb = log(mub * mub / Mb / Mb);
    
    switch (order) {
        case(LO):
            dcb = 0.;

            return (dcb);

        case(NLO):
            dcb = 0.;

            return (dcb);

        case(NNLO):
            dcb(0) = 0.;
            dcb(1) = -2. / 3. * Lb * ((238784. / 529. - 9608. / 1587 * etab) * pow(etab, 1. / 23.) - 1336. / 3. +
                    pow(etab, 24. / 23.) * CW1b(1));

            return (dcb);

        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::MatchingCp() " + out.str() + " wrong order assignment ");
    }
}

double Charm_Kpnunu::C_Be(orders order)
{
    switch (order) {
        case(LO):{ 
            return (C(LO,NO_QED,0)(1));
        }
        case(NLO): {
            double rhoB1_e = 5. + 4. * L - 4. * log(kc);
    
            return (C(NLO,NO_QED,0)(1) + 4. * rhoB1_e + xi1c * C(LO,NO_QED,0)(1));
        }
        case(NNLO): {
            double rhoB1_e = 5. + 4. * L - 4. * log(kc);
            double rhoB2_e = -2. * model.getCF()*(9. - L - 6. * L * L) - 8. / 3. * log(kc) + 16. * log(kc) * log(kc)
                - 32. * log(kc) * L - 4. * xi1c;

            return (C(NNLO,NO_QED,0)(1) + 4. * rhoB2_e + xi1c * C(NLO,NO_QED,0)(1) + 4. * xi1c * rhoB1_e  + xi2c * C(LO,NO_QED,0)(1));
        }
        default:{
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_Be() order" + out.str() + " not implemented");
        }
    }
}

double Charm_Kpnunu::C_Bt(orders order)
{
    
    switch (order) {
        case(LO):{
            return ( C(LO,NO_QED,0)(1));
        }
        case(NLO):{
            double x_t = model.getLeptons(model.TAU).getMass() * model.getLeptons(model.TAU).getMass() / mc_mc / mc_mc; 
            double rhoB1_t = 5. + 4. * L + 4. * x_t * log(x_t) / (1. - x_t) + 4. / (x_t - kc)*(kc * log(kc) -
                                x_t * (1. - kc) / (1. - x_t) * log(x_t));
            return (C(NLO,NO_QED,0)(1) + 4. * rhoB1_t + xi1c * C(LO,NO_QED,0)(1));
        }
        case(NNLO):{
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
            
            return (C(NNLO,NO_QED,0)(1) + 4. * rhoB2_t + xi1c * C(NLO,NO_QED,0)(1) + 4. * xi1c * rhoB1_t + xi2c * C(LO,NO_QED,0)(1));
        }
        default:{
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_Bt() order" + out.str() + " not implemented");
        }
    }
}





double Charm_Kpnunu::C_P_qed(orders_qed order_qed)
{
 
    switch (order_qed)
    {
    case(NO_QED):{
        return (0.);
    }
    case(LO_QED):{
        double rhop_e  = 0.;
        double rhom_e  = 0.;
      
        double xi_ce= 8. / 3. / (11. - 2. / 3. * 3.) * (etac - 1.);
        
        gslpp::vector<double> C0 = C(LO,NO_QED,1);
        gslpp::vector<double> Ce = C(LO,LO_QED,1);
        
        return (Ce(2) + C0(2) * xi_ce + 4. * C0(0) / 4. * rhop_e + 4. * C0(0) / 4. * rhom_e);
    }
    case(NLO_QED11):{
        double rhop_e  = 0.;
        double rhom_e  = 0.;
        double rhop_es = 4. * xice;
        double rhom_es = -2. * xice;
        double rhop_1 =  4. * (1. - L) + 4. * log(kc);
        double rhom_1 = -2. * (1. - L) - 2. * log(kc); 

        gslpp::vector<double> C0 = C(LO,NO_QED,1);
        gslpp::vector<double> C1 = C(NLO,NO_QED,1);
        gslpp::vector<double> Ce = C(LO,LO_QED,1);
        gslpp::vector<double> Ces = C(LO,NLO_QED11,1);
        
        return (Ces(2) + Ce(2) * xi1c + C1(2) * xice + C0(2) * xices +
                  4. * (rhop_es + rhop_e * xi1c + rhop_1 * xice) * C0(0) / 4. + 
                  4. * (rhom_es + rhom_e * xi1c + rhom_1 * xice) * C0(1) / 4. +
                  4. * rhop_1 * Ce(0) / 4. +
                  4. * rhom_1 * Ce(1) / 4. +
                  4. * rhop_e * C1(0) / 4. + 4. * rhom_e * C1(1) / 4.);
    }
    default:{
        std::stringstream out;
            out << order_qed;
            throw std::runtime_error("Charm_Kpnunu::C_P_qed() order" + out.str() + " not implemented");
        }
    }
}




double Charm_Kpnunu::C_Be_qed(orders_qed order_qed)
{  
    switch (order_qed)
    {
    case(NO_QED):{
        return (0.);
    }
    case(LO_QED):{
        double rho_e = 0.;
        
        gslpp::vector<double> C0 = C(LO,NO_QED,0);
        gslpp::vector<double> Ce = C(LO,LO_QED,0);
    
        return (Ce(1) + C0(1) * xice + 4. * C0(0) / 4. * rho_e);
    }
    case(NLO_QED11):{
        double rho_e = 0.;
        double rho_es = -4. * xice;
        double rho_1 = 5. + 4. * L - 4. * log(kc);
        
        gslpp::vector<double> C0 = C(LO,NO_QED,0);
        gslpp::vector<double> Ce = C(LO,LO_QED,0);
        gslpp::vector<double> C1 = C(NLO,NO_QED,0);
        gslpp::vector<double> Ces = C(LO,NLO_QED11,0);
  
        return (Ces(1) + Ce(1) * xi1c + C1(1) * xice + C0(1) * xices +
                     4. * C0(0) / 4. * rho_es + 4. * C0(0) / 4. * rho_e * xi1c +  
                     8. * Ce(0) / 4. * rho_1 + 4. * C0(0) / 4. * rho_1 * xice);
    }

    default:{
        std::stringstream out;
        out << order_qed;
        throw std::runtime_error("Charm_Kpnunu::C_Be_qed() order" + out.str() + " not implemented");
        }
    }

}

double Charm_Kpnunu::C_Bt_qed(orders_qed order_qed)
{
    
    switch (order_qed)
    {
    case(NO_QED):{
        return (0.);
    }
    case(LO_QED):{
        double rho_e = 0.;
        
        gslpp::vector<double> C0 = C(LO,NO_QED,0);
        gslpp::vector<double> Ce = C(LO,LO_QED,0);
    
        
        return (Ce(1) + C0(1) * xice + 4. * C0(0) / 4. * rho_e);
    }
    case(NLO_QED11):{
        double xt = model.getLeptons(model.TAU).getMass() * model.getLeptons(model.TAU).getMass() / mc_mc / mc_mc;
          
        double rho_e = 0.;
        double rho_es = -4. * kc * xice * (kc - xt * (1 - log(xt / kc) ) ) / (kc - xt) / (kc - xt) ;
        double rho_1 = 5. + 4. * L + 4. * xt / (1. - xt) + 4. / (xt - kc)*(kc * log(kc)
                        - xt * (1. - kc) / (1. - xt) * log(xt));
        
        gslpp::vector<double> C0 = C(LO,NO_QED,0);
        gslpp::vector<double> Ce = C(LO,LO_QED,0);
        gslpp::vector<double> C1 = C(NLO,NO_QED,0);
        gslpp::vector<double> Ces = C(LO,NLO_QED11,0);
  
        return (Ces(1) + Ce(1) * xi1c + C1(1) * xice + C0(1) * xices +
                     4. * C0(0) / 4. * rho_es + 4. * C0(0) / 4. * rho_e * xi1c +  
                     8. * Ce(0) / 4. * rho_1 + 4. * C0(0) / 4. * rho_1 * xice);
    }

    default:{
        std::stringstream out;
        out << order_qed;
        throw std::runtime_error("Charm_Kpnunu::C_Be_qed() order" + out.str() + " not implemented");
        }
    }
}

std::vector<WilsonCoefficient>& Charm_Kpnunu::EVOCkpnn()
{  
    double co = 4. * model.getGF() / sqrt(2.) * model.alphaMz() / 2. / M_PI / model.sW2_ND() ; //SM prefactor as in eq. (1.1) of arXiv:0805.4119
    double cpeng = kc * xc_mc_qed / 32. ;
    double cbox = 2. * cpeng ;
    gslpp::complex lam_c = model.getCKM().computelamc();
    
    vevoCkpnn.clear();
    
    evoCkpnn.setMu(model.getMuc());
    
    switch (evoCkpnn.getOrder()) {
        case NNLO:
            evoCkpnn.setCoeff(0, co * lam_c * model.Als(model.getMuc(),FULLNLO) / 4. / M_PI * (cpeng * C_P(NNLO) + cbox * (2. / 3. * C_Be(NNLO) + 1. / 3. * C_Bt(NNLO))) , NNLO);
        case NLO: 
            evoCkpnn.setCoeff(0, co * lam_c * (cpeng * C_P(NLO) + cbox * (2. / 3. * C_Be(NLO) + 1. / 3. * C_Bt(NLO))) , NLO);
        case LO:
            evoCkpnn.setCoeff(0, co * lam_c * 4. * M_PI / model.Als(model.getMuc(),FULLNLO) * (cpeng * C_P(LO) + cbox * (2. / 3. * C_Be(LO) + 1. / 3. * C_Bt(LO))) , LO);
            break;
        default:
            std::stringstream out;
            out << evoCkpnn.getOrder();
            throw std::runtime_error("Charm_Kpnunu::EVOCkpnn(): order " + out.str() + " not implemented"); 
    }
    
    switch (evoCkpnn.getOrder_qed()) {
        case NLO_QED11:
            evoCkpnn.setCoeff(0, co * lam_c * model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc(),FULLNLO) * (cpeng * C_P_qed(NLO_QED11) + cbox * (2. / 3. * C_Be_qed(NLO_QED11) + 1. / 3. * C_Bt_qed(NLO_QED11))) , NLO_QED11);
        case LO_QED:
            evoCkpnn.setCoeff(0, co * lam_c * 4. * M_PI * model.Ale(model.getMuc(),FULLNLO) / model.Als(model.getMuc(),FULLNLO) / model.Als(model.getMuc(),FULLNLO) * (cpeng * C_P_qed(LO_QED) + cbox * (2. / 3. * C_Be_qed(LO_QED) + 1. / 3. * C_Bt_qed(LO_QED))) , LO_QED);
        //case NO_QED:
        //    evoCkpnn.setCoeff(0, 0. , NO_QED);
            break;
        default:
            std::stringstream out;
            out << evoCkpnn.getOrder_qed();
            throw std::runtime_error("Charm_Kpnunu::EVOCkpnn(): qed order " + out.str() + " not implemented"); 
    }
    
    vevoCkpnn.push_back(evoCkpnn);
    return(vevoCkpnn);
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
                    case (LO): // 0 
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
        case (LO_QED): // e
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
        case (NLO_QED11): // es
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
                    case (LO): // 0 
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
        case (LO_QED): // e
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
        case (NLO_QED11): // es
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
            
                case(LO): // 0 
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
                double MW = model.Mw();  
                double LW = log(model.getMuw() * model.getMuw() / MW / MW); 
                
                gslpp::vector<double> CWin(2, 0.);
                double Cnu= -4. * (9. + 4. * LW) ; 
                
                CWin(0)= 0.  ;
                CWin(1)= Cnu ;

                return (CWin);
                break;
            
                }
                case (NNLO):
                {
                double MW = model.Mw();  
                double LW = log(model.getMuw() * model.getMuw() / MW / MW); 
                
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
        case(LO_QED):// e
            {
            gslpp::vector<double> CWin(2, 0.);
            double CW= 0.;
            double Cnu= 0. ; 

            CWin(0)= 4. * CW * CW ;
            CWin(1)= Cnu ;

            return (CWin);
            break;
            }
        case(NLO_QED11):// es
            {
            double MZ = model.getMz(); 
            double LZ = log(model.getMuw() * model.getMuw() / MZ / MZ); 
             
            gslpp::vector<double> CWin(2, 0.);
            double CW= -11. / 3. - 2. * LZ ;
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
                double MW = model.Mw();  
                double LW = log(model.getMuw() * model.getMuw() / MW / MW); 
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
                double MW = model.Mw();  
                double LW = log(model.getMuw() * model.getMuw() / MW / MW); 
                double x = model.getMatching().x_t(model.getMut());
                    
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
        case(LO_QED): // e
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
        case(NLO_QED11): // es
            {
            double mt = model.getMatching().getMt_mut();
            double MH = model.getMHl(); 
            double MW = model.Mw();  
            double MZ = model.getMz(); 
            double LZ = log(model.getMuw() * model.getMuw() / MZ / MZ); 
            
            double mt2 = mt * mt;
            double MW2 = MW * MW;
            double MZ2 = MZ * MZ;
            double MH2 = MH * MH;
            double MH4 = MH2 * MH2;
            double sw2 = model.sW2_MSbar_Approx();
            double sw4 = sw2 * sw2;
            double cw2 = 1. - sw2;
    
            gslpp::vector<double> CWin(3, 0.);
            double Cp= -22. / 9. - 4. / 3. * LZ ;
            double Cm= -22. / 9. - 4. / 3. * LZ ;
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
    gslpp::matrix<double> G1(dimension,0.) ;
    gslpp::matrix<double> S1(dimension,0.) ;
    gslpp::matrix<double> J1(dimension,0.) ;
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

    G1  = V_inv.real()  *  ADM_1  * V.real() ;
    
    //S1
    for(unsigned int i = 0; i < G1.size_i(); i++){
        for (unsigned int j = 0; j < G1.size_j(); j++) {
            S1.assign(i , j, beta_1 / 2. / beta_0 / beta_0 * e(i).real() * (i==j) - G1(i,j) / (2. * beta_0 + e(i).real() - e(j).real() ));
        }
    }

    J1 = V.real() * S1 * V_inv.real() ;
    
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
        return J1;
    } else if (order==NNLO){
        gslpp::matrix<double> S2(dimension,0.) , G2(dimension,0.) , J2(dimension,0.);
        G2  = V_inv.real()  *  ADM_2  * V.real() ;
        for(unsigned int i = 0; i < G1.size_i(); i++){
            for (unsigned int j = 0; j < G1.size_j(); j++) {
                double term=0. ;
                for (unsigned int k = 0; k < G1.size_j(); k++){
                    term += ( 1. + e(i).real() / 2. / beta_0 - e(k).real() / 2. / beta_0 ) / ( 2. + e(i).real() / 2. / beta_0 - e(j).real() / 2. / beta_0 ) * ( S1(i,k) * S1(k,j) - beta_1 / beta_0 * S1(i,j) * (j==k) );
                }
                S2.assign(i , j , beta_2 / 4. / beta_0 / beta_0 * e(i).real() * (i==j) + term - G2(i,j) / (4. * beta_0 + e(i).real() - e(j).real() )); 
            }
        }
        J2 = V.real() * S2 * V_inv.real() ;
        return J2;
    
    } else{
        std::stringstream out ;
        out << order;
        throw std::runtime_error("Charm_Kpnunu::RGevol_J : order " + out.str() +" not implemented\n");   
    }
}
                              

gslpp::matrix<double> Charm_Kpnunu::RGevol_R(orders_qed order_qed, int nf, int contribution)
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
        throw std::runtime_error("Charm_Kpnunu::RGevol_K : contribution " + out.str() +" not implemented\n");
        break;
    }
    
    
    gslpp::matrix<gslpp::complex> V(dimension,0.)  , V_inv(dimension,0.) ; 
    gslpp::matrix<double> M_0(dimension,0.) , K_0(dimension,0.) ; 
    gslpp::matrix<double> ADM_0(dimension,0.) , ADM_1(dimension,0.) , ADM_2(dimension,0.), ADM_e(dimension,0.), ADM_es(dimension,0.);
    gslpp::matrix<double> M_1(dimension,0.), S_1(dimension,0.) , K1_1(dimension,0.) , K2_1(dimension,0.) , K3_1(dimension,0.), G_1(dimension,0.) ;
                
    
    gslpp::vector<gslpp::complex> e(dimension,0.) ;
    
    double beta_0  = 11. - 2. / 3. * nf ; 
    double beta_1  = 102. - 38. / 3. * nf ;
    //double beta_2 = 2857. / 2. - 5033. / 18. * nf + 325. / 54. * nf * nf  ;
    
    
    ADM_0  = ADM(LO,NO_QED,nf,contribution); // remember that these matrices are all transposed as in the reference
    ADM_e  = ADM(LO,LO_QED,nf,contribution); // remember that these matrices are all transposed as in the reference
    ADM_es = ADM(LO,NLO_QED11,nf,contribution); // remember that these matrices are all transposed as in the reference
    ADM_1  = ADM(NLO,NO_QED,nf,contribution);
    ADM_2  = ADM(NNLO,NO_QED,nf,contribution);
    
    ADM_0.eigensystem(V,e);
    V_inv=V.inverse();
            
    
    M_0 = V_inv.real() * ADM_e * V.real();
    M_1 = V_inv.real() * ( (ADM_es - beta_1 / beta_0 * ADM_e) + (ADM_e * RGevol_J(NLO,nf,contribution) - RGevol_J(NLO,nf,contribution) * ADM_e)  ) * V.real();
    G_1 = V_inv.real()  *  ADM_1  * V.real();
    //S1
    for(unsigned int i = 0; i < G_1.size_i(); i++){
        for (unsigned int j = 0; j < G_1.size_j(); j++) {
            S_1.assign(i , j, beta_1 / 2. / beta_0 / beta_0 * e(i).real() * (i==j) - G_1(i,j) / (2. * beta_0 + e(i).real() - e(j).real() ));
        }
    } 
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
                return ((-2. * M_PI / beta_0) * (V.real() * K_0 * V_inv.real()) );
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
                
                return ( ( -0.5 / beta_0) * V.real() * (K1_1 + K2_1 + K3_1) * V_inv.real()); 
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
                return ((-2. * M_PI / beta_0) * (V.real() * K_0 * V_inv.real()) );
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
                
                return ( ( -0.5 / beta_0) * V.real() * (K1_1 + K2_1 + K3_1) * V_inv.real()); 
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
    
    switch(contribution){
        case(0):{ //box
            switch(order_qed){
                case (NO_QED):{
                    switch (order)
                    {
                        case (LO):{
                            return U0_4b * U0_5b * CW0b ;
                            break;
                        }
                        case (NLO):{
                            return J1_4b * U0_4b * U0_5b * CW0b + etacb * U0_4b * (J1_5b - J1_4b) * U0_5b * CW0b +
                                    etab * etacb * U0_4b * U0_5b * (CW1b - J1_5b * CW0b) ;

                            break;
                        }
                        case(NNLO):{ 
                            return J2_4b * U0_4b * U0_5b * CW0b + etacb * J1_4b * U0_4b * (J1_5b - J1_4b) * U0_5b * CW0b +
                                    etab * etacb * J1_4b * U0_4b * U0_5b * (CW1b - J1_5b * CW0b) + etacb * etacb * U0_4b *
                                    (J2_5b - J2_4b - J1_4b * (J1_5b - J1_4b)) * U0_5b * CW0b + etab * etacb * etacb *
                                    U0_4b * (J1_5b - J1_4b) * U0_5b * (CW1b - J1_5b * CW0b) + etab * etab * etacb * etacb * U0_4b * U0_5b *
                                    (CW2b - J1_5b * CW1b - (J2_5b - J1_5b * J1_5b) * CW0b) - ThresholdCb(NNLO);

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
                case(LO_QED):{
                    return ( (U0_4b*R0_5b/(4.*M_PI) + R0_4b*U0_5b/(4.*M_PI))*CW0b + U0_4b*U0_5b*CWeb );
                    break;
                }
                case(NLO_QED11):{
                    return ( (U0_4b*R1_5b + (J1_4b*U0_4b*R0_5b)/(4*M_PI) - (U0_4b*J1_4b*R0_5b)/(4*M_PI) + etacb/(4*M_PI)*(R0_4b*J1_5b*U0_5b) -
                            (etacb*etab)/(4*M_PI)*(R0_4b*U0_5b*J1_5b) + R1_4b*U0_5b)*CW0b + ( (etacb*etab)/(4*M_PI)*U0_4b*R0_5b + (etacb*etab)/(4*M_PI)*R0_4b*U0_5b)*CW1b +
                            (U0_4b*U0_5b)*CWesb + (J1_4b*U0_4b*U0_5b + etacb*U0_4b*J1_5b*U0_5b - etacb*U0_4b*J1_4b*U0_5b - etacb*etab*U0_4b*U0_5b*J1_5b)*CWeb  );
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
        case(1):{
            switch(order_qed){
                case (NO_QED):{
                    switch (order)
                    {
                        case (LO):{
                            return U0_4p * U0_5p * CW0p ;
                            break;
                        }
                        case (NLO):{
                            return J1_4p * U0_4p * U0_5p * CW0p + etacb * U0_4p * (J1_5p - J1_4p) * U0_5p * CW0p +
                                    etab * etacb * U0_4p * U0_5p * (CW1p - J1_5p * CW0p) ;

                            break;
                        }
                        case(NNLO):{ 
                            return J2_4p * U0_4p * U0_5p * CW0p + etacb * J1_4p * U0_4p * (J1_5p - J1_4p) * U0_5p * CW0p +
                                    etab * etacb * J1_4p * U0_4p * U0_5p * (CW1p - J1_5p * CW0p) + etacb * etacb * U0_4p *
                                    (J2_5p - J2_4p - J1_4p * (J1_5p - J1_4p)) * U0_5p * CW0p + etab * etacb * etacb *
                                    U0_4p * (J1_5p - J1_4p) * U0_5p * (CW1p - J1_5p * CW0p) + etab * etab * etacb * etacb * U0_4p * U0_5p *
                                    (CW2p - J1_5p * CW1p - (J2_5p - J1_5p * J1_5p) * CW0p) - ThresholdCp(NNLO);

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
                case(LO_QED):{
                    return ( (U0_4p*R0_5p/(4.*M_PI) + R0_4p*U0_5p/(4.*M_PI))*CW0p + U0_4p*U0_5p*CWep  );
                    break;
                }
                case(NLO_QED11):{
                    return ( (U0_4p*R1_5p + (J1_4p*U0_4p*R0_5p)/(4*M_PI) - (U0_4p*J1_4p*R0_5p)/(4*M_PI) + etacb/(4*M_PI)*(R0_4p*J1_5p*U0_5p) -
                            (etacb*etab)/(4*M_PI)*(R0_4p*U0_5p*J1_5p) + R1_4p*U0_5p)*CW0p + ( (etacb*etab)/(4*M_PI)*U0_4p*R0_5p + (etacb*etab)/(4*M_PI)*R0_4p*U0_5p)*CW1p +
                            (U0_4p*U0_5p)*CWesp + (J1_4p*U0_4p*U0_5p + etacb*U0_4p*J1_5p*U0_5p - etacb*U0_4p*J1_4p*U0_5p - etacb*etab*U0_4p*U0_5p*J1_5p)*CWep  );
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
        default:{
            std::stringstream out ;
            out << contribution;
            throw std::runtime_error("Charm_Kpnunu::C : contribution " + out.str() +" not implemented\n");
            break;
        }
    }
       
}