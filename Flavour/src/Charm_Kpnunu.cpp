/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Charm_Kpnunu.h"

Charm_Kpnunu::Charm_Kpnunu(const StandardModel& model_i) 
: model(model_i), modelmatching(*(model_i.getMyMatching())),
cp(3,0.), dcp(3,0.), c_p(3,0.), cpmuW0(3,0.), cpmuW1(3,0.), cpmuW2(3,0.),  
cb(2,0.), dcb(2,0.), c_b(2,0.), cbmuW0(2,0.), cbmuW1(2,0.), cbmuW2(2,0.),
U4p(3,0.), U5p(3,0.), J5p1(3,0.), J4p1(3,0.), J5p2(3,0.), J4p2(3,0.), dc_p(3,0.),
U4b(2,0.), U5b(2,0.), J5b1(2,0.), J4b1(2,0.), J5b2(2,0.), J4b2(2,0.), dc_b(2,0.) 
{

CP = 0.;    
CBe = 0.;
CBt = 0.;
etab = model.Als(model.getMuw())/model.Als(model.getMub());
etacb = model.Als(model.getMub())/model.Als(model.getMuc());
// WARNING: we need to put mc(mc) for xc, line right below to be checked!
mc = model.Mrun(model.getMuc(), model.getQuarks(QCD::CHARM).getMass_scale(),
     model.getQuarks(QCD::CHARM).getMass(), FULLNNLO);
etac = model.Als(model.getMuc())/model.Als(mc);
kc = pow(etac, 24./25.);
xi1c = 15212./1875.*(etac - 1.)/etac;
xi2c = 966966391./10546875. - 231404944./3515625.*(1./etac) - 272751559./10546875.*
       (1./etac)*(1./etac) - 128./5.*(1.-(1./etac)*(1./etac))*gsl_sf_zeta_int(3);
xc = kc*(1. + model.Als(mc)/4./M_PI*xi1c + (model.Als(mc)/4./M_PI)*(model.Als(mc)/4./M_PI)*xi2c);
}

Charm_Kpnunu::~Charm_Kpnunu() 
{}

gslpp::vector<double> Charm_Kpnunu::Cp(orders order)
{
    
    double x = modelmatching.x_t(model.getMuw());
    double L = log(model.getMuw()*model.getMuw()/model.Mw_tree()/model.Mw_tree());
    
    switch(order){
        case(LO):
            cp(0) = 4.;
            cp(1) = 4.;
            cp(2) = 0.;
            
            return(cp);
            
        case(NLO):
            cp(0) = 4./3.*(11. + 6.*L);
            cp(1) = -8./3.*(11. + 6.*L);
            cp(2) = 8.*(2. + L);
            
            return(cp);
            
        case(NNLO):
            cp(0) = 4.*(-(135677.-124095.)/3600. + 58./18.*M_PI*M_PI-0.5*(2./3.)*
                    (112./9. + 32.*x +(20./3. + 16.*x)*log(x) - (8.+16.*x)*
                    sqrt(4.*x - 1.)*gsl_sf_clausen(2.*asin(1./(2.*sqrt(x)))))
                    +(5./36.*238.*L) + 58./6.*L*L);
            cp(1) = 4.*(-(135677.+124095.)/3600. - 44./18.*M_PI*M_PI+0.5*(4./3.)*
                    (112./9. + 32.*x +(20./3. + 16.*x)*log(x) - (8.+16.*x)*
                    sqrt(4.*x - 1.)*gsl_sf_clausen(2.*asin(1./(2.*sqrt(x)))))
                    - (5./36.*260.*L) -44./6.*L*L);
            cp(2) = 4.*model.getCF()*(33.+4.*M_PI*M_PI+34.*L + 12.*L*L);
            
            return(cp);
            
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::Cp() order" + out.str() + " not implemented");
    }
}

gslpp::matrix<double> Charm_Kpnunu::RGevolP(int nf, orders order) 
{
    
    gslpp::matrix<double> evo(3,3, 0.);
    
    switch (order){
        case (LO):
            switch (nf){
                case(5):
                    evo(0,0) = pow(etab, 6./23.);
                    evo(1,1) = pow(etab, -12./23.);          
                    evo(2,0) = 12./5.*(pow(etab, 6./23.) - pow(etab, 1./23.));
                    evo(2,1) = 6./13.*(pow(etab, -12./23.) - pow(etab, 1./23.));
                    evo(2,2) = pow(etab, 1./23.);

                    return(evo);

                case(4):
                    evo(0,0) = pow(etacb, 6./25.);
                    evo(1,1) = pow(etacb, -12./25.);
                    evo(2,0) = 12./7.*(pow(etacb, 6./25.) - pow(etacb, -1./25.));
                    evo(2,1) = 6./11.*(pow(etacb, -12./25.) - pow(etacb, -1./25.));
                    evo(2,2) = pow(etacb, -1./25.);

                    return(evo); 
                    
                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolP() " + out.str() + " wrong number of falvours");
        }       
        case(NLO):
            switch(nf){
                case(5):
                    evo(0,0) =  5165./3174.;
                    evo(1,1) = -2267./1587.;   
                    evo(2,0) = -15857./1587.;
                    evo(2,1) = 15305./3174.;
                    evo(2,2) = -14924./1587.;
                    
                    return(evo);
                    
                case(4):
                    evo(0,0) =  6719./3750.;
                    evo(1,1) = -3569./1875.; 
                    evo(2,0) = -15931./1875.;
                    evo(2,1) = 5427./1250.;         
                    evo(2,2) = -15212./1875.;
                    
                    return(evo);
                    
                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolP() " + out.str() + " wrong number of flavours");
        }                
        case(NNLO):
            switch(nf){
                case(5):
                    evo(0,0) = -7.35665;
                    evo(1,1) = -54.9107;   
                    evo(2,0) = 17.7699;
                    evo(2,1) = -1.7514;
                    evo(2,2) = 18.3025;
                    
                    return(evo);
                    
                case(4):
                    evo(0,0) =  -10.2451;
                    evo(1,1) = -50.3422;
                    evo(2,0) = 8.0325;
                    evo(2,1) = -0.3657;
                    evo(2,2) = 4.91177;
                    
                    return(evo);
                    
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
    double L = log(mub*mub/Mb/Mb);
    
    switch (order){
        case(LO):
            dcp = 0.;
            
            return(dcp);
        
        case(NLO):
            dcp = 0.;
            
            return(dcp);
        
        case(NNLO):
            dcp(0) = -pow(etab, 6./23.)*(2./3.*L*((631.+9699.)/6348.*(1-etab)+etab*Cp(NLO)(0))
                     -2./3.*(59./36.+1./3.*L+L*L));
            dcp(1) = -pow(etab, 8./23.)*(2./3.*L*((631.-9699.)/6348.*(1-etab)+etab*Cp(NLO)(1))
                     +4./3.*(59./36.+1./3.*L+L*L));
            dcp(2) = (-2./3.)*L*((284704./2645. + 694522./20631.*etab)*pow(etab, 1./23.) 
                     -(1033492./7935.+ 8264./529.*etab)*pow(etab, 6./23.) + (3058./1587.+ 18136./6877.*etab)
                     *pow(etab, -12./23.) + etab*(pow(etab, 1./23.)*Cp(NLO)(2) + 48./5.*(pow(etab, 6./23.)
                     -pow(etab, 1./23.))*Cp(NLO)(0) + 24./13.*(pow(etab, -12./23.)-pow(etab, 1./23.))*Cp(NLO)(1)));
            
            return(dcp);
            
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
    
    U4p = RGevolP(4,LO);
    U5p = RGevolP(5,LO);
    J4p1 = RGevolP(4,NLO);
    J5p1 = RGevolP(5,NLO);
    J4p2 = RGevolP(4,NNLO);
    J5p2 = RGevolP(5,NLO);
    
    for(int i=0; i<3; i++){
    dc_p(i,i) =  ThresholdCp(NNLO)(i);
    }
    
    switch(order){
        case(LO):
            c_p = U4p*U5p*cpmuW0;
            
            return(c_p);
            
        case(NLO):
            c_p = J4p1*U4p*U5p*cpmuW0 + etacb*U4p*(J5p1 - J4p1)*U5p*cpmuW0 +
                  etab*etacb*U4p*U5p*(cpmuW1-J5p1*cpmuW0);
                    
            return(c_p);
            
        case(NNLO):
            //threshold term need to be checked (not really clear in (86) of hep-ph/0603079v2!)
            c_p = J4p2*U4p*U5p*cpmuW0 + etacb*J4p1*U4p*(J5p1 - J4p1)*U5p*cpmuW0 +
                  etab*etacb*J4p1*U4p*U5p*(cpmuW1-J5p1*cpmuW0) + etacb*etacb*U4p*
                  (J5p2 - J4p2 - J4p1*(J5p1-J4p1) - dc_p)*U5p*cpmuW0 + etab*etacb*etacb*
                  U4p*(J5p1-J4p1)*U5p*(cpmuW1-J5p1*cpmuW0) + etab*etab*etacb*etacb*U4p*U5p*
                  (cpmuW2-J5p1*cpmuW1-(J5p2-J5p1*J5p1)*cpmuW0);
            
            return(c_p);
            
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_p() order" + out.str() + " not implemented");
    }
    
}

double Charm_Kpnunu::C_P(orders order)
{
    
    double L = log(model.getMuc()*model.getMuc()/mc/mc);
    double C_P0 = C_p(LO)(2);
    double rhoP1p =  4*(1.- L) + 4.*log(kc);
    double rhoP1m = -2.*(1.- L) - 2*log(kc);
    double C_P1 = C_p(NLO)(2) + 4.*(C_p(LO)(0)*rhoP1p + C_p(LO)(1)*rhoP1m) + xi1c*C_p(LO)(2);
    double rhoP2p = 11. + 20.*L - 12.*L*L - 20.*log(kc) - 12.*log(kc)*log(kc) + 
            24.*log(kc)*L + 4.*xi1c;
    double rhoP2m = -7. + 13.*L + 12.*L*L - 12.*log(kc) + 12.*log(kc)*log(kc)
            - 24.*log(kc)*L - 2.*xi1c;
    double C_P2 = C_p(NNLO)(2) + 4.*(C_p(NLO)(0)*rhoP1p + C_p(NLO)(1)*rhoP1p + C_p(LO)(0)*rhoP2p 
          + C_p(LO)(1)*rhoP2m) + xi1c*(C_p(NLO)(2) + 4.*(C_p(LO)(0)*rhoP1p + C_p(LO)(1)
          *rhoP1m)) + xi2c*C_p(LO)(2);
    
    switch(order){
        case(LO):
            CP = kc*xc/32.*(4.*M_PI/model.Als(model.getMuc())*C_P0);
            
            return(CP);
            
        case(NLO):
            CP = kc*xc/32.*(4.*M_PI/model.Als(model.getMuc())*C_P0 + C_P1);
            
            return(CP);
            
        case(NNLO):
            CP = kc*xc/32.*(4.*M_PI/model.Als(model.getMuc())*C_P0 + C_P1 + model.Als(model.getMuc())
                /4./M_PI*C_P2);
            
            return(CP); 
            
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_P() order" + out.str() + " not implemented");
    }
}
 
gslpp::vector<double> Charm_Kpnunu::Cb(orders order)
{
    
    double L = log(model.getMuw()*model.getMuw()/model.Mw_tree()/model.Mw_tree());
    
    switch(order){
        case(LO):
            cb(0) = 4.;
            cb(1) = 0.;
            
            return(cb);
            
        case(NLO):
            cb(0) = 0.;
            cb(1) = -4.*(9. + 4.*L);
            
            return(cb);
            
        case(NNLO):
            cb(0) = 0.;
            cb(1) = -8.*model.getCF()*(20. + 2.*M_PI*M_PI + 25.*L + 6.*L*L);
            
            return(cb);
            
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::Cb() order" + out.str() + " not implemented");
    }
}

gslpp::matrix<double> Charm_Kpnunu::RGevolB(int nf, orders order)
{
    
    gslpp::matrix<double> evo(2,2, 0.);
    
    switch (order){
        case (LO):
            switch (nf){
                case(5):
                    evo(0,0) = 1.;
                    evo(0,1) = 0.;          
                    evo(1,0) = 12.*(1. - pow(etab, 1./23.));
                    evo(1,1) = pow(etab, 1./23.);

                    return(evo);

                case(4):
                    evo(0,0) = 1.;
                    evo(0,1) = 0.;
                    evo(1,0) = -12.*(1. - pow(etacb, -1./25.));
                    evo(1,1) = pow(etacb, -1./25.);

                    return(evo); 
                    
                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolB() " + out.str() + " wrong number of falvours");
        }       
        case(NLO):
            switch(nf){
                case(5):
                    evo(0,0) = 0.;
                    evo(0,1) = 0.;   
                    evo(1,0) = 2402./1587.;
                    evo(1,1) = -14924./1587.;
                    
                    return(evo);
                    
                case(4):
                    evo(0,0) = 0.;
                    evo(0,1) = 0.; 
                    evo(1,0) = 581./1875.;
                    evo(1,1) = -15212./1875.;
                    
                    return(evo);
                    
                default:
                    std::stringstream out;
                    out << nf;
                    throw std::runtime_error("Charm_Kpnunu::RgevolB() " + out.str() + " wrong number of flavours");
        }                
        case(NNLO):
            switch(nf){
                case(5):
                    evo(0,0) = 0.;
                    evo(0,1) = 0.;   
                    evo(1,0) = 1296371522./39457581. - 34624./1081.*gsl_sf_zeta_int(3);
                    evo(1,1) = -177621017./7555707. + 800./23.*gsl_sf_zeta_int(3);
                    
                    return(evo);
                    
                case(4):
                    evo(0,0) = 0.;
                    evo(0,1) = 0.;
                    evo(1,0) = 684990354./19140625. - 6976./245.*gsl_sf_zeta_int(3);
                    evo(1,1) = -272751559./10546875. + 128./5.*gsl_sf_zeta_int(3);
                    
                    return(evo);
                    
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
    double L = log(mub*mub/Mb/Mb);
    
    switch (order){
        case(LO):
            dcb = 0.;
            
            return(dcb);
        
        case(NLO):
            dcb = 0.;
            
            return(dcb);
        
        case(NNLO):
            dcb(0) = 0.;
            dcb(1) = -2./3.*L*((238784./529.-9608./1587*etab)*pow(etab,1./23.) -1336./3. +
                     pow(etab,24./23.)*Cb(NLO)(1));
                  
            return(dcb);
            
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
    
    U4b = RGevolB(4,LO);
    U5b = RGevolB(5,LO);
    J4b1 = RGevolB(4,NLO);
    J5b1 = RGevolB(5,NLO);
    J4b2 = RGevolB(4,NNLO);
    J5b2 = RGevolB(5,NLO);
    
    for(int i=0; i<3; i++){
    dc_b(i,i) =  ThresholdCb(NNLO)(i);
    }
    
    switch(order){
        case(LO):
            c_b = U4b*U5b*cbmuW0;
            
            return(c_b);
            
        case(NLO):
            c_b = J4b1*U4b*U5b*cbmuW0 + etacb*U4b*(J5b1 - J4b1)*U5b*cbmuW0 +
                  etab*etacb*U4b*U5b*(cbmuW1-J5b1*cbmuW0);
                    
            return(c_b);
            
        case(NNLO):
            //threshold term need to be checked (not really clear in (86) of hep-ph/0603079v2!)
            c_b = J4b2*U4b*U5b*cbmuW0 + etacb*J4b1*U4b*(J5b1 - J4b1)*U5b*cbmuW0 +
                  etab*etacb*J4b1*U4b*U5b*(cbmuW1-J5b1*cbmuW0) + etacb*etacb*U4b*
                  (J5b2 - J4b2 - J4b1*(J5b1-J4b1) - dc_b)*U5b*cbmuW0 + etab*etacb*etacb*
                  U4b*(J5b1-J4b1)*U5b*(cbmuW1-J5b1*cbmuW0) + etab*etab*etacb*etacb*U4b*U5b*
                  (cbmuW2-J5b1*cbmuW1-(J5b2-J5b1*J5b1)*cbmuW0);
            
            return(c_b);
            
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_p() order" + out.str() + " not implemented");
    }
    
}

double Charm_Kpnunu::C_Be(orders order)
{
    
    double L = log(model.getMuc()*model.getMuc()/mc/mc);
    double rhoB1_e = 5. + 4.*L - 4.*log(kc);
    double rhoB2_e = -2.*model.getCF()*(9. - L - 6.*L*L) -8./3.*log(kc) + 16.*log(kc)*log(kc) 
                     - 32.*log(kc)*L - 4.*xi1c;
    double C_B0e = C_b(LO)(1);
    double C_B1e = C_b(NLO)(1) + 4.*rhoB1_e + xi1c*C_b(LO)(1);
    double C_B2e = C_b(NNLO)(1) + 4.*rhoB2_e + xi1c*C_b(NLO)(1) + xi2c*C_b(LO)(1);

    switch(order){
        case(LO):
            CBe = kc*xc/32.*(4.*M_PI/model.Als(model.getMuc())*C_B0e);
            
            return(CBe);
            
        case(NLO):
            CBe = kc*xc/32.*(4.*M_PI/model.Als(model.getMuc())*C_B0e + C_B1e);
            
            return(CBe); 
            
        case(NNLO):
            CBe = kc*xc/32.*(4.*M_PI/model.Als(model.getMuc())*C_B0e + C_B1e + model.Als(model.getMuc())/4./M_PI*C_B2e);
            
            return(CBe); 
            
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_Be() order" + out.str() + " not implemented");
    }
}

double Charm_Kpnunu::C_Bt(orders order)
{
     
    double L = log(model.getMuc()*model.getMuc()/mc/mc);
    double x_t = model.getLeptons(model.TAU).getMass()*model.getLeptons(model.TAU).getMass()/mc/mc;
   
    double C_B0t = C_b(LO)(1);
    double rhoB1_t = 5. + 4.*L - 4.*log(kc) + 4./(x_t - kc)*(kc*log(kc) 
                     - x_t*(1. - kc)/(1. - x_t)*log(x_t));
    double C_B1t = C_b(NLO)(1) + 4.*rhoB1_t + xi1c*C_b(LO)(1);
    double rhoB2_t = -2.*model.getCF()*((9. + 7.*x_t)/(1. - x_t) + x_t*(3. + 13.*x_t)/
                 (1. - x_t)/(1. - x_t) - 12.*x_t/(1. - x_t)*gsl_sf_dilog(1.-x_t) 
                 - ((1. - 13.*x_t)/(1. - x_t) - 12.*x_t*x_t/(1. - x_t)/(1. - x_t)*log(x_t))
                 *L - 6.*L*L) + 32./(x_t - kc)*(4.*x_t*(1 - kc)/3./(1.-x_t) - x_t
                 *((x_t*(13. - 29.*x_t)) + kc*(3. + 29.*x_t*x_t) -kc*kc*(3. + 13.*x_t))
                 /12./(x_t - kc)/(1. - x_t)/(1. - x_t)*log(x_t) + kc*(17.*x_t - kc)/12.
                 /(x_t - kc)*log(kc) + x_t*x_t/(x_t - kc)*log(x_t)*log(kc) - (x_t*x_t + 
                 2.*x_t*kc - kc*kc)/2./(x_t - kc)*log(kc)*log(kc) - x_t*(x_t - kc)/2.
                 /(1. - x_t)*gsl_sf_dilog(1. - x_t) - x_t*gsl_sf_dilog(1. - x_t/kc)) 
                 + 32./(x_t - kc)/(1. - x_t)*(x_t*(1. - kc) + kc*(1. - x_t)*(2.*x_t - kc)
                 /(x_t - kc)*log(x_t) + x_t*x_t*(1. - kc)*(1. - 2.*x_t + kc)/(x_t - kc)
                 /(1. - x_t)*log(kc))*L + 4.*kc/(x_t - kc)*(1. - x_t/(x_t - kc)*log(x_t)
                 + x_t/(x_t - kc)*log(kc))*xi1c;
    double C_B2t = C_b(NNLO)(1) + 4.*rhoB2_t + xi1c*C_b(NLO)(1) + xi2c*C_b(LO)(1);
            
    switch(order){
        case(LO):
            CBt = kc*xc/16.*(4.*M_PI/model.Als(model.getMuc())*C_B0t);
            
            return(CBt);
            
        case(NLO):
            CBt = kc*xc/16.*(4.*M_PI/model.Als(model.getMuc())*C_B0t + C_B1t);
            
            return(CBt);
            
        case(NNLO):
            CBt = kc*xc/16.*(4.*M_PI/model.Als(model.getMuc())*C_B0t + C_B1t 
                  + model.Als(model.getMuc())/4./M_PI*C_B2t);
            
            return(CBt);
            
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Charm_Kpnunu::C_Bt() order" + out.str() + " not implemented");
    }  
}

double Charm_Kpnunu::P_C(orders order)
{
    double Xe = C_P(order) + C_Be(order);
    double Xt = C_P(order) + C_Bt(order);
    double lambda4 = model.getLambda()*model.getLambda()*model.getLambda()*model.getLambda();
    double pc = 1./lambda4*(2./3.*Xe + 1./3.*Xt);
    
    return(pc);    
}

double Charm_Kpnunu::C_TOT(orders order, orders_ew order_ew)
{
    
    double xt = modelmatching.x_t(model.getMut());
    double Muw = model.getMuw();
    double Ale = model.getAle();
    double a = 1./modelmatching.mt2omh2(Muw);
    gslpp::complex lambdat = model.computelamt();
    gslpp::complex lambdac = model.computelamc();
    double lambda = model.getLambda();
    double lambda5 = model.getLambda()*model.getLambda()*model.getLambda()
                     *model.getLambda()*model.getLambda();
    double IBT = model.getDeltaP_cu();
    double X = 0.;
    
    if((order == NNLO) && (order_ew == NLO_ew)){
            X = lambdat.real()/lambda5*(modelmatching.X0t(xt) + model.Als(model.getMut(), FULLNLO)
                /4./M_PI* modelmatching.X1t(xt) + Ale/4./M_PI*modelmatching.Xewt(xt, a, Muw));
            return(X + (lambdac.real()/lambda)*(P_C(NNLO) + IBT));
    }
    
    if((order == NLO) && (order_ew == LO_ew)){
            X = lambdat.real()/lambda5*(modelmatching.X0t(xt)/lambda5 + model.Als(model.getMut(), FULLNLO)/4./M_PI*modelmatching.X1t(xt));
            return(X + (lambdac.real()/lambda)*(P_C(NLO) + IBT));
    }
    
    if((order == LO) && (order_ew == LO_ew)){
            X = lambdat.real()*modelmatching.X0t(xt)/lambda5;
            return(X + (lambdac.real()/lambda)*(P_C(LO) + IBT));
    }

    else{
         std::stringstream out;
         out << order_ew;
         throw std::runtime_error("BRKppnunu::BRKppnunu(): order_ew " + out.str() + "not implemented\n");
         out << order;
         throw std::runtime_error("BRKppnunu::BRKppnunu(): order " + out.str() + "not implemented");
    }
}
  