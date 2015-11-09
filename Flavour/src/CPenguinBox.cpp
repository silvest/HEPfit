/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CPenguinBox.h"

CPenguinBox::CPenguinBox(const StandardModel& model_i) 
: model(model_i), modelmatching(*(model_i.getMyMatching())),
        up(3,0.), j5p1(3,0.), j4p1(3,0.), j5p2(3,0.), j4p2(1,0.),
        ub(2,0.), j5b1(2,0.), j4b1(2,0.), j5b2(2,0.), j4b2(2,0.),
        up4(3,0.), up5(3,0.), ub4(2,0.), ub5(2,0.),
        p(3,0.), b(2,0.),
        p0(3,0.), p1(3,0.), p2(3,0.), b0(2,0.), b1(2,0.), b2(2,0.),
        r1(3,0.), r2(3,0.), CPM0(3,0.), CPM1(3,0.), CPM2(3,0.),
        CBM0(3,0.), CBM1(3,0.), CBM2(3,0.) 
{
    
    j5p1(0,0) =  5165./3174.;
    j5p1(1,1) = -2267./1587.;          
    j5p1(2,2) = -14924./1587.;
    j5p1(2,0) = -15857./1587.;
    j5p1(2,1) = 15305./3174.;
    
    j5p2(0,0) = -7.35665;
    j5p2(1,1) = -54.9107;          
    j5p2(2,2) = 18.3025;
    j5p2(2,0) = 17.7699;
    j5p2(2,1) = -1.7514;
    
    j4p1(0,0) =  6719./3750.;
    j4p1(1,1) = -3569./1875.;          
    j4p1(2,2) = -15212./1875.;
    j4p1(2,0) = -15931./1875.;
    j4p1(2,1) = 5427./1250.;
    
    j4p2(0,0) =  -10.2451;
    j4p2(1,1) = -50.3422;
    j4p2(2,2) = 4.91177;
    j4p2(2,0) = 8.0325;
    j4p2(2,1) = -0.3657;
    
    j5b1(1,0) = 2402./1587.;
    j5b1(1,1) = -14924./1587.;
    
    j5b2(1,0) = 1296371522./39457581. - 34624./1081.*gsl_sf_zeta(3);
    j5b2(1,1) = -177621017./7555707. +800./23.*gsl_sf_zeta(3);
    
    j5b1(1,0) = 581./1875.;
    j5b1(1,1) = -15212./1875.;
    
    j5b2(1,0) = 684990354./19140625. - 6976./245.*gsl_sf_zeta(3);
    j5b2(1,1) = -272751559./10546875. + 128./5.*gsl_sf_zeta(3);
}

CPenguinBox::~CPenguinBox(){
    
}

gslpp::vector<double> CPenguinBox::WMatchP(orders order){
    double x = modelmatching.x_t(model.getMuw());
    
    double l = log(model.getMuw()*model.getMuw()/model.Mw_tree()/model.Mw_tree());
    
    switch(order){
        case(LO):
            p(0) = 4.;
            p(1) = 4.;
            p(2) = 0.;
            
            return(p);
        case(NLO):
            p(0) = 4.*0.5*2./3.*(11. + 6.*l);
            p(1) =-4.*0.5*4./3.*(11. + 6.*l);
            p(2) = 8.*(2. + l);
            
            return(p);
        case(NNLO):
            p(0) = 4.*(-(135677.-124095.)/3600. + 58./18.*M_PI*M_PI-0.5*(2./3.)*
                    ( 112./9. + 34.*x +(20./3. + 16.*x)*log(x) - (8.+16.*x)*
                    sqrt(4.*x - 1.)*gsl_sf_clausen(2.*asin(1./(2.*sqrt(x))))))
                    + (5./36.*238.*l) +58./6.*l*l;
            p(1) = -4.*((135677.+124095.)/3600. - 44./18.*M_PI*M_PI+0.5*(4./3.)*
                    ( 112./9. + 34.*x +(20./3. + 16.*x)*log(x) - (8.+16.*x)*
                    sqrt(4.*x - 1.)*gsl_sf_clausen(2.*asin(1./(2.*sqrt(x))))))
                    + (-5./36.*260.*l) -44./6.*l*l;
            p(2) = 4.*model.getCF()*(33.+4.*M_PI*M_PI+34.*l + 12.*l*l);
            
            return(p);
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("CPenguinBox::WMatchB() order" + out.str() + " not implemented");
    }
}

gslpp::matrix<double> CPenguinBox::RGevolP(int nf) 
{
    
    switch (nf){
        case(5):
            eta = model.Als(model.getMuw())/model.Als(model.getMub());
            
            up(0,0) = pow(eta, 6./23.);
            up(1,1) = pow(eta, -12./23.);          
            up(2,2) = pow(eta, 1./23.);
            up(2,0) = 12./5.*(pow(eta, 6./23.) - pow(eta, 1./23.));
            up(2,1) = 6./13.*(pow(eta, -12./23.) - pow(eta, 1./23.));
            
            return(up);
        case(4):
            eta = model.Als(model.getMub())/model.Als(model.getMuc());
            
            up(0,0) = pow(eta, 6./25.);
            up(1,1) = pow(eta, -12./25.);          
            up(2,2) = pow(eta, -1./25.);
            up(2,0) = 12./7.*(pow(eta, 6./25.) - pow(eta, -1./25.));
            up(2,1) = 6./11.*(pow(eta, -12./25.) - pow(eta, -1./25.));
            
            return(up);
        default:
            std::stringstream out;
            out << nf;
            throw std::runtime_error("CPenguinBox::RgevolP() " + out.str() + " wrong number of falvours");
    }
}

double CPenguinBox::BmatchP(orders order)
{
    
    double mub = model.getMub();
    double Mb = model.Mrun(model.getMub(), model.getQuarks(QCD::BOTTOM).getMass_scale(),
                           model.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    eta = model.Als(model.getMuw())/mub;

    switch (order){
        case(NNLO):
            return((-2./3.)*log(mub*mub/Mb/Mb)
                    *((284704./2645. + 694522./20631.*eta)*pow(eta, 1./23.) - (1033492./7935. 
                    + 8264./529.*eta)*pow(eta, 6./23.) + (3058./1587. 
                    + 18136./6877.*eta)*pow(eta, 12./23.) + eta*(pow(eta, 1./23.)*WMatchP(NLO)(2)
                    + 48./5.*(pow(eta, 6./23.)-pow(eta, 1./23.))*WMatchP(NLO)(0) 
                    + 24./13.*(pow(eta, -12./23.)-pow(eta, 1./23.))*WMatchP(NLO)(1))));
        case(NLO):
            return(0.);
        case(LO):
            return(0.);
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("CPenguinBox::BmatchB() order " + out.str() + " not implemented");
    }
}

gslpp::vector<double> CPenguinBox::WMatchB(orders order)
{
    double l = log(model.getMuw()*model.getMuw()/model.Mw_tree()/model.Mw_tree());
    
    switch(order){
        case(LO):
            b(0) = 0.;
            b(1) = 4.;
            
            return(b);
        case(NLO):
            b(0) = -4.*(9.+4.*l);
            b(1) = 4.;
            
            return(b);
        case(NNLO):
            b(0) = -8.*model.getCF()*(20.+2.*M_PI*M_PI+25.*l+6.*l*l);
            b(1) = 4.;
            
            return(b);
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("CPenguinBox::WMatchB() order" + out.str() + " not implemented");
    }
}

gslpp::matrix<double> CPenguinBox::RGevolB(int nf) 
{
    
    switch (nf){
        case(5):
            eta = model.Als(model.getMuw())/model.Als(model.getMub());
            
            ub(0,0) = 1.;
            ub(1,1) = pow(eta, 1./23.);          
            ub(1,0) = 12.*(1.-pow(eta, 1./23.));
            
            return(ub);
        case(4):
            eta = model.Als(model.getMub())/model.Als(model.getMuc());
            
            ub(0,0) = 1.;
            ub(1,1) = pow(eta, -1./25.);          
            ub(1,0) = -12.*(1.-pow(eta, -1./25.));
            
            return(ub);
        default:
            std::stringstream out;
            out << nf;
            throw std::runtime_error("CPenguinBox::RgevolP() " + out.str() + " wrong number of falvours");
    }
}

double CPenguinBox::BmatchB(orders order)
{
    
    double mub = model.getMub();
    double Mb = model.Mrun(model.getMub(), model.getQuarks(QCD::BOTTOM).getMass_scale(),
                           model.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    eta = model.Als(model.getMuw())/mub;            
    
    switch (order){
        case(NNLO):
            return( -2./3.*log(mub*mub/Mb/Mb)*(
                    (238784./529. - 9608./1587.*eta)*pow(eta,1./23.) - 1336./3. 
                    + pow(eta,24./23.)*WMatchB(NLO)(1) ) );
        case(NLO):
            return(0.);
        case(LO):
            return(0.);
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("CPenguinBox::BmatchB() order " + out.str() + " not implemented");
    }
}

double CPenguinBox::Cmatch(orders order)
{
    
    up5 = RGevolP(5);
    up4 = RGevolP(4);
    ub5 = RGevolB(5);
    ub4 = RGevolB(4);
    
    p0 = WMatchP(LO);
    b0 = WMatchB(LO);
    
    double Muc = model.getMuc();
    double Muw = model.getMuw();
    double mc = model.Mrun(model.getMuc(), model.getQuarks(QCD::CHARM).getMass_scale(),
                           model.getQuarks(QCD::CHARM).getMass(), FULLNNLO);
    double l=log(Muc*Muc/mc/mc);
    double A = model.Als(Muc, FULLNNLO);
    
    //double etac = model.Als(Muc)/model.Als(mc);
    double etab = model.Als(Muw)/model.Als(model.getMub());
    double etacb = model.Als(model.getMub())/model.Als(Muc);
    double xtau = pow(1.7769/Muw,2.);
    //double kc = pow(etac, 24./25.);
    
    r1(0) = 4.*(1.-l);
    r1(1) = -2.*(1.-l);
    
    r2(0) = 11.+20.*l-12.*l*l;
    r2(0) = -7.-12.*l+12.*l*l;
    
    r1tau = 5. + 4.*xtau/(1-xtau)*log(xtau) + 4.*l;
    r2tau = -2.*model.getCF()*((9.+7.*xtau)/(1.-xtau) + (3.*xtau+13.*xtau*xtau)/(1.-xtau)/(1.-xtau)
            -12.*xtau*gsl_sf_dilog(1.-xtau)/(1.-xtau)
            -((1.-13.*xtau)/(1.-xtau) - 12.*xtau*xtau/(1.-xtau)/(1.-xtau)*log(xtau))*l-6.*l*l);
    
    r1e = 5.+4.*l;
    r2e = -2.*model.getCF()*(9.-l-6.*l*l);
    
    /*
     * this commented lines gives the matching in the approximation given in
     * hep-ph/060379
     * 
     * 
    rho1(0) = (4.*(1-2.*log(Muc/mc)) 
              + 4.*log(kc));
    rho1(1) = (-2.*(1-2.*log(Muc/mc)) 
              - 2.*log(kc));
    
    rho2(0) = 11.-40.*log(Muc/mc)
              + 24.*Muc/mc
                   *Muc/mc
              -20.*log(kc)-12.*log(kc)*log(kc)
              +48.*log(kc)*log(Muc/mc)
              +4.*15212./1875.*(1.-1/etac);
    
    rho2(1) = -7.-24.*log(Muc/mc)
              - 24.*Muc/mc
                   *Muc/mc
              +12.*log(kc)+12.*log(kc)*log(kc)
              -48.*log(kc)*log(Muc/mc)
              -2.*15212./1875.*(1.-1/etac);
    
    rho1tau = 5. + 4.*xtau/(1-xtau)*log(xtau) 
                 + 16.*log(Muc/mc) 
                 + 4./(xtau-kc)*(kc*log(kc)-xtau*(1.-kc)/(1.-xtau)*log(xtau));
    
    rho2tau = -2.*model.getCF()*((9.+7.*xtau)/(1.-xtau) + (3.*xtau+13.*xtau*xtau)/(1.-xtau)/(1.-xtau)
              -12.*xtau*gsl_sf_dilog(1.-xtau)/(1.-xtau)
              -((1.-13.*xtau)/(1.-xtau) - 12.*xtau*xtau/(1.-xtau)/(1.-xtau)*log(xtau))*2.*log(Muc/mc)
              -6.*log(Muc*Muc/mc/mc)*log(Muc*Muc/mc/mc))
              + 32./(xtau-kc)*(4.*xtau*(1.-kc)/3.*(1.-xtau) - 
                - xtau*(xtau*13.-29.*xtau*xtau+kc*(3.+29*xtau*xtau)-kc*kc*(3.+13.*xtau))/
                (12.*(xtau-kc)*(1.-xtau)*(1.-xtau))*log(xtau) 
                + kc*(17.*xtau-kc)/12./(xtau-kc)*log(kc)
                + xtau*xtau/(xtau-kc)*log(xtau)+log(kc)
                - (xtau*xtau+2.*xtau*kc-kc*kc)/(2.*(xtau-kc)*log(kc)*log(kc))
                - (xtau*(xtau-kc))/2./(1.-xtau)*gsl_sf_dilog(1.-xtau)
                - xtau*gsl_sf_dilog(1.-xtau/kc))
              + 32./(xtau-kc)/(1.-xtau)*log(Muc*Muc/mc/mc)*(
                xtau*(1.-kc)
                + kc*(1.-xtau)*(2.*xtau-kc)/(xtau-kc)*log(xtau)
                + xtau*xtau*(1.-kc)*(1.-2.*xtau+kc)/(xtau-kc)/(1.-xtau)*log(kc))
              +4.*kc/(xtau-kc)*15212./1875.*(1.-1./etac)*(
                1.-xtau/(xtau-kc)*log(xtau)+xtau/(xtau-kc)*log(kc));
    
    rho1e = 5.+4.*log(Muc+Muc/mc/mc)-4.*log(kc);
    
    rho2e = -2.*model.getCF()*(9.-log(Muc/Muc/mc/mc)-6.*log(Muc/Muc/mc/mc)*log(Muc/Muc/mc/mc))
            -8./3.*log(kc)+16.*log(kc)*log(kc) - 32.*log(kc)*log(Muc/Muc/mc/mc)
            -4.*15212./1875.*(1.-1./etac);
    */
    
    switch(order){
        case(NNLO):
            p2 = WMatchP(NNLO);
            p1 = WMatchP(NLO);
            
            b2 = WMatchB(NNLO);
            b1 = WMatchB(NLO);
            
            CPM0 = up4 * up5 * p0;
            CBM0 = ub4 * ub5 * b0;
            
            CPM1 = j4p1 * up4 * up5 * p0
                + etacb * up4 * (j5p1 -j4p1) * up5 *  p0
                + etab * etacb * up4 * up5 * (p1 - j5p1 * p0) ;
            CBM1 = j4b1 * ub4 * ub5 * b0
                + etacb * ub4 * (j5b1 -j4b1) * ub5 *  b0
                + etab * etacb * ub4 * ub5 * (b1 - j5b1 * b0) ;
            
            CPM2 = j4p2 * up4 * up5 * p0
                + etacb * j4p1 * up4 * (j5p1 - j4p1) * up5 * p0
                + etab * etacb * j4p1 * up4 * up5 * (p1 - j5p1 * p0)
                + etacb * etacb * up4 * (j5p2 - j4p2 - j4p1 * (j5p1 - j4p1) - BmatchP(NNLO)) * up5 * p0
                + etab * etacb * etacb * up4 * (j5p1 - j4p1) * up5 * (p1 - j5p1 * p0)
                + etab * etab * etacb * etacb * up4 * up5 * (p2 - j5p1 * p1 - (j5p2 - j5p1 * j5p1) * p0);
            CBM2 = j4b2 * ub4 * ub5 * b0
                + etacb * j4b1 * ub4 * (j5b1 - j4b1) * ub5 * b0
                + etab * etacb * j4b1 * ub4 * ub5 * (b1 - j5b1 * b0)
                + etacb * etacb * ub4 * (j5b2 - j4b2 - j4b1 * (j5b1 - j4b1) - BmatchB(NNLO)) * ub5 * b0
                + etab * etacb * etacb * ub4 * (j5b1 - j4b1) * ub5 * (b1 - j5b1 * b0)
                + etab * etab * etacb * etacb * ub4 * ub5 * (b2 - j5b1 * b1 - (j5b2 - j5b1 * j5b1) * b0);
            
            cpm = pow(mc/Muw, 2.)/32. * 4.*M_PI/A * (
                    CPM0(2) +

                    A/4./M_PI * (CPM1(2) + 4.*CPM0(0)*r1(0) + 4.*CPM0(1)*r1(1) ) +
                    
                    A/4./M_PI*A/4./M_PI*(CPM2(2) + 4.*CPM0(0)*r2(0) + 4.*CPM1(0)*r1(0)
                             + 4.*CPM0(1)*r2(1) + 4.*CPM1(1)*r1(1))
                    );
            
            cbmt = pow(mc/Muw, 2.)/16. * 4.*M_PI/A * (
                    CBM0(1) +
                    
                    A/4./M_PI * (CBM1(1) + 4.*CBM0(0)*r1tau) +
                    
                    A/4./M_PI*A/4./M_PI * (CBM2(1) + 4.*CBM0(0)*r2tau + 4.*CBM1(0)*r1tau)
                    );
            
            cbme = pow(mc/Muw, 2.)/16. * 4.*M_PI/A * (
                    CBM0(1) +
                    
                    A/4./M_PI * (CBM1(1) + 4.*CBM0(0)*r1e) +
                    
                    A/4./M_PI*A/4./M_PI*(CBM2(1) + 4.*CBM0(0)*r2e + 4.*CBM1(0)*r1e));
                    
            /*
            cpm = kc * pow(mc/Muw, 2.) / 32. *(
                    
                    4.*M_PI/model.Als(model.getMuc(), 4, FULLNNLO) *CPM0(2) +
                    
                    CPM1(2) + 
                    4.*CPM0(0)*r1(0) +
                    4.*CPM1(0)*r1(1) +
                    15212./1875.*(1.-1/etac)*CPM0(2) +
                    
                    model.Als(Muc)/4./M_PI*(
                    CPM2(2) + 4.*CPM1(0)*r1(0) + 4*CPM0(0)*r2(0) + 4.*CPM1(1)*r1(1) + 4*CPM0(1)*r2(1) +
                    15212./1872.*(1-1./kc)*CPM1(2) +
                    4.*15212./1872.*(1.-1./kc)*(CPM0(0)*r1(0) + CPM0(1)*r1(1)) +
                    CPM0(2)*(966966391./10546875. -231404944./3515625./etac-272751559./10546875./etac/etac -
                    128./5.*(1.-1./etac/etac)*gsl_sf_zeta(3))));
            
            cbmt = kc * pow(mc/Muw, 2.) / 16. * (
                    
                   4.*M_PI/model.Als(Muc)*CBM0(1)
                    
                   + CBM1(1) + 4.*r1tau*CBM1(0) + 15212./1875.*(1.-1./etac)*CBM0(1)
                    
                   + model.Als(Muc)/4./M_PI * (CBM2(1) + 4.*r2tau*CBM0(0)
                   + 15212./1875.*(1.-1./etac)*CBM1(1)
                   + ( 966966391./10546875. - 231404944./3515625./etac - 272751559./10546875./etac/etac
                       -128./5.*(1.-1/etac/etac)*gsl_sf_zeta(3))*CBM0(1) ) );
            
            cbme = kc * pow(model.Mrun(Muc, mc)/Muw, 2.) / 16. * (
                    
                   4.*M_PI/model.Als(Muc)*CBM0(1)
                    
                   + CBM1(1) + 4.*r1e*CBM1(0) + 15212./1875.*(1.-1./etac)*CBM0(1)
                    
                   + model.Als(Muc)/4./M_PI * (CBM2(1) + 4.*r2e*CBM0(0)
                   + 15212./1875.*(1.-1./etac)*CBM1(1)
                   + ( 966966391./10546875. - 231404944./3515625./etac - 272751559./10546875./etac/etac
                       -128./5.*(1.-1/etac/etac)*gsl_sf_zeta(3))*CBM0(1) ) );           
            */
            break;
        case(NLO):
            p1 = WMatchP(NLO);
            b1 = WMatchB(NLO);
            
            CPM0 = up4 * up5 * p0;
            CBM0 = ub4 * ub5 * b0;
            
            CPM1 = j4p1 * up4 * up5 * p0
                + etacb * up4 * (j5p1 -j4p1) * up5 *  p0
                + etab * etacb * up4 * up5 * (p1 - j5p1 * p0) ;
            CBM1 = j4b1 * ub4 * ub5 * b0
                + etacb * ub4 * (j5b1 -j4b1) * ub5 *  b0
                + etab * etacb * ub4 * ub5 * (b1 - j5b1 * b0) ;
            
            cpm = pow(mc/Muw, 2.)/32. * 4.*M_PI/A * (
                    
                    CPM0(2) +
                    
                    A/4./M_PI * (CPM1(2) + 4.*CPM0(0)*r1(0) + 4.*CPM0(1)*r1(1) ) );
            
            cbmt = pow(mc/Muw, 2.)/16. * 4.*M_PI/A * (
                    CBM0(1) +
                    
                    A/4./M_PI * (CBM1(1) + 4.*CBM0(0)*r1tau) );
            
            cbme = pow(mc/Muw, 2.)/16. * 4.*M_PI/A * (
                    CBM0(1) +
                    
                    A/4./M_PI * (CBM1(1) + 4.*CBM0(0)*r1e) );
            /*
            cpm = kc * pow(mc/Muw, 2.) / 32. *(
                    
                    4.*M_PI/model.Als(Muc) *CPM0(2) +
                    
                    CPM1(2) + 4.*CPM0(0)*(4.*(1-2.*log(Muc/
                    mc/Muw)) + 4.*log(kc)) +
                    4.*CPM1(0)*(-2.*(1-2.*log(Muc/
                    mc/Muw)) - 2.*log(kc)) +
                    15212./1875.*(1.-1/etac)*CPM0(2));
            
            cbmt = kc * pow(model.Mrun(Muc, mc)/Muw, 2.) / 16. * (
                    
                   4.*M_PI/model.Als(Muc)*CBM0(1)
                    
                   + CBM1(1) + 4.*r1tau*CBM1(0) + 15212./1875.*(1.-1./etac)*CBM0(1) );
            
            cbme = kc * pow(model.Mrun(Muc, mc)/Muw, 2.) / 16. * (
                    
                   4.*M_PI/model.Als(Muc)*CBM0(1)
                    
                   + CBM1(1) + 4.*r1e*CBM1(0) + 15212./1875.*(1.-1./etac)*CBM0(1) );
             */
            break;
        case(LO):
            CPM0 = up4 * up5 * p0;
            
            cpm = pow(mc/Muw, 2.)/32. * 4.*M_PI/A * CPM0(2);
            
            cbmt = pow(mc/Muw, 2.)/16. * 4.*M_PI/A * CBM0(1);
            
            cbme = pow(mc/Muw, 2.)/16. * 4.*M_PI/A * CBM0(1);
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("CPenguinBox::Cmatch() order " + out.str() + " not implemented");
    }
    
    return(3.*cpm + 2./3.*cbme + 1./3.*cbmt);
}

double CPenguinBox::CT_tot(orders order, orders_ew order_ew)
{
    
    double IBT = model.getDeltaP_cu();
    double xt = modelmatching.x_t(model.getMuw());
    double a = modelmatching.mt2omh2(model.getMuw());
    
    switch(order){
        case(NNLO):
            switch(order_ew){
                case(NLO_ew):
                    return(model.computelamc().real()/model.getLambda()*(Cmatch(NNLO) + IBT) +
                           model.computelamt().real()/model.getLambda()*(model.Als(model.getMuw())/4./M_PI*
                                 modelmatching.X1t(xt) + modelmatching.X1t(xt) + model.getAle()/4./M_PI*
                                 modelmatching.Xewt(xt, a, model.getMuc())));
                case(LO_ew):
                    return(model.computelamc().real()/model.getLambda()*(Cmatch(NNLO) + IBT) +
                           model.computelamt().real()/model.getLambda()*(model.Als(model.getMuw())/4./M_PI*
                                modelmatching.X1t(xt) + modelmatching.X1t(xt)));
                default:
                    std::stringstream out;
                    out << order_ew;
                    throw std::runtime_error("CPenguinBox::CT_tot() order_ew " + out.str() + " not implemented");
            }           
        case(NLO):
            switch(order_ew){
                case(NLO_ew):
                    return(model.computelamc().real()/model.getLambda()*(Cmatch(NLO) + IBT) +
                           model.computelamt().real()/model.getLambda()*(model.Als(model.getMuw())/4./M_PI*
                                 modelmatching.X1t(xt) + modelmatching.X1t(xt) + model.getAle()/4./M_PI*
                                 modelmatching.Xewt(xt, a, model.getMuc())));
                case(LO_ew):
                    return(model.computelamc().real()/model.getLambda()*(Cmatch(NLO) + IBT) +
                           model.computelamt().real()/model.getLambda()*(model.Als(model.getMuw())/4./M_PI*
                                modelmatching.X1t(xt) + modelmatching.X1t(xt)));
                default:
                    std::stringstream out;
                    out << order_ew;
                    throw std::runtime_error("CPenguinBox::CT_tot() order_ew " + out.str() + " not implemented");
            }
        case(LO):
            switch(order_ew){
                case(NLO_ew):
                    return(model.computelamc().real()/model.getLambda()*(Cmatch(LO) + IBT) +
                           model.computelamt().real()/model.getLambda()*(model.Als(model.getMuw())/4./M_PI*
                                 modelmatching.X1t(xt) + model.getAle()/4./M_PI*
                                 modelmatching.Xewt(xt, a, model.getMuc())));
                case(LO_ew):
                    return(model.computelamc().real()/model.getLambda()*(Cmatch(NLO) + IBT) +
                           model.computelamt().real()/model.getLambda()*(model.Als(model.getMuw())/4./M_PI*
                                modelmatching.X1t(xt)));
                default:
                    std::stringstream out;
                    out << order_ew;
                    throw std::runtime_error("CPenguinBox::CT_tot() order_ew " + out.str() + " not implemented");
            }
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("CPenguinBox::CT_tot() order " + out.str() + " not implemented");            
    }
}