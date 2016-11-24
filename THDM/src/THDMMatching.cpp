/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMMatching.h"
#include "THDM.h"
#include <gsl/gsl_sf_dilog.h>
#include <math.h>
#include <stdexcept>

THDMMatching::THDMMatching(const THDM & THDM_i) :

    StandardModelMatching(THDM_i),
    myTHDM(THDM_i),
    myCKM(3, 3, 0.),
    mcdbs2(5, NDR, NLO),
    mcbtaunu(3, NDR, LO),
    mcbsg(8, NDR, NNLO),
    mcprimebsg(8, NDR, NNLO)
{}

std::vector<WilsonCoefficient>& THDMMatching::CMdbs2() {

    double Mut = myTHDM.getMut();
    double xt = x_t(Mut);
    double GF=myTHDM.getGF();
    double MW=myTHDM.Mw();
    gslpp::complex co = GF / 4. / M_PI * MW * myTHDM.computelamt_s();
    double tanb = myTHDM.gettanb();
    double mHp2=myTHDM.getmHp2();
    double xHW=mHp2/(MW*MW);
    double xtH=xt/xHW;
    double SWH=xtH*((2.0*xHW-8.0)*log(xtH)/((1.0-xHW)*(1.0-xtH)*(1.0-xtH))+6.0*xHW*log(xt)/((1.0-xHW)*(1.0-xt)*(1.0-xt))-(8.0-2.0*xt)/((1.0-xt)*(1.0-xtH)))/(tanb*tanb);
    double SHH=xtH*((1.0+xtH)/((1.0-xtH)*(1.0-xtH))+2.0*xtH*log(xtH)/((1.0-xtH)*(1.0-xtH)*(1.0-xtH)))/(tanb*tanb*tanb*tanb);

    vmcds = StandardModelMatching::CMdbs2();
    mcdbs2.setMu(Mut);

    switch (mcdbs2.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcdbs2.setCoeff(0, co * co * xt * (SWH+SHH), LO);
            break;
        default:
            std::stringstream out;
            out << mcdbs2.getOrder();
            throw std::runtime_error("THDMMatching::CMdbs2(): order " + out.str() + "not implemented");
    }

    vmcds.push_back(mcdbs2);
    return(vmcds);
}

std::vector<WilsonCoefficient>& THDMMatching::CMbtaunu() {

    double Muw = myTHDM.getMuw();
    double GF = myTHDM.getGF();
    myCKM = myTHDM.getVCKM();
    double mB = myTHDM.getMesons(QCD::B_P).getMass();
    double tanb = myTHDM.gettanb();
    double mHp2=myTHDM.getmHp2();

    vmcbtaunu = StandardModelMatching::CMbtaunu();
    mcbtaunu.setMu(Muw);
 
    switch (mcbtaunu.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcbtaunu.setCoeff(0, -4.*GF * myCKM(0,2) / sqrt(2.) * mB*mB*tanb*tanb/mHp2, LO);
            break;
        default:
            std::stringstream out;
            out << mcbtaunu.getOrder();
            throw std::runtime_error("THDMMatching::CMbtaunu(): order " + out.str() + "not implemented");
    }

    vmcbtaunu.push_back(mcbtaunu);
    return(vmcbtaunu);

}

std::vector<WilsonCoefficient>& THDMMatching::CMbsg() 
{    
    double Muw = myTHDM.getMuw();
    
    double Mut = myTHDM.getMut();
    double mt = myTHDM.Mrun(Mut, myTHDM.getQuarks(QCD::TOP).getMass_scale(), 
                        myTHDM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    double mHp=myTHDM.getmHp();
    
    double tanb = myTHDM.gettanb();
    
    gslpp::complex co = 1.; // (- 4. * GF / sqrt(2)) * SM.computelamt_s(); THIS SHOULD ALREADY BE IMPLEMENTED IN THE OBSERVABLE 
    
    vmcbsg = StandardModelMatching::CMbsg();
    mcbsg.setMu(Muw);
    
    switch (mcbsg.getOrder()) {
        case NNLO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, co * myTHDM.Alstilde5(Muw) * myTHDM.Alstilde5(Muw) * setWCbsg(j, tanb, mt, mHp, Muw, NNLO), NNLO);
            }
        case NLO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, co * myTHDM.Alstilde5(Muw) * setWCbsg(j, tanb, mt, mHp, Muw, NLO), NLO);
            }
        case LO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, co * setWCbsg(j, tanb, mt, mHp, Muw, LO), LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("THDMMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }
    
    vmcbsg.push_back(mcbsg);
    return(vmcbsg);
}

 
 std::vector<WilsonCoefficient>& THDMMatching::CMprimebsg() 
{    
    double Muw = myTHDM.getMuw();
    
    vmcprimebsg = StandardModelMatching::CMprimebsg();
    mcprimebsg.setMu(Muw);
    
    switch (mcprimebsg.getOrder()) {
        case NNLO:
            for (int j=0; j<8; j++){
            mcprimebsg.setCoeff(j, 0., NNLO);//* CHECK ORDER *//
            }
        case NLO:
            for (int j=0; j<8; j++){
            mcprimebsg.setCoeff(j, 0., NLO);//* CHECK ORDER *//
            }
        case LO:
            for (int j=0; j<8; j++){
            mcprimebsg.setCoeff(j, 0., LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getOrder();
            throw std::runtime_error("THDMMatching::CMprimebsg(): order " + out.str() + "not implemented"); 
    }
    
    vmcprimebsg.push_back(mcprimebsg);
    return(vmcprimebsg);
}



/*******************************************************************************
 * Wilson coefficients calcoulus, misiak base for b -> s gamma                  *  
 * ****************************************************************************/

double THDMMatching::setWCbsg(int i, double tan, double mt, double mhp, double mu, orders order)
{  
    if ( tanbsg == tan && mtbsg == mt  && mhpbsg == mhp && mubsg == mu){
        switch (order){
        case NNLO:
            return ( CWbsgArrayNNLO[i] );
        case NLO:
            return ( CWbsgArrayNLO[i] );
            break;
        case LO:
            return ( CWbsgArrayLO[i] );
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
        }
    }
    
    tanbsg = tan; mtbsg = mt; mhpbsg = mhp; mubsg = mu;
    
    double x = mt*mt/mhp/mhp;
    
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    double x5 = x4*x;
    double xm = 1. - x;
    double xm2 = xm*xm;
    double xm3 = xm2*xm;
    double xm4 = xm3*xm;
    double xm6 = xm4*xm2;
    double xm8 = xm4*xm4;
    double xo = 1. - 1./x;
    double xo2 = xo*xo;
    double xo4 = xo*xo2*xo;
    double xo6 = xo4*xo2;
    double xo8 = xo4*xo4;
    double Lx = log(x);
    double Lx2 = Lx*Lx;
    double Lx3 = Lx2*Lx;
    double tan2 = tan*tan;
    double Li2 = gsl_sf_dilog(1.-1./x);
    
    double lstmu = 2. * log(mu/mt);
    
    double n70ct = - ( (7. - 5.*x - 8.*x2)/(36.*xm3) + (2.*x - 3.*x2)*Lx/(6.*xm4) ) * x/2.;
    double n70fr = ( (3.*x - 5.*x2)/(6.*xm2) + (2.*x - 3.*x2)*Lx/(3.*xm3) ) / 2.;
    
    double n80ct = - ( (2. + 5.*x - x2)/(12.*xm3) + (x*Lx)/(2.*xm4) ) * x/2.;
    double n80fr = ( (3.*x - x2)/(2.*xm2) + (x*Lx)/xm3 ) / 2.;
    
    double n41ct = (-16.*x + 29.*x2 - 7.*x3)/(36.*xm3) + (-2.*x + 3.*x2)*Lx/(6.*xm4);
    
    
    double n71ct = (797.*x - 5436.*x2 + 7569.*x3 - 1202.*x4)/(486.*xm4) + 
                    (36.*x2 - 74.*x3 + 16.*x4)*Li2/(9.*xm4) + 
                    ((7.*x - 463.*x2 + 807.*x3 - 63.*x4)*Lx)/(81.*xm3*xm2);
    double cd71ct = (-31.*x - 18.*x2 + 135.*x3 - 14.*x4)/(27.*xm4) + (-28.*x2 + 46.*x3 + 6.*x4)*Lx/(9.*xm3*xm2);
    double n71fr = (28.*x - 52.*x2 + 8.*x3)/(3.*xm3) + (-48.*x + 112.*x2 - 32.*x3)*Li2/(9.*xm3) + 
                    (66.*x - 128.*x2 + 14.*x3)*Lx/(9.*xm4);
    double cd71fr = (42.*x - 94.*x2 + 16.*x3)/(9.*xm3) + (32.*x - 56.*x2 - 12.*x3)*Lx/(9.*xm4);
    
    double n81ct = (1130.*x - 18153.*x2 + 7650.*x3 - 4451.*x4)/(1296.*xm4) + 
                    (30.*x2 - 17.*x3 + 13.*x4)*Li2/(6.*xm4) + 
                    (-2.*x - 2155.*x2 + 321.*x3 - 468.*x4)*Lx/(216.*xm3*xm2);
    double cd81ct = (-38.*x - 261.*x2 + 18.*x3 - 7.*x4)/(36.*xm4) + (-31.*x2 - 17.*x3)*Lx/(6.*xm3*xm2);
    double n81fr = (143.*x - 44.*x2 + 29.*x3)/(8.*xm3) + (-36.*x + 25.*x2 - 17.*x3)*Li2/(6.*xm3) + 
                    (165.*x - 7.*x2 + 34.*x3)*Lx/(12.*xm4);
    double cd81fr = (81.*x - 16.*x2 + 7.*x3)/(6.*xm3) + (19.*x + 17.*x2)*Lx/(3.*xm4);
    
    double n32ct = (10.*x4 + 30.*x2 - 20.*x)/(27.*xm4) * Li2 + 
                    (30.*x3 - 66.*x2 - 56.*x)/(81.*xm4) * Lx + (6.*x3 - 187.*x2 + 213.*x)/(-81.*xm3);
    double cd32ct = (-30.*x2 + 20.*x)/(27.*xm4)*Lx + (-35.*x3 + 145.*x2 - 80.*x)/(-81.*xm3);
    
    double n42ct = (515.*x4 - 906.*x3 + 99.*x2 + 182.*x)/(54.*xm4) * Li2 + 
                    (1030.*x4 - 2763.*x3 - 15.*x2 + 980.*x)/(-108.*xm3*xm2)*Lx + 
                    (-29467.*x4 + 68142.*x3 - 6717.*x2 - 18134.*x)/(1944.*xm4);
    double cd42ct = (-375.*x3 - 95.*x2 + 182.*x)/(-54.*xm3*xm2)*Lx + 
                    (133.*x4 - 108.*x3 + 4023.*x2 - 2320.*x)/(324.*xm4);
    
    double cd72ct = -(x * (67930.*x4 - 470095.*x3 + 1358478.*x2 - 700243.*x + 54970.))/(-2187.*xm3*xm2) + 
                    (x * (10422.*x4 - 84390.*x3 + 322801.*x2 - 146588.*x + 1435.))/(729.*xm6)*Lx +
                    (2.*x2 * (260.*x3 - 1515.*x2 + 3757.*x - 1446.))/(-27. * xm3*xm2) * Li2;
    double ce72ct = (x * (-518.*x4 + 3665.*x3 - 17397.*x2 + 3767.*x + 1843.))/(-162.*xm3*xm2) +
                    (x2 * (-63.*x3 + 532.*x2 + 2089.*x - 1118.))/(27.*xm6)*Lx;
    double cd72fr = -(x * (3790.*x3 - 22511.*x2 + 53614.*x - 21069.))/(81.*xm4) -
                    (2.*x * (-1266.*x3 + 7642.*x2 - 21467.*x + 8179.))/(-81.*xm3*xm2)*Lx +
                    (8.*x * (139.*x3 - 612.*x2 + 1103.*x - 342.))/(27.*xm4) * Li2;
    double ce72fr = -(x * (284.*x3 - 1435.*x2 + 4304.*x - 1425.))/(27.*xm4) -
                    (2.*x * (63.*x3 - 397.*x2 - 970.*x + 440.))/(-27.*xm3*xm2)*Lx;

    double cd82ct = -(x * (522347.*x4 - 2423255.*x3 + 2706021.*x2 - 5930609.*x + 148856))/(-11664.*xm3*xm2) + 
                    (x * (51948.*x4 - 233781.*x3 + 48634.*x2 - 698693.*x + 2452.))/(1944.*xm6)*Lx + 
                    (x2 * (481.*x3 - 1950.*x2 + 1523.*x - 2550.))/(-18.*xm3*xm2) * Li2;
    double ce82ct = (x * (-259.*x4 + 1117.*x3 + 2925.*x2 + 28411.*x + 2366.))/(-216.*xm3*xm2) -
                    (x2 * (139.*x2 + 2938.*x + 2683.))/(36.*xm6)*Lx;
    double cd82fr = -(x * (1463.*x3 - 5794.*x2 + 5543.*x - 15036.))/(27.*xm4) - 
                    (x * (-1887.*x3 + 7115.*x2 + 2519.*x + 19901.))/(-54.*xm3*xm2)*Lx -
                    (x * (-629.*x3 + 2178.*x2 - 1729.*x + 2196.))/(18.*xm4) * Li2;
    double ce82fr = -(x * (259.*x3 - 947.*x2 - 251.*x - 5973.))/(36.*xm4)-
                    (x * (139.*x2 + 2134.*x + 1183.))/(-18.*xm3*xm2)*Lx;
    
    double n72ct = 0.;
    if (mhp < 50.) 
        n72ct = -274.2/x5 - 72.13*Lx/x5 + 24.41/x4 - 168.3*Lx/x4 + 79.15/x3 - 
                103.8*Lx/x3 + 47.09/x2 - 38.12*Lx/x2 + 15.35/x - 8.753*Lx/x + 3.970;
    else if (mhp < mt)
        n72ct = 1.283 + 0.7158 * xo + 0.4119 * xo2 + 0.2629 * xo*xo2 + 0.1825 * xo4 + 
                0.1347 * xo*xo4 + 0.1040 * xo6 + 0.0831 * xo*xo6 + 0.06804 * xo8 + 
                0.05688 * xo*xo8 + 0.04833 * xo2*xo8 + 0.04163 * xo*xo2*xo8 + 0.03625 * xo4*xo8 + 
                0.03188 * xo*xo4*xo8 + 0.02827 * xo6*xo8 + 0.02525 * xo*xo6*xo8 + 0.02269 * xo8*xo8;
    else if (mhp < 520.)
        n72ct = 1.283 - 0.7158 * xm - 0.3039 * xm2 - 0.1549 * xm3 - 0.08625 * xm4 - 
                0.05020 * xm3*xm2 - 0.02970 * xm6 - 0.01740 * xm3*xm4 - 0.009752 * xm8 - 
                0.004877 * xm3*xm6 - 0.001721 * xm2*xm8 + 0.0003378 * xm3*xm8 + 0.001679 * xm4*xm8 + 
                0.002542 * xm*xm4*xm8 + 0.003083 * xm6*xm8 + 0.003404 * xm*xm6*xm8 + 0.003574 * xm8*xm8;
    else n72ct = -823.0*x5 + 42.30*x5*Lx3 - 412.4*x5*Lx2 - 3362*x5*Lx -
                1492*x4 - 23.26*x4*Lx3 - 541.4*x4*Lx2 - 2540*x4*Lx -
                1158*x3 - 34.50*x3*Lx3 - 348.2*x3*Lx2 - 1292*x3*Lx -
                480.9*x2 - 20.73*x2*Lx3 - 112.4*x2*Lx2 - 396.1*x2*Lx -
                8.278*x + 0.9225*x*Lx2 + 4.317*x*Lx;
    
    double n72fr = 0.;
    if (mhp < 50.) 
        n72fr = -( 194.3/x5 + 101.1*Lx/x5 - 24.97/x4 + 168.4*Lx/x4 - 78.90/x3 + 
                106.2*Lx/x3 - 49.32/x2 + 38.43*Lx/x2 - 12.91/x + 9.757*Lx/x + 8.088 );
    else if (mhp < mt)
        n72fr = -( 12.82 - 1.663 * xo - 0.8852 * xo2 - 0.4827 * xo*xo2 - 0.2976 * xo4 - 
                0.2021 * xo*xo4 - 0.1470 * xo6 - 0.1125 * xo*xo6 - 0.08931 * xo8 - 
                0.07291 * xo*xo8 - 0.06083 * xo2*xo8 - 0.05164 * xo*xo2*xo8 - 0.04446 * xo4*xo8 - 
                0.03873 * xo*xo4*xo8 - 0.03407 * xo6*xo8 - 0.03023 * xo*xo6*xo8 - 0.02702 * xo8*xo8 );
    else if (mhp < 400.)
        n72fr = -( 12.82 + 1.663 * xm + 0.7780 * xm2 + 0.3755 * xm3 + 0.1581 * xm4 + 
                0.03021 * xm3*xm2 - 0.04868 * xm6 - 0.09864 * xm3*xm4 - 0.1306 * xm8 - 
                0.1510 * xm3*xm6 - 0.1637 * xm2*xm8 - 0.1712 * xm3*xm8 - 0.1751 * xm4*xm8 - 
                0.1766 * xm*xm4*xm8 - 0.1763 * xm6*xm8 - 0.1748 * xm*xm6*xm8 - 0.1724 * xm8*xm8 );
    else n72fr = -( 2828 * x5 - 66.63 * x5*Lx3 + 469.4 * x5*Lx2 + 1986 * x5*Lx + 
                1480 * x4 + 36.08 * x4*Lx3 + 323.2 * x4*Lx2 + 169.9 * x4*Lx + 
                166.7 * x3 + 19.73 * x3*Lx3 - 46.61 * x3*Lx2 - 826.2 * x3*Lx - 
                524.1 * x2 - 8.889 * x2*Lx3 - 195.7 * x2*Lx2 - 870.3 * x2*Lx - 
                572.2 * x - 20.94 * x*Lx3 - 123.5 * x*Lx2 - 453.5 * x*Lx );
    
    double n82ct = 0.;
    if (mhp < 50.) 
        n82ct = 826.2/x5 - 300.7*Lx/x5 + 96.35/x4 + 91.89*Lx/x4 - 66.39/x3 +
                78.58*Lx/x3 - 39.76/x2 + 20.02*Lx/x2 - 5.214/x + 2.278;
    else if (mhp < mt)
        n82ct = 1.188 + 0.4078 * xo + 0.2002 * xo2 + 0.1190 * xo*xo2 + 0.07861 * xo4 + 
                0.05531 * xo*xo4 + 0.04061 * xo6 + 0.03075 * xo*xo6 + 0.02386 * xo8 + 
                0.01888 * xo*xo8 + 0.01520 * xo2*xo8 + 0.01241 * xo*xo2*xo8 + 0.01026 * xo4*xo8 + 
                0.008575 * xo*xo4*xo8 + 0.007238 * xo6*xo8 + 0.006164 * xo*xo6*xo8 + 0.005290 * xo8*xo8;
    else if (mhp < 600.)
        n82ct = 1.188 - 0.4078 * xm - 0.2076 * xm2 - 0.1265 * xm3 - 0.08570 * xm4 -
                0.06204 * xm3*xm2 - 0.04689 * xm6 - 0.03652 * xm3*xm4 - 0.02907 * xm8 -
                0.02354 * xm3*xm6 - 0.01933 * xm2*xm8 - 0.01605 * xm3*xm8 - 0.01345 * xm4*xm8 -
                0.01137 * xm*xm4*xm8 - 0.009678 * xm6*xm8 - 0.008293 * xm*xm6*xm8 - 0.007148 * xm8*xm8;
    else n82ct = -19606 * x5 - 226.7 * x5*Lx3 - 5251 * x5*Lx2 - 26090 * x5*Lx - 
                9016 * x4 - 143.4 * x4*Lx3 - 2244 * x4*Lx2 - 10102 * x4*Lx - 
                3357 * x3 - 66.32 * x3*Lx3 - 779.6 * x3*Lx2 - 3077 * x3*Lx - 
                805.5 * x2 - 22.98 * x2*Lx3 - 169.1 * x2*Lx2 - 602.7 * x2*Lx + 
                0.7437  * x + 0.6908 * x*Lx2 + 3.238 * x*Lx;
    
    double n82fr = 0.;
    if (mhp < 50.) 
        n82fr = -( -1003/x5 + 476.9*Lx/x5 - 205.7/x4 - 71.62*Lx/x4 + 62.26/x3 -
                110.7*Lx/x3 + 63.74/x2 - 35.42*Lx/x2 + 10.89/x - 3.174 );
    else if (mhp < mt)
        n82fr = -( -0.6110 - 1.095 * xo - 0.4463 * xo2 - 0.2568 * xo*xo2 - 0.1698 * xo4 - 
                0.1197 * xo*xo4 - 0.08761 * xo6 - 0.06595 * xo*xo6 - 0.05079 * xo8 - 
                0.03987 * xo*xo8 - 0.03182 * xo2*xo8 - 0.02577 * xo*xo2*xo8 - 0.02114 * xo4*xo8 - 
                0.01754 * xo*xo4*xo8 - 0.01471 * xo6*xo8 - 0.01244 * xo*xo6*xo8 - 0.01062 * xo8*xo8 );
    else if (mhp < 520.)
        n82fr = -( -0.6110 + 1.095 * xm + 0.6492 * xm2 + 0.4596 * xm3 + 0.3569 * xm4 +
                0.2910 * xm3*xm2 + 0.2438 * xm6 + 0.2075 * xm3*xm4 + 0.1785 * xm8 +
                0.1546 * xm3*xm6 + 0.1347 * xm2*xm8 + 0.1177 * xm3*xm8 + 0.1032 * xm4*xm8 +
                0.09073 * xm*xm4*xm8 + 0.07987 * xm6*xm8 + 0.07040 * xm*xm6*xm8 + 0.06210 * xm8*xm8 );
    else n82fr = -( -15961 * x5 + 1003 * x5*Lx3 - 2627 * x5*Lx2 - 29962 * x5*Lx - 
                11683 * x4 + 54.66 * x4*Lx3 - 2777 * x4*Lx2 - 17770 * x4*Lx - 
                6481 * x3 - 40.68 * x3*Lx3 - 1439 * x3*Lx2 - 7906 * x3*Lx - 
                2943 * x2 - 31.83 * x2*Lx3 - 612.6 * x2*Lx2 - 2770 * x2*Lx -
                929.8 * x - 19.80 * x*Lx3 - 174.7 * x*Lx2 - 658.4 * x*Lx );
        
    
    switch (order){
        case NNLO:
            CWbsgArrayNNLO[2] = ( n32ct + cd32ct * lstmu ) / tan2;
            CWbsgArrayNNLO[3] = ( n42ct + cd42ct * lstmu ) / tan2;
            CWbsgArrayNNLO[4] = 2./15. * n41ct / tan2 - CWbsgArrayNNLO[2] / 10.;
            CWbsgArrayNNLO[5] = 1./4. * n41ct / tan2 - CWbsgArrayNNLO[2] * 3./16.;
            CWbsgArrayNNLO[6] = ( n72ct + cd72ct * lstmu + ce72ct * lstmu * lstmu)/tan2 +
                                n72fr + cd72fr * lstmu + ce72fr * lstmu * lstmu
                     - 1./3.*CWbsgArrayNNLO[2] - 4./9.*CWbsgArrayNNLO[3] - 20./3.*CWbsgArrayNNLO[4] - 80./9.*CWbsgArrayNNLO[5];
            CWbsgArrayNNLO[7] = ( n82ct + cd82ct * lstmu + ce82ct * lstmu * lstmu)/tan2 + 
                                n82fr + cd82fr * lstmu + ce82fr * lstmu * lstmu
                    + CWbsgArrayNNLO[2] - 1./6.*CWbsgArrayNNLO[3] - 20.*CWbsgArrayNNLO[4] - 10./3.*CWbsgArrayNNLO[5];
        case NLO:
            CWbsgArrayNLO[3] = n41ct / tan2;
            CWbsgArrayNLO[6] = ( n71ct + cd71ct * lstmu ) / tan2 + n71fr + cd71fr * lstmu  - 4./9.*CWbsgArrayNLO[3];
            CWbsgArrayNLO[7] = ( n81ct + cd81ct * lstmu ) / tan2 + n81fr + cd81fr * lstmu  - 1./6.*CWbsgArrayNLO[3];
        case LO:
            CWbsgArrayLO[6] = n70ct /tan2 + n70fr ;
            CWbsgArrayLO[7] = n80ct /tan2 + n80fr ;
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
            }
    
    /*std::cout << "CWbsgArrayLO[6] = " << CWbsgArrayLO[6] << std::endl;
    std::cout << "CWbsgArrayLO[7] = " << CWbsgArrayLO[7] << std::endl << std::endl;
    std::cout << "CWbsgArrayNLO[3] = " << CWbsgArrayNLO[3] << std::endl;
    std::cout << "CWbsgArrayNLO[6] = " << CWbsgArrayNLO[6] << std::endl;
    std::cout << "CWbsgArrayNLO[7] = " << CWbsgArrayNLO[7] << std::endl << std::endl;
    std::cout << "CWbsgArrayNNLO[2] = " << CWbsgArrayNNLO[2] << std::endl;
    std::cout << "CWbsgArrayNNLO[3] = " << CWbsgArrayNNLO[3] << std::endl;
    std::cout << "CWbsgArrayNNLO[4] = " << CWbsgArrayNNLO[4] << std::endl;
    std::cout << "CWbsgArrayNNLO[5] = " << CWbsgArrayNNLO[5] << std::endl;
    std::cout << "CWbsgArrayNNLO[6] = " << CWbsgArrayNNLO[6] << std::endl;
    std::cout << "CWbsgArrayNNLO[7] = " << CWbsgArrayNNLO[7] << std::endl << std::endl;*/
    
    
    switch (order){
        case NNLO:
            return ( CWbsgArrayNNLO[i] );
        case NLO:
            return ( CWbsgArrayNLO[i] );
            break;
        case LO:
            return ( CWbsgArrayLO[i] );
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
        }
}
