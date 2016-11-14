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
    mcbsg(8, NDR, NLO),
    mcprimebsg(8, NDR, NLO)
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
    double xt = x_t(Mut);
    double mHp2=myTHDM.getmHp2();
    double MW=myTHDM.Mw();
    double xHW=mHp2/(MW*MW);
    double xtH=xt/xHW;
    
    double tanb = myTHDM.gettanb();
    
    gslpp::complex co = 1.; // (- 4. * GF / sqrt(2)) * SM.computelamt_s(); THIS SHOULD ALREADY BE IMPLEMENTED IN THE OBSERVABLE 
    
    vmcbsg = StandardModelMatching::CMbsg();
    mcbsg.setMu(Muw);
    
    switch (mcbsg.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, co * myTHDM.Alstilde5(Muw) * setWCbsg(j, xtH, tanb, Muw * Muw / (xt * MW * MW), NLO), NLO);
            }
            //break;
        case LO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, co * setWCbsg(j, xtH, tanb, Muw * Muw / (xt * MW * MW), LO), LO);
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

double THDMMatching::setWCbsg(int i, double x, double tan, double mu, orders order)
{  
    if ( tanbsg == tan && xbsg == x && mubsg == mu){
        switch (order){
        case NNLO:
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
    
    tanbsg = tan; xbsg = x; mubsg = mu;
    
    double xm = 1. - x;
    double xm2 = xm*xm;
    double xm3 = xm2*xm;
    double xm4 = xm3*xm;
    double xm5 = xm4*xm;
    double x2 = x*x;
    double x3 = x2*x;
    double x4 = x3*x;
    //double x5 = x4*x;
    double Lx = log(x);
    double tan2 = tan*tan;
    double Li2 = gsl_sf_dilog(1.-1./x);
    
    double lstmu = log (mu);
    
    switch (order){
        case NNLO:
        case NLO:
            CWbsgArrayNLO[0] = 0.;
            CWbsgArrayNLO[3] = ( (-16.*x + 29.*x2 - 7.*x3)/(36.*xm3) + (-2.*x + 3.*x2)*Lx/(6.*xm4) ) / tan2;
            CWbsgArrayNLO[6] = ( ( (797.*x - 5436.*x2 + 7569.*x3 - 1202.*x4)/(486.*xm4) + 
                    (36.*x2 - 74.*x3 + 16.*x4)*Li2/(9.*xm4) + ((7.*x - 463.*x2 + 807.*x3 - 63.*x4)*Lx)/(81.*xm5) ) +
                    ( (-31.*x - 18.*x2 + 135.*x3 - 14.*x4)/(27.*xm4) + (-28.*x2 + 46.*x3 + 6.*x4)*Lx/(9.*xm5) ) * lstmu ) / tan2
                    - 4./9. * CWbsgArrayNLO[3] +
                    ( (28.*x - 52.*x2 + 8.*x3)/(3.*xm3) + (-48.*x + 112.*x2 - 32.*x3)*Li2/(9.*xm3) + 
                    (66.*x - 128.*x2 + 14.*x3)*Lx/(9.*xm4) ) +
                    ( (42.*x - 94.*x2 + 16.*x3)/(9.*xm3) + (32.*x - 56.*x2 - 12.*x3)*Lx/(9.*xm4) ) * lstmu ;
            CWbsgArrayNLO[7] = ( ( (1130.*x - 18153.*x2 + 7650.*x3 - 4451.*x4)/(1296.*xm4) + 
                    (30.*x2 - 17.*x3 + 13.*x4)*Li2/(6.*xm4) + (-2.*x - 2155.*x2 + 321.*x3 - 468.*x4)*Lx/(216.*xm5) ) +
                    ( (-38.*x - 261.*x2 + 18.*x3 - 7.*x4)/(36.*xm4) + (-31.*x2 - 17.*x3)*Lx/(6.*xm5) ) * lstmu ) / tan2
                    - 1./6. * CWbsgArrayNLO[3] +
                    ( (143.*x - 44.*x2 + 29.*x3)/(8.*xm3) + (-36.*x + 25.*x2 - 17.*x3)*Li2/(6.*xm3) + 
                    (165.*x - 7.*x2 + 34.*x3)*Lx/(12.*xm4) ) +
                    ( (81.*x - 16.*x2 + 7.*x3)/(6.*xm3) + (19.*x + 17.*x2)*Lx/(3.*xm4) ) * lstmu ;
        case LO:
            CWbsgArrayLO[0] = 0.;
            CWbsgArrayLO[6] = - ( (7. - 5.*x - 8.*x2)/(36.*xm3) + (2.*x - 3.*x2)*Lx/(6.*xm4) ) * x/2./tan2 
                    + ( (3.*x - 5.*x2)/(6.*xm2) + (2.*x - 3.*x2)*Lx/(3.*xm3) ) / 2. ;
            CWbsgArrayLO[7] = - ( (2. + 5.*x - x2)/(12.*xm3) + (x*Lx)/(2.*xm4) ) * x/2./tan2 
                    + ( (3.*x - x2)/(2.*xm2) + (x*Lx)/xm3 ) / 2. ;
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
            }
    
    /*std::cout << "CWbsgArrayLO[6] = " << CWbsgArrayLO[6] << std::endl;
    std::cout << "CWbsgArrayLO[7] = " << CWbsgArrayLO[7] << std::endl;
    std::cout << "CWbsgArrayNLO[3] = " << CWbsgArrayNLO[3] << std::endl;
    std::cout << "CWbsgArrayNLO[6] = " << CWbsgArrayNLO[6] << std::endl;
    std::cout << "CWbsgArrayNLO[7] = " << CWbsgArrayNLO[7] << std::endl << std::endl;*/
    
    switch (order){
        case NNLO:
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
