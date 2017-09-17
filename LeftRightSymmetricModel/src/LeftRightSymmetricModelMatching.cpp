/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LeftRightSymmetricModelMatching.h"
#include "LeftRightSymmetricModel.h"
#include <stdexcept>

LeftRightSymmetricModelMatching::LeftRightSymmetricModelMatching(const LeftRightSymmetricModel & LeftRightSymmetricModel_i) :

    StandardModelMatching(LeftRightSymmetricModel_i),
    myLeftRightSymmetricModel(LeftRightSymmetricModel_i),
    myCKM(3, 3, 0.),
    myCKMR(3, 3, 0.),
    mcbsg(8, NDR, NLO),
    mcprimebsg(8, NDR, NLO),
    mcBMll(13, NDR, NLO),
    mcprimeBMll(13, NDR, NLO)
{}

void LeftRightSymmetricModelMatching::updateLeftRightSymmetricModelParameters()
{
    myCKM = myLeftRightSymmetricModel.getVCKM();
    StandardModelMatching::updateSMParameters();
    Muw = myLeftRightSymmetricModel.getMuw();
    Mut = myLeftRightSymmetricModel.getMut();
    mW = myLeftRightSymmetricModel.Mw();
    mtop = myLeftRightSymmetricModel.Mrun(Mut, myLeftRightSymmetricModel.getQuarks(QCD::TOP).getMass_scale(),
                        myLeftRightSymmetricModel.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    mbottom = myLeftRightSymmetricModel.getQuarks(QCD::BOTTOM).getMass();
    vev = myLeftRightSymmetricModel.v();
    gW = sqrt(8.0*myLeftRightSymmetricModel.getGF() / sqrt(2.0)) * mW;
    myCKMR = myLeftRightSymmetricModel.getVCKMR();
    mWR = myLeftRightSymmetricModel.getmWR();
    mH2psq = myLeftRightSymmetricModel.getmH2p_2();
    xi = myLeftRightSymmetricModel.getxi_LRSM();
    alpha = myLeftRightSymmetricModel.getalpha_LRSM();
}

gslpp::complex LeftRightSymmetricModelMatching::setWCbsg(int i, double mu, orders order)
{
//    if ( tanbsg == tan && mtbsg == mt  && mhpbsg == mhp && mubsg == mu){
//        switch (order){
//        case NNLO:
//            return ( CWbsgArrayNNLO[i] );
//        case NLO:
//            return ( CWbsgArrayNLO[i] );
//            break;
//        case LO:
//            return ( CWbsgArrayLO[i] );
//            break;
//        default:
//            std::stringstream out;
//            out << order;
//            throw std::runtime_error("order" + out.str() + "not implemeted"); 
//        }
//    }
    
//    tanbsg = tan; mtbsg = mt; mhpbsg = mhp; mubsg = mu;
    updateLeftRightSymmetricModelParameters();
    double xt = mtop*mtop/(mW*mW);
    double yt = mtop*mtop/mH2psq;
    double xtm1sq = (xt-1.0)*(xt-1.0);
    double ytm1sq = (yt-1.0)*(yt-1.0);
    double kappa = vev*sqrt(0.5/(1.0+xi*xi));
    double kappaprime = xi*kappa;
    double sinb = kappaprime*sqrt(2.0)/vev;
    double kappaR = mWR/gW;
    gslpp::complex Atb = (mtop/mbottom)* kappaprime*kappa/(kappaR*kappaR)*(cos(alpha)+gslpp::complex::i()*sin(alpha))*(myCKMR(2,2)/myCKM(2,2));
    double A1Hp = ((yt*yt-2.0/3.0*yt)*log(yt)/((1.0-yt)*ytm1sq) + (5.0*yt*yt-3.0*yt)/(6.0*ytm1sq));
    double A2Hp = -((0.5*yt*yt*yt-yt*yt/3.0)*log(yt)/(ytm1sq*ytm1sq)+(-8.0*yt*yt*yt-5.0*yt*yt+7.0*yt)/(36.0*ytm1sq*(yt-1.0)))-A1Hp;
    gslpp::complex DeltaC7_LR = Atb * ((1.5*xt*xt-xt)*log(xt)/((1.0-xt)*xtm1sq) + (-5.0*xt*xt+31.0*xt-20.0)/(12.0*xtm1sq));
    gslpp::complex DeltaC7_Hp = -(sinb*kappa*sqrt(2.0)/vev * mtop/mbottom * (cos(alpha)+gslpp::complex::i()*sin(alpha))*(myCKMR(2,2)/myCKM(2,2))*A1Hp 
                            +4.0*sinb*sinb*kappa*kappa/(vev*vev)*A2Hp)
                          /(1.0-4.0*sinb*sinb+4.0*sinb*sinb*sinb*sinb);
    gslpp::complex DeltaC8_LR = Atb * (-1.5*xt*log(xt)/((1.0-xt)*xtm1sq) + (-xt*xt-xt-4.0)/(4.0*xtm1sq));
    gslpp::complex DeltaC8_Hp = 0.0;

    switch (order){
        case NNLO:
        case NLO:
        case LO:

            CWbsgArrayLO[6] = DeltaC7_LR + DeltaC7_Hp;
            CWbsgArrayLO[7] = DeltaC8_LR + DeltaC8_Hp;

            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
            }
    
    /*std::cout << "CWbsgArrayLO[6] = " << CWbsgArrayLO[6] << std::endl;
    std::cout << "CWbsgArrayLO[7] = " << CWbsgArrayLO[7] << std::endl << std::endl;*/
    
    
    switch (order){
//        case NNLO:
//        case NLO:
        case LO:
            return ( CWbsgArrayLO[i] );
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
        }
}

std::vector<WilsonCoefficient>& LeftRightSymmetricModelMatching::CMbsg()
{

    vmcbsg = StandardModelMatching::CMbsg();

    switch (mcbsg.getScheme()) {
        case NDR:
            break;
        default:
            std::stringstream out;
            out << mcbsg.getScheme();
            throw std::runtime_error("LeftRightSymmetricModelMatching::CMbsg(): scheme " + out.str() + "not implemented"); 
    }

    mcbsg.setMu(myLeftRightSymmetricModel.getMuw());

    switch (mcbsg.getOrder()) {
        case NNLO:
            mcbsg.setCoeff(6, 0., NNLO);
        case NLO:
            mcbsg.setCoeff(6, 0., NLO);
        case LO: {
                std::cout<<"LRSM matching"<<std::endl;
            gslpp::complex dummy7 = setWCbsg(6, 100.0, LO);
            gslpp::complex dummy8 = setWCbsg(7, 100.0, LO);
            std::cout<<"C7 = "<<dummy7<<std::endl;
            std::cout<<"C8 = "<<dummy8<<std::endl;
            mcbsg.setCoeff(6, 0., LO);
            break;
        }
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("LeftRightSymmetricModelMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }

    vmcbsg.push_back(mcbsg);
    return (vmcbsg);
}
 
std::vector<WilsonCoefficient>& LeftRightSymmetricModelMatching::CMprimebsg()
{

    vmcprimebsg = StandardModelMatching::CMprimebsg();

    switch (mcprimebsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getScheme();
            throw std::runtime_error("LeftRightSymmetricModelMatching::CMprimebsg(): scheme " + out.str() + "not implemented"); 
    }

    mcprimebsg.setMu(myLeftRightSymmetricModel.getMuw());

    switch (mcprimebsg.getOrder()) {
        case NNLO:
        case NLO:
            mcprimebsg.setCoeff(6, 0., NLO);
        case LO:
            mcprimebsg.setCoeff(6, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getOrder();
            throw std::runtime_error("LeftRightSymmetricModelMatching::CMprimebsg(): order " + out.str() + "not implemented"); 
    }

    vmcprimebsg.push_back(mcprimebsg);
    return (vmcprimebsg);
}

std::vector<WilsonCoefficient>& LeftRightSymmetricModelMatching::CMBMll(QCD::lepton lepton)
{

    vmcBMll = StandardModelMatching::CMBMll(lepton);

    switch (mcbsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcBMll.getScheme();
            throw std::runtime_error("LeftRightSymmetricModelMatching::CMBMll(): scheme " + out.str() + "not implemented"); 
    }

    mcBMll.setMu(myLeftRightSymmetricModel.getMuw());

    switch (mcBMll.getOrder()) {
        case NNLO:
        case NLO:
            //mcBMll.setCoeff(6, 0., NLO);
        case LO:
            mcBMll.setCoeff(6, 0., LO);
            mcBMll.setCoeff(8, 0., LO);
            mcBMll.setCoeff(9, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcBMll.getOrder();
            throw std::runtime_error("LeftRightSymmetricModelMatching::CMBMll(): order " + out.str() + "not implemented"); 
    }

    vmcBMll.push_back(mcBMll);
    return (vmcBMll);
}

std::vector<WilsonCoefficient>& LeftRightSymmetricModelMatching::CMprimeBMll(QCD::lepton lepton)
{

    vmcprimeBMll = StandardModelMatching::CMprimeBMll(lepton);

    switch (mcprimebsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getScheme();
            throw std::runtime_error("LeftRightSymmetricModelMatching::CMprimeBMll(): scheme " + out.str() + "not implemented"); 
    }

    mcprimeBMll.setMu(myLeftRightSymmetricModel.getMuw());

    switch (mcprimeBMll.getOrder()) {
        case NNLO:
        case NLO:
            //mcprimeBMll.setCoeff(6, 0., NLO);
        case LO:
            mcprimeBMll.setCoeff(6, 0., LO);
            mcprimeBMll.setCoeff(8, 0., LO);
            mcprimeBMll.setCoeff(9, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getOrder();
            throw std::runtime_error("LeftRightSymmetricModelMatching::CMprimeBMll(): order " + out.str() + "not implemented"); 
    }

    vmcprimeBMll.push_back(mcprimeBMll);
    return (vmcprimeBMll);
}
