/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "FlavourWilsonCoefficientMatching.h"
#include "FlavourWilsonCoefficient.h"
#include <stdexcept>

FlavourWilsonCoefficientMatching::FlavourWilsonCoefficientMatching(const FlavourWilsonCoefficient & FlavourWilsonCoefficient_i) :

    StandardModelMatching(FlavourWilsonCoefficient_i),
    myFlavourWilsonCoefficient(FlavourWilsonCoefficient_i),
    myCKM(3, 3, 0.),
    mcbsg(8, NDR, NLO),
    mcprimebsg(8, NDR, NLO),
    mcbsmm(8, NDR, NNLO, NLO_ewt4),
    mcbdmm(8, NDR, NNLO, NLO_ewt4),
    mcBMll(13, NDR, NLO),
    mcprimeBMll(13, NDR, NLO)
{}

void FlavourWilsonCoefficientMatching::updateFlavourWilsonCoefficientParameters()
{
    myCKM = myFlavourWilsonCoefficient.getVCKM();
    
    DC1 = myFlavourWilsonCoefficient.getDC1();
    DC2 = myFlavourWilsonCoefficient.getDC2();
    DC3 = myFlavourWilsonCoefficient.getDC3();
    DC4 = myFlavourWilsonCoefficient.getDC4();
    DC5 = myFlavourWilsonCoefficient.getDC5();
    DC6 = myFlavourWilsonCoefficient.getDC6();
    DC7 = myFlavourWilsonCoefficient.getDC7();
    DC8 = myFlavourWilsonCoefficient.getDC8();
    DC9 = myFlavourWilsonCoefficient.getDC9();
    DC10 = myFlavourWilsonCoefficient.getDC10();
    DC7g = myFlavourWilsonCoefficient.getDC7g();
    DC8g = myFlavourWilsonCoefficient.getDC8g();
    
    DC7p = myFlavourWilsonCoefficient.getDC7p();
    DC8p = myFlavourWilsonCoefficient.getDC8p();
    DC9p = myFlavourWilsonCoefficient.getDC9p();
    DC10p = myFlavourWilsonCoefficient.getDC10p();
    DC7gp = myFlavourWilsonCoefficient.getDC7gp();
    DC8gp = myFlavourWilsonCoefficient.getDC8gp();
    
    StandardModelMatching::updateSMParameters();
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficientMatching::CMbsg()
{

    vmcbsg = StandardModelMatching::CMbsg();

    switch (mcbsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcbsg.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsg(): scheme " + out.str() + "not implemented"); 
    }

    mcbsg.setMu(myFlavourWilsonCoefficient.getWCscale());

    switch (mcbsg.getOrder()) {
        case NNLO:
        case NLO:
            mcbsg.setCoeff(6, 0., NLO);
        case LO:
            mcbsg.setCoeff(6, DC7g, LO);
            break;
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }

    vmcbsg.push_back(mcbsg);
    return (vmcbsg);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficientMatching::CMbsmm()
{

    vmcbsmm = StandardModelMatching::CMbsmm();

    switch (mcbsmm.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcbsmm.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsmm(): scheme " + out.str() + "not implemented"); 
    }

    mcbsmm.setMu(myFlavourWilsonCoefficient.getWCscale());

    switch (mcbsmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbsmm.setCoeff(7, 0., NLO);
        case LO:
            mcbsmm.setCoeff(7, DC10, LO);
            break;
        default:
            std::stringstream out;
            out << mcbsmm.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsmm(): order " + out.str() + "not implemented"); 
    }

    vmcbsmm.push_back(mcbsmm);
    return (vmcbsmm);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficientMatching::CMbdmm()
{

    vmcbdmm = StandardModelMatching::CMbdmm();

    switch (mcbdmm.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcbdmm.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbdmm(): scheme " + out.str() + "not implemented"); 
    }

    mcbdmm.setMu(myFlavourWilsonCoefficient.getWCscale());

    switch (mcbdmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbdmm.setCoeff(7, 0., NLO);
        case LO:
            mcbdmm.setCoeff(7, DC10, LO);
            break;
        default:
            std::stringstream out;
            out << mcbdmm.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbdmm(): order " + out.str() + "not implemented"); 
    }

    vmcbdmm.push_back(mcbdmm);
    return (vmcbdmm);
}
 
std::vector<WilsonCoefficient>& FlavourWilsonCoefficientMatching::CMprimebsg()
{

    vmcprimebsg = StandardModelMatching::CMprimebsg();

    switch (mcprimebsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimebsg(): scheme " + out.str() + "not implemented"); 
    }

    mcprimebsg.setMu(myFlavourWilsonCoefficient.getWCscale());

    switch (mcprimebsg.getOrder()) {
        case NNLO:
        case NLO:
            mcprimebsg.setCoeff(6, 0., NLO);
        case LO:
            mcprimebsg.setCoeff(6, DC7gp, LO);
            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimebsg(): order " + out.str() + "not implemented"); 
    }

    vmcprimebsg.push_back(mcprimebsg);
    return (vmcprimebsg);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficientMatching::CMBMll(QCD::lepton lepton)
{

    vmcBMll = StandardModelMatching::CMBMll(lepton);

    switch (mcbsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcBMll.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMBMll(): scheme " + out.str() + "not implemented"); 
    }

    mcBMll.setMu(myFlavourWilsonCoefficient.getWCscale());

    switch (mcBMll.getOrder()) {
        case NNLO:
        case NLO:
            //mcBMll.setCoeff(6, 0., NLO);
        case LO:
            mcBMll.setCoeff(6, DC7g, LO);
            mcBMll.setCoeff(8, DC9, LO);
            mcBMll.setCoeff(9, DC10, LO);
            break;
        default:
            std::stringstream out;
            out << mcBMll.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMBMll(): order " + out.str() + "not implemented"); 
    }

    vmcBMll.push_back(mcBMll);
    return (vmcBMll);
}

std::vector<WilsonCoefficient>& FlavourWilsonCoefficientMatching::CMprimeBMll(QCD::lepton lepton)
{

    vmcprimeBMll = StandardModelMatching::CMprimeBMll(lepton);

    switch (mcprimebsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimeBMll(): scheme " + out.str() + "not implemented"); 
    }

    mcprimeBMll.setMu(myFlavourWilsonCoefficient.getWCscale());

    switch (mcprimeBMll.getOrder()) {
        case NNLO:
        case NLO:
            //mcprimeBMll.setCoeff(6, 0., NLO);
        case LO:
            mcprimeBMll.setCoeff(6, DC7gp, LO);
            mcprimeBMll.setCoeff(8, DC9p, LO);
            mcprimeBMll.setCoeff(9, DC10p, LO);
            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimeBMll(): order " + out.str() + "not implemented"); 
    }

    vmcprimeBMll.push_back(mcprimeBMll);
    return (vmcprimeBMll);
}