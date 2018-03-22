/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "NPSMEFTd6Matching.h"
#include "NPSMEFTd6.h"
#include <stdexcept>

NPSMEFTd6Matching::NPSMEFTd6Matching(const NPSMEFTd6 & NPSMEFTd6_i) :

    StandardModelMatching(NPSMEFTd6_i),
    myNPSMEFTd6(NPSMEFTd6_i),
    myCKM(3, 3, 0.),
    mcbsg(8, NDR, NNLO),
    mcprimebsg(8, NDR, NNLO),
    mcbsmm(8, NDR, NNLO, NLO_QED22),
    mcbdmm(8, NDR, NNLO, NLO_QED22),
    mcBMll(13, NDR, NLO),
    mcprimeBMll(13, NDR, NLO)
{}

void NPSMEFTd6Matching::updateNPSMEFTd6Parameters()
{
    myCKM = myNPSMEFTd6.getVCKM();
    
    C9 = 0.;
    C10 = 0.;
    C7g = 0.;
    C8g = 0.;
    
    C9p = 0.;
    C10p = 0.;
    C7gp = 0.;
    C8gp = 0.;
    
    Muw = myNPSMEFTd6.getMuw();
    
    StandardModelMatching::updateSMParameters();
}

NPSMEFTd6Matching::~NPSMEFTd6Matching()
{}

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMbsmm()
{
    vmcbsmm.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMbsmm().begin(); it != StandardModelMatching::CMbsmm().end(); it++ ) vmcbsmm.push_back(*it);

    switch (mcbsmm.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcbsmm.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsmm(): scheme " + out.str() + "not implemented"); 
    }

    mcbsmm.setMu(Muw);

    switch (mcbsmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbsmm.setCoeff(7, 0., NLO);
        case LO:
            mcbsmm.setCoeff(7, C10, LO);
            break;
        default:
            std::stringstream out;
            out << mcbsmm.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsmm(): order " + out.str() + "not implemented"); 
    }

    vmcbsmm.push_back(mcbsmm);
    return (vmcbsmm);
}

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMbdmm()
{
    vmcbdmm.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMbdmm().begin(); it != StandardModelMatching::CMbdmm().end(); it++ ) vmcbdmm.push_back(*it);

    switch (mcbdmm.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcbdmm.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbdmm(): scheme " + out.str() + "not implemented"); 
    }

    mcbdmm.setMu(Muw);

    switch (mcbdmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbdmm.setCoeff(7, 0., NLO);
        case LO:
            mcbdmm.setCoeff(7, C10, LO);
            break;
        default:
            std::stringstream out;
            out << mcbdmm.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbdmm(): order " + out.str() + "not implemented"); 
    }

    vmcbdmm.push_back(mcbdmm);
    return (vmcbdmm);
}

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMbsg()
{
    vmcbsg.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMbsg().begin(); it != StandardModelMatching::CMbsg().end(); it++ ) vmcbsg.push_back(*it);

    switch (mcbsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcbsg.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsg(): scheme " + out.str() + "not implemented"); 
    }

    mcbsg.setMu(Muw);

    switch (mcbsg.getOrder()) {
        case NNLO:
            mcbsg.setCoeff(6, 0., NNLO);
        case NLO:
            mcbsg.setCoeff(6, 0., NLO);
        case LO:
            mcbsg.setCoeff(6, C7g, LO);
            break;
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }

    vmcbsg.push_back(mcbsg);
    return (vmcbsg);
}
 
std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMprimebsg()
{
    vmcprimebsg.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMprimebsg().begin(); it != StandardModelMatching::CMprimebsg().end(); it++ ) vmcprimebsg.push_back(*it);

    switch (mcprimebsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimebsg(): scheme " + out.str() + "not implemented"); 
    }

    mcprimebsg.setMu(Muw);

    switch (mcprimebsg.getOrder()) {
        case NNLO:
            mcprimebsg.setCoeff(6, 0., NNLO);
        case NLO:
            mcprimebsg.setCoeff(6, 0., NLO);
        case LO:
            mcprimebsg.setCoeff(6, C7gp, LO);
            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimebsg(): order " + out.str() + "not implemented"); 
    }

    vmcprimebsg.push_back(mcprimebsg);
    return (vmcprimebsg);
}

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMBMll(QCD::lepton lepton)
{
    vmcBMll.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMBMll(lepton).begin(); it != StandardModelMatching::CMBMll(lepton).end(); it++ ) vmcBMll.push_back(*it);

    switch (mcbsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcBMll.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMBMll(): scheme " + out.str() + "not implemented"); 
    }

    mcBMll.setMu(Muw);

    switch (mcBMll.getOrder()) {
        case NNLO:
        case NLO:
            //mcBMll.setCoeff(6, 0., NLO);
        case LO:
            mcBMll.setCoeff(6, C7g, LO);
            mcBMll.setCoeff(8, C9, LO);
            mcBMll.setCoeff(9, C10, LO);
            break;
        default:
            std::stringstream out;
            out << mcBMll.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMBMll(): order " + out.str() + "not implemented"); 
    }

    vmcBMll.push_back(mcBMll);
    return (vmcBMll);
}

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMprimeBMll(QCD::lepton lepton)
{
    vmcprimeBMll.clear();
    for (std::vector<WilsonCoefficient>::iterator it = StandardModelMatching::CMprimeBMll(lepton).begin(); it != StandardModelMatching::CMprimeBMll(lepton).end(); it++ ) vmcprimeBMll.push_back(*it);

    switch (mcprimebsg.getScheme()) {
        case NDR:

            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getScheme();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimeBMll(): scheme " + out.str() + "not implemented"); 
    }

    mcprimeBMll.setMu(Muw);

    switch (mcprimeBMll.getOrder()) {
        case NNLO:
        case NLO:
            //mcprimeBMll.setCoeff(6, 0., NLO);
        case LO:
            mcprimeBMll.setCoeff(6, C7gp, LO);
            mcprimeBMll.setCoeff(8, C9p, LO);
            mcprimeBMll.setCoeff(9, C10p, LO);
            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimeBMll(): order " + out.str() + "not implemented"); 
    }

    vmcprimeBMll.push_back(mcprimeBMll);
    return (vmcprimeBMll);
}