/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RealWeakEFTLFVMatching.h"
#include "RealWeakEFTLFV.h"
#include <stdexcept>

RealWeakEFTLFVMatching::RealWeakEFTLFVMatching(const RealWeakEFTLFV & RealWeakEFTLFV_i) :

    StandardModelMatching(RealWeakEFTLFV_i),
    myRealWeakEFTLFV(RealWeakEFTLFV_i),
    mcbsg(8, NDR, NNLO),
    mcprimebsg(8, NDR, NNLO),
    mcbsmm(8, NDR, NNLO, NLO_QED22),
    mcbdmm(8, NDR, NNLO, NLO_QED22),
    mcBMll(13, NDR, NLO),
    mcprimeBMll(13, NDR, NLO)
{}

void RealWeakEFTLFVMatching::updateRealWeakEFTLFVParameters()
{    
    C7NP = myRealWeakEFTLFV.getC7();
    C7pNP = myRealWeakEFTLFV.getC7p();
    C8NP = myRealWeakEFTLFV.getC8();
    C8pNP = myRealWeakEFTLFV.getC8p();
    
    C9NPe = myRealWeakEFTLFV.getC9_11();
    C9pNPe  = myRealWeakEFTLFV.getC9p_11();
    C10NPe = myRealWeakEFTLFV.getC10_11();
    C10pNPe  = myRealWeakEFTLFV.getC10p_11();
    CSNPe = myRealWeakEFTLFV.getCS_11();
    CSpNPe = myRealWeakEFTLFV.getCSp_11();
    CPNPe = myRealWeakEFTLFV.getCP_11();
    CPpNPe = myRealWeakEFTLFV.getCPp_11();
    
    C9NPmu = myRealWeakEFTLFV.getC9_22();
    C9pNPmu  = myRealWeakEFTLFV.getC9p_22();
    C10NPmu = myRealWeakEFTLFV.getC10_22();
    C10pNPmu  = myRealWeakEFTLFV.getC10p_22();
    CSNPmu = myRealWeakEFTLFV.getCS_22();
    CSpNPmu = myRealWeakEFTLFV.getCSp_22();
    CPNPmu = myRealWeakEFTLFV.getCP_22();
    CPpNPmu = myRealWeakEFTLFV.getCPp_22();
    
    WCscale = myRealWeakEFTLFV.getWCscale();

    StandardModelMatching::updateSMParameters();
}

RealWeakEFTLFVMatching::~RealWeakEFTLFVMatching()
{}

std::vector<WilsonCoefficient>& RealWeakEFTLFVMatching::CMbsmm()
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

    mcbsmm.setMu(WCscale);

    switch (mcbsmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbsmm.setCoeff(7, 0., NLO);
        case LO:
            mcbsmm.setCoeff(7, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcbsmm.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsmm(): order " + out.str() + "not implemented"); 
    }

    vmcbsmm.push_back(mcbsmm);
    return (vmcbsmm);
}

std::vector<WilsonCoefficient>& RealWeakEFTLFVMatching::CMbdmm()
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

    mcbdmm.setMu(WCscale);

    switch (mcbdmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbdmm.setCoeff(7, 0., NLO);
        case LO:
            mcbdmm.setCoeff(7, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcbdmm.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbdmm(): order " + out.str() + "not implemented"); 
    }

    vmcbdmm.push_back(mcbdmm);
    return (vmcbdmm);
}

std::vector<WilsonCoefficient>& RealWeakEFTLFVMatching::CMbsg()
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

    mcbsg.setMu(WCscale);

    switch (mcbsg.getOrder()) {
        case NNLO:
            mcbsg.setCoeff(6, 0., NNLO);
        case NLO:
            mcbsg.setCoeff(6, 0., NLO);
        case LO:
            mcbsg.setCoeff(6, C7NP, LO);
            break;
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }

    vmcbsg.push_back(mcbsg);
    return (vmcbsg);
}
 
std::vector<WilsonCoefficient>& RealWeakEFTLFVMatching::CMprimebsg()
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

    mcprimebsg.setMu(WCscale);

    switch (mcprimebsg.getOrder()) {
        case NNLO:
            mcprimebsg.setCoeff(6, 0., NNLO);
        case NLO:
            mcprimebsg.setCoeff(6, 0., NLO);
        case LO:
            mcprimebsg.setCoeff(6, C7pNP, LO);
            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimebsg(): order " + out.str() + "not implemented"); 
    }

    vmcprimebsg.push_back(mcprimebsg);
    return (vmcprimebsg);
}

std::vector<WilsonCoefficient>& RealWeakEFTLFVMatching::CMBMll(QCD::lepton lepton)
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

    mcBMll.setMu(WCscale);

    switch (mcBMll.getOrder()) {
        case NNLO:
        case NLO:
                mcBMll.setCoeff(6, 0., NLO);
                mcBMll.setCoeff(8, 0., NLO);
                mcBMll.setCoeff(9, 0., NLO);
                mcBMll.setCoeff(10, 0., NLO);
                mcBMll.setCoeff(11, 0., NLO);
        case LO:
            mcBMll.setCoeff(6, C7NP, LO);
            if(lepton == RealWeakEFTLFV::ELECTRON){
                mcBMll.setCoeff(8, C9NPe, LO);
                mcBMll.setCoeff(9, C10NPe, LO);
                mcBMll.setCoeff(10, CSNPe, LO);
                mcBMll.setCoeff(11, CPNPe, LO);
            }
            else{
                mcBMll.setCoeff(8, C9NPmu, LO);
                mcBMll.setCoeff(9, C10NPmu, LO);
                mcBMll.setCoeff(10, CSNPmu, LO);
                mcBMll.setCoeff(11, CPNPmu, LO);               
            }
            break;
        default:
            std::stringstream out;
            out << mcBMll.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMBMll(): order " + out.str() + "not implemented"); 
    }

    vmcBMll.push_back(mcBMll);
    return (vmcBMll);
}

std::vector<WilsonCoefficient>& RealWeakEFTLFVMatching::CMprimeBMll(QCD::lepton lepton)
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

    mcprimeBMll.setMu(WCscale);

    switch (mcprimeBMll.getOrder()) {
        case NNLO:
        case NLO:
                mcprimeBMll.setCoeff(6, 0., NLO);
                mcprimeBMll.setCoeff(8, 0., NLO);
                mcprimeBMll.setCoeff(9, 0., NLO);
                mcprimeBMll.setCoeff(10, 0., NLO);
                mcprimeBMll.setCoeff(11, 0., NLO);
        case LO:
            mcprimeBMll.setCoeff(6, C7pNP, LO);
            if(lepton == RealWeakEFTLFV::ELECTRON){
                mcprimeBMll.setCoeff(8, C9pNPe, LO);
                mcprimeBMll.setCoeff(9, C10pNPe, LO);
                mcprimeBMll.setCoeff(10, CSpNPe, LO);
                mcprimeBMll.setCoeff(11, CPpNPe, LO);
            }
            else{
                mcprimeBMll.setCoeff(8, C9pNPmu, LO);
                mcprimeBMll.setCoeff(9, C10pNPmu, LO);
                mcprimeBMll.setCoeff(10, CSpNPmu, LO);
                mcprimeBMll.setCoeff(11, CPpNPmu, LO);              
            }
            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getOrder();
            throw std::runtime_error("FlavourWilsonCoefficientMatching::CMprimeBMll(): order " + out.str() + "not implemented"); 
    }

    vmcprimeBMll.push_back(mcprimeBMll);
    return (vmcprimeBMll);
}

