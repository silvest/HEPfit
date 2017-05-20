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
    mcbsg(8, NDR, NLO),
    mcprimebsg(8, NDR, NLO),
    mcBMll(13, NDR, NLO),
    mcprimeBMll(13, NDR, NLO)
{}

void LeftRightSymmetricModelMatching::updateLeftRightSymmetricModelParameters()
{
    myCKM = myLeftRightSymmetricModel.getVCKM();
    
    StandardModelMatching::updateSMParameters();
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
        case NLO:
            mcbsg.setCoeff(6, 0., NLO);
        case LO:
            mcbsg.setCoeff(6, 0., LO);
            break;
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