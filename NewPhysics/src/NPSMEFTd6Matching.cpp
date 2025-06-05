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
    mcbsg(8, NDR, NNLO),
    mcprimebsg(8, NDR, NNLO),
    mcbsmm(8, NDR, NNLO, NLO_QED22),
    mcbdmm(8, NDR, NNLO, NLO_QED22),
    mcBMll(13, NDR, NLO),
    mcprimeBMll(13, NDR, NLO),
    mcbsnn(2, NDR, NLO, NLO_QED11)
{}

void NPSMEFTd6Matching::updateNPSMEFTd6Parameters()
{
    
    SMEFT_factor = (M_PI/myNPSMEFTd6.getAle())*(myNPSMEFTd6.v()/myNPSMEFTd6.getLambda_NP())*(myNPSMEFTd6.v()/myNPSMEFTd6.getLambda_NP())/myNPSMEFTd6.getCKM().computelamt_s();
    Muw = myNPSMEFTd6.getMuw();
    
    double ytop = myNPSMEFTd6.getQuarks(QCD::TOP).getMass()/myNPSMEFTd6.v();
    loop_factor = myNPSMEFTd6.getCKM().computelamt_s()*ytop*ytop/(16.*M_PI*M_PI);
    
    C7NP = 0.; // to be implemented
    C8NP = 0.; // to be implemented 
    C7pNP = 0.; // to be implemented
    C8pNP = 0.; // to be implemented
    
    double logLambdaomu = log(myNPSMEFTd6.getLambda_NP()/Muw);
    
    CLQ1_1123_tot = myNPSMEFTd6.getCLQ1_1123();
    CLQ1_1123_tot -= loop_factor*logLambdaomu*(myNPSMEFTd6.getCHL1_11()-myNPSMEFTd6.getCLu_1133());
    CLQ1_2223_tot = myNPSMEFTd6.getCLQ1_2223();
    CLQ1_2223_tot -= loop_factor*logLambdaomu*(myNPSMEFTd6.getCHL1_22()-myNPSMEFTd6.getCLu_2233());
    CLQ1_3323_tot = myNPSMEFTd6.getCLQ1_3323();
    //CLQ1_3323_tot -= loop_factor*logLambdaomu*(myNPSMEFTd6.getCHL1_33()-myNPSMEFTd6.getCLu_3333());
    CLQ3_1123_tot = myNPSMEFTd6.getCLQ3_1123();
    CLQ3_1123_tot += loop_factor*logLambdaomu*myNPSMEFTd6.getCHL3_11();
    CLQ3_2223_tot = myNPSMEFTd6.getCLQ3_2223();
    CLQ3_2223_tot += loop_factor*logLambdaomu*myNPSMEFTd6.getCHL3_22();
    CQe_2311_tot = myNPSMEFTd6.getCQe_2311();
    CQe_2311_tot -= loop_factor*logLambdaomu*(myNPSMEFTd6.getCHe_11()-myNPSMEFTd6.getCeu_1133());
    CQe_2322_tot = myNPSMEFTd6.getCQe_2322();
    CQe_2322_tot -= loop_factor*logLambdaomu*(myNPSMEFTd6.getCHe_22()-myNPSMEFTd6.getCeu_2233());
    
    C9NPmu = SMEFT_factor*(CQe_2322_tot+CLQ1_2223_tot+CLQ3_2223_tot);
    C10NPmu = SMEFT_factor*(CQe_2322_tot-CLQ1_2223_tot-CLQ3_2223_tot);
    C9pNPmu = SMEFT_factor*(myNPSMEFTd6.getCed_2223()+myNPSMEFTd6.getCLd_2223());
    C10pNPmu = SMEFT_factor*(myNPSMEFTd6.getCed_2223()-myNPSMEFTd6.getCLd_2223());
    CSNPmu = SMEFT_factor*myNPSMEFTd6.getCLedQ_22();
    CSpNPmu = SMEFT_factor*myNPSMEFTd6.getCpLedQ_22();
    CPNPmu = -SMEFT_factor*CSNPmu;
    CPpNPmu = SMEFT_factor*CSpNPmu;
    
    C9NPe = SMEFT_factor*(CQe_2311_tot+CLQ1_1123_tot+CLQ3_1123_tot);
    C10NPe = SMEFT_factor*(CQe_2311_tot-CLQ1_1123_tot-CLQ3_1123_tot);
    C9pNPe  = SMEFT_factor*(myNPSMEFTd6.getCed_1123()+myNPSMEFTd6.getCLd_1123());
    C10pNPe  = SMEFT_factor*(myNPSMEFTd6.getCed_1123()-myNPSMEFTd6.getCLd_1123());
    CSNPe = SMEFT_factor*myNPSMEFTd6.getCLedQ_11();
    CSpNPe = SMEFT_factor*myNPSMEFTd6.getCpLedQ_11();
    CPNPe = -SMEFT_factor*CSNPe;
    CPpNPe = SMEFT_factor*CSpNPe; 

    CLnunuNPe = SMEFT_factor*(CLQ1_1123_tot-CLQ3_1123_tot);
    CRnunuNPe = SMEFT_factor*myNPSMEFTd6.getCLd_1123();
    CLnunuNPmu = SMEFT_factor*(CLQ1_2223_tot-CLQ3_2223_tot);
    CRnunuNPmu = SMEFT_factor*myNPSMEFTd6.getCLd_2223();
    CLnunuNPtau = SMEFT_factor*(CLQ1_3323_tot-CLQ3_3323_tot);
    CRnunuNPtau = SMEFT_factor*myNPSMEFTd6.getCLd_3323();
    
    StandardModelMatching::updateSMParameters();
}

NPSMEFTd6Matching::~NPSMEFTd6Matching()
{}

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMbsmm()
{
    vmcbsmm.clear();
    for (std::vector<WilsonCoefficient>::const_iterator it = StandardModelMatching::CMbsmm().begin(); it != StandardModelMatching::CMbsmm().end(); it++ ) vmcbsmm.push_back(*it);

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

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMbdmm()
{
    vmcbdmm.clear();
    for (std::vector<WilsonCoefficient>::const_iterator it = StandardModelMatching::CMbdmm().begin(); it != StandardModelMatching::CMbdmm().end(); it++ ) vmcbdmm.push_back(*it);

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

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMbsg()
{
    vmcbsg.clear();
    for (std::vector<WilsonCoefficient>::const_iterator it = StandardModelMatching::CMbsg().begin(); it != StandardModelMatching::CMbsg().end(); it++ ) vmcbsg.push_back(*it);

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
 
std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMprimebsg()
{
    vmcprimebsg.clear();
    for (std::vector<WilsonCoefficient>::const_iterator it = StandardModelMatching::CMprimebsg().begin(); it != StandardModelMatching::CMprimebsg().end(); it++ ) vmcprimebsg.push_back(*it);

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

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMBMll(QCD::lepton lepton)
{
    vmcBMll.clear();
    for (std::vector<WilsonCoefficient>::const_iterator it = StandardModelMatching::CMBMll(lepton).begin(); it != StandardModelMatching::CMBMll(lepton).end(); it++ ) vmcBMll.push_back(*it);

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
                mcBMll.setCoeff(6, 0., NLO);
                mcBMll.setCoeff(8, 0., NLO);
                mcBMll.setCoeff(9, 0., NLO);
                mcBMll.setCoeff(10, 0., NLO);
                mcBMll.setCoeff(11, 0., NLO);
        case LO:
            mcBMll.setCoeff(6, C7NP, LO);
            if(lepton == NPSMEFTd6::ELECTRON){
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

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMprimeBMll(QCD::lepton lepton)
{
    vmcprimeBMll.clear();
    for (std::vector<WilsonCoefficient>::const_iterator it = StandardModelMatching::CMprimeBMll(lepton).begin(); it != StandardModelMatching::CMprimeBMll(lepton).end(); it++ ) vmcprimeBMll.push_back(*it);

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
                mcprimeBMll.setCoeff(6, 0., NLO);
                mcprimeBMll.setCoeff(8, 0., NLO);
                mcprimeBMll.setCoeff(9, 0., NLO);
                mcprimeBMll.setCoeff(10, 0., NLO);
                mcprimeBMll.setCoeff(11, 0., NLO);
        case LO:
            mcprimeBMll.setCoeff(6, C7pNP, LO);
            if(lepton == NPSMEFTd6::ELECTRON){
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

std::vector<WilsonCoefficient>& NPSMEFTd6Matching::CMBXsnn(QCD::lepton lepton) {

    vmcbsnn.clear();
    for (std::vector<WilsonCoefficient>::const_iterator it = StandardModelMatching::CMBXsnn(lepton).begin(); it != StandardModelMatching::CMBXsnn(lepton).end(); it++ ) vmcbsnn.push_back(*it);

    mcbsnn.setMu(Muw);

    switch (mcbsnn.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff(0, 0., NLO);
            mcbsnn.setCoeff(1, 0., NLO);
        case LO:
            if(lepton == NPSMEFTd6::NEUTRINO_1){
                mcbsnn.setCoeff(0, CLnunuNPe, LO);
                mcbsnn.setCoeff(1, CRnunuNPe, LO);
            }
            else if(lepton == NPSMEFTd6::NEUTRINO_2){
                mcbsnn.setCoeff(0, CLnunuNPmu, LO);
                mcbsnn.setCoeff(1, CRnunuNPmu, LO);
            }
            else if(lepton == NPSMEFTd6::NEUTRINO_3){
                mcbsnn.setCoeff(0, CLnunuNPtau, LO);
                mcbsnn.setCoeff(1, CRnunuNPtau, LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcbsnn.getOrder();
            throw std::runtime_error("NPSMEFTd6Matching::CMBXsnn(): order " + out.str() + "not implemented"); 
    }

    switch (mcbsnn.getOrder_qed()) {
        case NLO_QED11:
            mcbsnn.setCoeff(0, 0., NLO_QED11);
            mcbsnn.setCoeff(1, 0., NLO_QED11);
        case LO_QED:
            mcbsnn.setCoeff(0, 0., LO_QED);
            mcbsnn.setCoeff(1, 0., LO_QED);
            break; 
        default:
            std::stringstream out;
            out << mcbsnn.getOrder_qed();
            throw std::runtime_error("NPSMEFTd6Matching::CMXsnn(): qed order " + out.str() + " not implemented"); 
    }

    vmcbsnn.push_back(mcbsnn);
    return (vmcbsnn);
}