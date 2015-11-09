/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDLi3j.h"
#include "gslpp_complex.h"

using namespace gslpp;

HeffDLi3j::HeffDLi3j(const StandardModel & SM_i) :
        model(SM_i),
        coeffDLi3j_1(12, NDR , LO),coeffDLi3j_2(12, NDR , LO),coeffDLi3j_3(12, NDR , LO){
}

HeffDLi3j::~HeffDLi3j() {
}

vector<complex>** HeffDLi3j::ComputeCoeffDLi3j(int li_lj) {

//    switch (li_lj) {
//        case 1:
//        {
            std::vector<WilsonCoefficient>& mcb4 = model.getMyMatching() -> CMDLi3j(1);
            orders ordDLi3j_1 = coeffDLi3j_1.getOrder();
            for (unsigned int i = 0; i < mcb4.size(); i++){
                for (int j = LO; j <= ordDLi3j_1; j++){
                    coeffDLi3j_1.setCoeff(*coeffDLi3j_1.getCoeff(orders(j))
                                            + *mcb4[i].getCoeff(orders(j)), orders(j));
                }
            }
        std::cout<<"mcb4 = "<<**coeffDLi3j_1.getCoeff()<<std::endl;
        return coeffDLi3j_1.getCoeff();
//        }
//        case 2:
//        {
//            std::vector<WilsonCoefficient>& mcb2 = model.getMyMatching() -> CMDLijjj(2);
//            orders ordDL1_2 = coeffDL1_2.getOrder();
//            for (unsigned int i = 0; i < mcb2.size(); i++){
//                for (int j = LO; j <= ordDL1_2; j++){
//                    coeffDL1_2.setCoeff(*coeffDL1_2.getCoeff(orders(j))
//                                            + *mcb2[i].getCoeff(orders(j)), orders(j));
//                }
//            }
//        std::cout<<"mcb2 = "<<**coeffDL1_2.getCoeff()<<std::endl;
//        return coeffDL1_2.getCoeff();
//        }
//        case 3:
//        {
//            std::vector<WilsonCoefficient>& mcb3 = model.getMyMatching() -> CMDLijjj(3);
//            orders ordDL1_3 = coeffDL1_3.getOrder();
//            for (unsigned int i = 0; i < mcb3.size(); i++){
//                for (int j = LO; j <= ordDL1_3; j++){
//                    coeffDL1_3.setCoeff(*coeffDL1_3.getCoeff(orders(j))
//                                            + *mcb3[i].getCoeff(orders(j)), orders(j));
//                }
//            }
//        std::cout<<"mcb3 = "<<**coeffDL1_3.getCoeff()<<std::endl;
//        return coeffDL1_3.getCoeff();
//        }
//        default:
//        {
//            throw std::runtime_error("wrong input for the lepton flag in li->lj gamma decays");
//        }
//    }
}
