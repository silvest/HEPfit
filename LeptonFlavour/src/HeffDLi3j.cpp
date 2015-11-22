/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDLi3j.h"

HeffDLi3j::HeffDLi3j(const StandardModel & SM_i) :
        model(SM_i),
        coeffDLi3j_1(20, NDR , LO),coeffDLi3j_2(20, NDR , LO),coeffDLi3j_3(20, NDR , LO),coeffDLi3j_4(20, NDR , LO)
{
}

HeffDLi3j::~HeffDLi3j() {
}

gslpp::vector<gslpp::complex>** HeffDLi3j::ComputeCoeffDLi3j(int li_lj) {

    switch (li_lj) {
        case 1:
        {
            std::vector<WilsonCoefficient>& mcb4 = model.getMyMatching() -> CMDLi3j(1);
            orders ordDLi3j_1 = coeffDLi3j_1.getOrder();
            coeffDLi3j_1.resetCoefficient();
            for (unsigned int i = 0; i < mcb4.size(); i++){
                for (int j = LO; j <= ordDLi3j_1; j++){
                    coeffDLi3j_1.setCoeff(*coeffDLi3j_1.getCoeff(orders(j))
                                            + *mcb4[i].getCoeff(orders(j)), orders(j));
                }
            }
        return coeffDLi3j_1.getCoeff();
        }
        case 2:
        {
            std::vector<WilsonCoefficient>& mcb5 = model.getMyMatching() -> CMDLi3j(2);
            orders ordDLi3j_2 = coeffDLi3j_2.getOrder();
            coeffDLi3j_2.resetCoefficient();
            for (unsigned int i = 0; i < mcb5.size(); i++){
                for (int j = LO; j <= ordDLi3j_2; j++){
                    coeffDLi3j_2.setCoeff(*coeffDLi3j_2.getCoeff(orders(j))
                                            + *mcb5[i].getCoeff(orders(j)), orders(j));
                }
            }
        return coeffDLi3j_2.getCoeff();
        }
        case 3:
        {
            std::vector<WilsonCoefficient>& mcb6 = model.getMyMatching() -> CMDLi3j(3);
            orders ordDLi3j_3 = coeffDLi3j_3.getOrder();
            coeffDLi3j_3.resetCoefficient();
            for (unsigned int i = 0; i < mcb6.size(); i++){
                for (int j = LO; j <= ordDLi3j_3; j++){
                    coeffDLi3j_3.setCoeff(*coeffDLi3j_3.getCoeff(orders(j))
                                            + *mcb6[i].getCoeff(orders(j)), orders(j));
                }
            }
        return coeffDLi3j_3.getCoeff();
        }
        case 4:
        {
            std::vector<WilsonCoefficient>& mcb7 = model.getMyMatching() -> CMDLi3j(4);
            orders ordDLi3j_4 = coeffDLi3j_4.getOrder();
            coeffDLi3j_4.resetCoefficient();
            for (unsigned int i = 0; i < mcb7.size(); i++){
                for (int j = LO; j <= ordDLi3j_4; j++){
                    coeffDLi3j_4.setCoeff(*coeffDLi3j_4.getCoeff(orders(j))
                                            + *mcb7[i].getCoeff(orders(j)), orders(j));
                }
            }
        return coeffDLi3j_4.getCoeff();
        }
        default:
        {
            throw std::runtime_error("wrong input for the lepton flag in li->3lj decays");
        }
    }
}
