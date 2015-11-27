/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "HeffDLij.h"

HeffDLij::HeffDLij(const StandardModel & SM_i) :
        model(SM_i),
        coeffDLij_1(2, NDR , LO),coeffDLij_2(2, NDR , LO),coeffDLij_3(2, NDR , LO){
}

HeffDLij::~HeffDLij() {
}

gslpp::vector<gslpp::complex>** HeffDLij::ComputeCoeffDLij(int li_lj) {

    switch (li_lj) {
        case 1:
        {
            std::vector<WilsonCoefficient>& mcb1 = model.getMyMatching() -> CMDLij(1);
            orders ordDLij_1 = coeffDLij_1.getOrder();
            coeffDLij_1.resetCoefficient();
            for (unsigned int i = 0; i < mcb1.size(); i++){
                for (int j = LO; j <= ordDLij_1; j++){
                    coeffDLij_1.setCoeff(*coeffDLij_1.getCoeff(orders(j))
                                            + *mcb1[i].getCoeff(orders(j)), orders(j));
                }
            }
        return coeffDLij_1.getCoeff();
        }
        case 2:
        {
            std::vector<WilsonCoefficient>& mcb2 = model.getMyMatching() -> CMDLij(2);
            orders ordDLij_2 = coeffDLij_2.getOrder();
            coeffDLij_2.resetCoefficient();
            for (unsigned int i = 0; i < mcb2.size(); i++){
                for (int j = LO; j <= ordDLij_2; j++){
                    coeffDLij_2.setCoeff(*coeffDLij_2.getCoeff(orders(j))
                                            + *mcb2[i].getCoeff(orders(j)), orders(j));
                }
            }
        return coeffDLij_2.getCoeff();
        }
        case 3:
        {
            std::vector<WilsonCoefficient>& mcb3 = model.getMyMatching() -> CMDLij(3);
            orders ordDLij_3 = coeffDLij_3.getOrder();
            coeffDLij_3.resetCoefficient();
            for (unsigned int i = 0; i < mcb3.size(); i++){
                for (int j = LO; j <= ordDLij_3; j++){
                    coeffDLij_3.setCoeff(*coeffDLij_3.getCoeff(orders(j))
                                            + *mcb3[i].getCoeff(orders(j)), orders(j));
                }
            }
        return coeffDLij_3.getCoeff();
        }
        default:
        {
            throw std::runtime_error("wrong input for the lepton flag in li->lj gamma decays");
        }
    }
}
