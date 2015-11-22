/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "Heffgminus2.h"

Heffgminus2::Heffgminus2(const StandardModel & SM_i) :
        model(SM_i),
        coeffgminus2mu(2, NDR , LO){
}

Heffgminus2::~Heffgminus2() {
}

gslpp::vector<gslpp::complex>** Heffgminus2::ComputeCoeffgm2mu() {

    std::vector<WilsonCoefficient>& mcb9 = model.getMyMatching() -> CMgminus2mu();
    orders ordgminus2mu = coeffgminus2mu.getOrder();
    coeffgminus2mu.resetCoefficient();
    for (unsigned int i = 0; i < mcb9.size(); i++){
        for (int j = LO; j <= ordgminus2mu; j++){
            coeffgminus2mu.setCoeff(*coeffgminus2mu.getCoeff(orders(j))
                                    + *mcb9[i].getCoeff(orders(j)), orders(j));
        }
    }

    return coeffgminus2mu.getCoeff();

}
