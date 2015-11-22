/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */
    
#include "Heffmueconv.h"

Heffmueconv::Heffmueconv(const StandardModel & SM_i) :
        model(SM_i),
        coeffmueconv(8, NDR , LO){
}

Heffmueconv::~Heffmueconv() {
}

gslpp::vector<gslpp::complex>** Heffmueconv::ComputeCoeffmueconv() {

    std::vector<WilsonCoefficient>& mcb8 = model.getMyMatching() -> CMmueconv();
    orders ordmueconv = coeffmueconv.getOrder();
    coeffmueconv.resetCoefficient();
    for (unsigned int i = 0; i < mcb8.size(); i++){
        for (int j = LO; j <= ordmueconv; j++){
            coeffmueconv.setCoeff(*coeffmueconv.getCoeff(orders(j))
                                    + *mcb8[i].getCoeff(orders(j)), orders(j));
        }
    }

    return coeffmueconv.getCoeff();

}
