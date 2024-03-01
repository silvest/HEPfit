/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */



#include "HeffDF1_diujlknu.h"
#include "StandardModel.h"
#include "EvolDF1_diujlknu.h"

HeffDF1_diujlknu::HeffDF1_diujlknu(const StandardModel & SM)
: model(SM), coeffdiujleptonknu(5, NDR, LO), evolDF1(new EvolDF1_diujlknu(5, NDR, NLO, SM)) {


}

HeffDF1_diujlknu::~HeffDF1_diujlknu() {
}

/******************************************************************************
 * evolution of Wilson Coefficients for M -> lepton nu  
 * The WC are written in the LEFT basis of arxiv:1709.04486, ordered as CnueduVLLkkji^+, CnueduVLRkkji^+, CnueduSRRkkji^+, CnueduSRLkkji^+, CnueduTRRkkji^+
 * the expressions can be found in arxiv:1706.00410 and arxiv:1605.07114
 *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDF1_diujlknu::ComputeCoeffdiujleptonknu(int i, int j, int k, double mu) {

    const std::vector<WilsonCoefficient>& mcu = model.getMatching().CMdiujleptonknu(i, j, k);
    coeffdiujleptonknu.resetCoefficient();
    coeffdiujleptonknu.setMu(mu);
    coeffdiujleptonknu.setScheme(mcu[0].getScheme());

    //Hard-coding the Sirlin factor for the moment, Mauro's student can implement the full calculation
    double Sew;
    if(i==2)
        Sew = 1.0066;
    else
        Sew = 1.0232;

    orders ordDF1 = coeffdiujleptonknu.getOrder();
    for (unsigned int i = 0; i < mcu.size(); i++) {
        if(i==0)
            coeffdiujleptonknu.setCoeff(*coeffdiujleptonknu.getCoeff(LO)*Sew, orders(LO));
        else
        for (int j = LO; j <= ordDF1; j++) {
            for (int k = LO; k <= j; k++) {
                coeffdiujleptonknu.setCoeff(*coeffdiujleptonknu.getCoeff(orders(j))
                        + evolDF1->Df1diujlknuEvol(mu, mcu[i].getMu(), orders(k), mcu[i].getScheme()) *
                        (*mcu[i].getCoeff(orders(j - k))), orders(j));
            }
        }
    }
    return coeffdiujleptonknu.getCoeff();
}