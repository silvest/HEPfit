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
: model(SM), evolDF1(new EvolDF1_diujlknu(5, NDR, NLO, SM)) {
    for(int i = 0; i < 27; i++)
        coeffdiujleptonknu.push_back(WilsonCoefficient(5,NDR,LO));
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
    int index = i+3*j+9*k;
    coeffdiujleptonknu[index].resetCoefficient();
    coeffdiujleptonknu[index].setMu(mu);
    coeffdiujleptonknu[index].setScheme(mcu[0].getScheme());

    //Hard-coding the Sirlin factor for the moment, Mauro's student can implement the full calculation. 
    //Since it includes the QED running, it should be applied to the SM only. The NP instead will be run using the QED and QCD anomalous dimensions
    //The SM V_LL coefficient returned here therefore includes the leading-log short-distance electroweak correction
    //eta_EW = 1 + alpha/pi ln(M_Z/mu) (A. Sirlin, Nucl. Phys. B196 (1982) 83), the same factor that the LO QED running
    //in EvolDF1_diujlknu (gamma = -4 alpha/(4 pi)) gives the NP V_LL coefficient, so SM and NP are normalised alike.
    //Observables whose inputs already contain S_EW = eta_EW^2 (Pmunu: delPl2em) must divide it out.

    orders ordDF1 = coeffdiujleptonknu[index].getOrder();
    for (unsigned int i = 0; i < mcu.size(); i++) {
        if(i==0) {
            gslpp::vector<gslpp::complex> coeffSM(*mcu[0].getCoeff(LO));
            coeffSM.assign(0, coeffSM(0) * (1. + model.getAle() / M_PI * log(model.getMz() / mu)));
            coeffdiujleptonknu[index].setCoeff(coeffSM, orders(LO));
        }
        else
        for (int j = LO; j <= ordDF1; j++) {
            for (int k = LO; k <= j; k++) {
                coeffdiujleptonknu[index].setCoeff(*coeffdiujleptonknu[index].getCoeff(orders(j))
                        + evolDF1->Df1diujlknuEvol(mu, mcu[i].getMu(), orders(k), mcu[i].getScheme()) *
                        (*mcu[i].getCoeff(orders(j - k))), orders(j));
            }
        }
    }
    return coeffdiujleptonknu[index].getCoeff();
}