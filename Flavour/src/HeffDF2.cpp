/* 
 * File:   HeffDF2.cpp
 * Author: silvest
 * 
 * Created on April 28, 2011, 4:33 PM
 */

#include "HeffDF2.h"
#include <sstream>
#include <QCD.h>

HeffDF2::HeffDF2(const StandardModel & SM, StandardModelMatching & SM_Matching) :
        model(SM), modelmatching(SM_Matching), coeff(5, NDR, NLO),
        u(5, NDR, NLO, SM), drNDRLRI(5, 5, 0) {
    double Nc = SM.getNc();
    drNDRLRI(0, 0) = -(((-1. + Nc)*(-7. + log(4096.))) / Nc);
    drNDRLRI(1,1) = (-2.*(-1. + 6.*Nc*Nc - 8.*log(2.) + Nc*(-13. + log(1024.))))/(3.*Nc);
    drNDRLRI(1,2) = (-2.*(13. - 10.*log(2.) + Nc*(-5. + log(256.))))/(3.*Nc);
    drNDRLRI(2,1) = (-8. + 6.*Nc*Nc + 20.*log(2.) - 8.*Nc*(1. + log(4.)))/(3.*Nc);
    drNDRLRI(2,2) = (2.*(4. + Nc - 10.*Nc*log(2.) + log(256.)))/(3.*Nc);
    drNDRLRI(3,3) = (2. - 4.*Nc*Nc + log(4.))/Nc;
    drNDRLRI(3,4) = 2. - log(4.);
    drNDRLRI(4,3) = -2.*(1. + log(2.));
    drNDRLRI(4,4) = (2. + log(4.))/Nc;
}

HeffDF2::~HeffDF2() {
}

vector<complex>** HeffDF2::ComputeCoeff(double mu, schemes scheme) {

    const std::vector<WilsonCoefficient>& mc = modelmatching.CMdf2(model);

    coeff.setMu(mu); // also assign coeff to zero

    orders ordDF2 = coeff.getOrder();
    for (int i = 0; i < mc.size(); i++)
        for (int j = LO; j <= ordDF2; j++)
            for (int k = LO; k <= j; k++)
                coeff.setCoeff(*coeff.getCoeff(orders(j)) +
                    u.Df2Evol(mu, mc[i].getMu(), orders(k), mc[i].getScheme()) *
                    (*(mc[i].getCoeff(orders(j - k)))), orders(j));

    ChangeScheme(scheme, mc[0].getScheme(), ordDF2);

    coeff.setScheme(scheme);

    return coeff.getCoeff();
}

void HeffDF2::ChangeScheme(schemes schout, schemes schin, orders order) {
    if (schout == schin) return;
    switch(schin) {
        case NDR:
            switch(schout) {
                case LRI:
                    coeff.setCoeff(*coeff.getCoeff(NLO) -
                    model.Als(coeff.getMu())/4./M_PI*drNDRLRI.transpose()*
                    (*coeff.getCoeff(LO)), NLO);
                    break;
                default:
                    throw "HeffDF2::ChangeScheme(): out scheme not implemented";
            }
        default:
            throw "HeffDF2::ChangeScheme(): in scheme not implemented";
    }
}

