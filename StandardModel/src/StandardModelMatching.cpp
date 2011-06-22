/* 
 * File:   StandardModelMatching.cpp
 * Author: silvest
 * 
 * Created on June 9, 2011, 2:16 PM
 */

#include "StandardModelMatching.h"
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>


StandardModelMatching::StandardModelMatching(const StandardModel& SM_i) : mcdf2((unsigned int) 5, NDR, NLO) {
    SM = SM_i;
}

double StandardModelMatching::S0(double x) const {
    return S0(x, x);
}

double StandardModelMatching::S0(double x, double y) const { // Buras 2000 Appendix
    if (fabs(1. - y / x) < LEPS)
        return ((x * (-4. + 15. * x - 12. * x * x + x*x*x +
            6. * x * x * log(x))) / (4. * pow(-1. + x, 3.)));
    else
        return (x * y * ((1. / 4. + 3. / 2. / (1. - x) - 3. / 4. / pow(1. - x, 2.)) *
            log(x) / (x - y) +
            (1. / 4. + 3. / 2. / (1. - y) - 3. / 4. / pow(1. - y, 2.)) *
            log(y) / (y - x) -
            3. / 4. / (1. - x) / (1. - y)));
}

double StandardModelMatching::S0p( double x) const {
    return (x * (-4. + 18. * x + 3. * x*x + x*x*x) / 4. / pow(x - 1., 3.)
            - 9. * x*x*x / 2. / pow(x - 1., 4.) * log(x));
}

double StandardModelMatching::S11(double x) const {
    return (x * (4. - 39. * x + 168. * x*x + 11. * x*x*x) / 4. / pow(x - 1., 3.)
            + 3. * x*x*x * gsl_sf_dilog(1. - x)*(5. + x) / pow(x - 1., 3.)
            + 3. * x * log(x)*(-4. + 24. * x - 36. * x*x - 7. * x*x*x - x*x*x*x) / 2.
            / pow(x - 1., 4.) + 3. * x*x*x * pow(log(x), 2.)*(13. + 4. * x + x*x) / 2.
            / pow(x - 1., 4.));
}

double StandardModelMatching::S18(double x) const {
    return ((-64. + 68. * x + 17. * x*x - 11. * x*x*x) / 4. / pow(x - 1., 2.)
            + pow(M_PI, 2.)*8. / 3. / x + 2. * gsl_sf_dilog(1. - x)*(8. - 24. * x
            + 20. * x*x - x*x*x + 7. * x*x*x*x - pow(x, 5.)) / x / pow(x - 1., 3.)
            + log(x)*(-32. + 68. * x - 32. * x*x + 28. * x*x*x - 3. * x*x*x*x)
            / 2. / pow(x - 1., 3.) + x*x * pow(log(x), 2.)*(4. - 7. * x + 7. * x*x
            - 2. * x*x*x) / 2. / pow(x - 1., 4.));
}

double StandardModelMatching::S1(double x) const {
    return (SM.getCF() * S11(x) + (SM.getNc()-1.)/2./SM.getNc() * S18(x)); 
}

const std::vector<WilsonCoefficient>& StandardModelMatching::CMdf2(const StandardModel& SM_i) {
   
//    if(SM_i == SM)
//        return(vmc);
    
    SM = SM_i;
    
    double gammam = 8.;
    double Bt;
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.mW(), 2.); // always FULLNLO
    complex co = SM.getGF() / 4. / M_PI * SM.mW() * SM.getlamt_d();
    double Nc = SM.getNc();

    vmc.clear();

    switch (mcdf2.getScheme()) {
        case NDR:
            Bt = 5. * (Nc - 1.) / 2. / Nc + 3. * SM.getCF();
            break;
        case HV:
        case LRI:
        default:
            std::stringstream out;
            out << mcdf2.getScheme();
            throw "StandardModel::CMdf2(): scheme " + out.str() + "not implemented";
    }

    mcdf2.setMu(SM.getMuw());
 
    switch (mcdf2.getOrder()) {
        case NNLO:
        case NLO:
            mcdf2.setCoeff(0, co * co * 4. * (SM.Als(SM.getMuw()) / 4. / M_PI * (S1(xt) +
                    Bt * S0(xt, xt) + 2. * gammam * S0p(xt) * log(SM.getMuw() / SM.mW()))), NLO);
        case LO:
            mcdf2.setCoeff(0, co * co * 4. * S0(xt, xt), LO);
            break;
        default:
            std::stringstream out;
            out << mcdf2.getOrder();
            throw "StandardModelMatching::CMdf2(): order " + out.str() + "not implemented";
    }

    vmc.push_back(mcdf2);
    return(vmc);
}

