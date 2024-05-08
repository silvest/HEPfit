/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Mll.h"
#include "StandardModel.h"
#include "HeffDB1.h"

Mll::Mll(const StandardModel& SM_i, int obsFlag, QCD::meson meson_i, QCD::lepton lep_i)
: ThObservable(SM_i), dgs(SM_i)
{
    meson = meson_i;
    lep = lep_i;
    if (obsFlag > 0 and obsFlag < 5) obs = obsFlag;
    else throw std::runtime_error("obsFlag in Bsmumu(myFlavour, obsFlag) called from ThFactory::ThFactory() can only be 1 (BR) or 2 (BRbar) or 3 (Amumu) or 4 (Smumu)");
    SM.initializeMeson(meson);
    FixedWCbtos = SM.getFlavour().getFlagFixedWCbtos();
    std::vector<std::string> pars = dgs.getParametersForObservable();
    if (FixedWCbtos) pars.push_back("C10_SM");
    setParametersForObservable(pars);
};

double Mll::computeThValue()
{
    computeObs(FULLNLO, FULLNLO_QED);
    double FBs = SM.getMesons(meson).getDecayconst();
    
//    double coupling = SM.getGF() * SM.getGF() * SM.Mw() * SM.Mw() /M_PI /M_PI ; /* Double GF for including EW corrections*/
    double coupling = SM.getGF() * SM.getAle() / 4. / M_PI; /* Single GF for excluding EW corrections*/
 
//    double PRF = pow(coupling, 2.) / M_PI /8. / SM.getMesons(meson).computeWidth() * pow(FBs, 2.) * pow(mlep, 2.) * mBs * beta; /* Double GF for including EW corrections*/
    double PRF = pow(coupling, 2.) / M_PI / SM.getMesons(meson).computeWidth() * pow(FBs, 2.) * pow(mlep, 2.) * mBs * beta; /* Single GF for excluding EW corrections*/
    timeInt = (1. + Amumu * ys) / (1. - ys * ys); // Note modification in form due to algorithm
     
    if (obs == 1) return( PRF * ampSq);
    if (obs == 2) return( PRF * ampSq * timeInt);
    if (obs == 3) return( Amumu );
    if (obs == 4) return( Smumu );
    
    throw std::runtime_error("Bsmumu::computeThValue(): Observable type not defined. Can be only any of (1,2,3,4)");
    return (EXIT_FAILURE);
}

void Mll::computeObs(orders order, orders_qed order_qed)
{   
    double mu = SM.getMub();  
    
    mlep = SM.getLeptons(lep).getMass();
    mBs = SM.getMesons(meson).getMass();
    mW = SM.Mw();
    mb = SM.getQuarks(QCD::BOTTOM).getMass();
    if (meson == QCD::B_S) {
        ms = SM.getQuarks(QCD::STRANGE).getMass();
        CKM_factor = SM.getCKM().computelamt_s();
        ys = dgs.computeThValue()*SM.getMesons(QCD::B_S).getLifetime()/2.;
    } else if (meson == QCD::B_D) {
        ms = SM.getQuarks(QCD::DOWN).getMass();
        CKM_factor = SM.getCKM().computelamt_d();
        ys = 0.;
    }
    chiral = pow(mBs, 2.) / 2. / mlep * mb / (mb + ms);
    beta = sqrt(1. - pow(2. * mlep / mBs, 2.));
    computeAmpSq(order, order_qed, mu);
    Amumu = (absP * absP * cos(2. * argP - phiNP) -  absS * absS * cos(2. * argS - phiNP)) / (absP * absP + absS * absS);
    Smumu = (absP * absP * sin(2. * argP - phiNP) -  absS * absS * sin(2. * argS - phiNP)) / (absP * absP + absS * absS);
}

double Mll::computeAmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Amumu);
}

double Mll::computeSmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Smumu);
}

void Mll::computeAmpSq(orders order, orders_qed order_qed, double mu)
{  
    if (SM.getFlavour().getHDB1().getCoeffsmumu().getOrder() < order % 3){
        std::stringstream out;
        out << order;
        throw std::runtime_error("Bsmumu::computeAmpSq(): required cofficient of "
                                 "order " + out.str() + " not computed");
    }
    /* Temporary usage of MVll class here below  */
//    gslpp::vector<gslpp::complex> ** allcoeffmumu; /* Double GF for including EW corrections*/
//    if (meson == QCD::B_S) allcoeffmumu = SM.getFlavour().ComputeCoeffsmumu(mu, NDR); /* Double GF for including EW corrections*/
//    if (meson == QCD::B_D) allcoeffmumu = SM.getFlavour().ComputeCoeffdmumu(mu, NDR); /* Double GF for including EW corrections*/
    
    allcoeff = SM.getFlavour().ComputeCoeffBMll(mu, lep); /* Single GF for excluding EW corrections*/
    allcoeffprime = SM.getFlavour().ComputeCoeffprimeBMll(mu, lep); /* Single GF for excluding EW corrections*/

//    double alsmu = SM.Als(mu, FULLNNLO, true)/4./M_PI; /* tilde */ /* Double GF for including EW corrections*/
//    double alemu = SM.Ale(mu, FULLNNLO)/4./M_PI; /* tilde */ /* Double GF for including EW corrections*/

//    double sw = sqrt( (M_PI * SM.getAle() ) / ( sqrt(2.) * SM.getGF() * SM.Mw() * SM.Mw()) ); /* Spurious sw */

    C_10p = (*(allcoeffprime[LO]))(9) + (*(allcoeffprime[NLO]))(9);
    C_S = (*(allcoeff[LO]))(10) + (*(allcoeff[NLO]))(10);
    C_Sp = (*(allcoeffprime[LO]))(10) + (*(allcoeffprime[NLO]))(10);
    C_P = (*(allcoeff[LO]))(11) + (*(allcoeff[NLO]))(11);
    C_Pp = (*(allcoeffprime[LO]))(11) + (*(allcoeffprime[NLO]))(11);
    
    if (FixedWCbtos) {
        allcoeff_noSM = SM.getFlavour().ComputeCoeffBMll(mu, lep, true); /* Single GF for excluding EW corrections*/
        C_10 = SM.getOptionalParameter("C10_SM") + ((*(allcoeff_noSM[LO]))(9) + (*(allcoeff_noSM[NLO]))(9));
    }
    else C_10 = ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9));

    if ((order == FULLNLO) && (order_qed == FULLNLO_QED)) {

        switch (order_qed) {
            case FULLNLO_QED:
            {
                /* Implementation to be corrected and updated with new EVO: At present better to use MVll!*/ /* Double GF for including EW corrections*/
                //            gslpp::complex C10_SM = (*(allcoeffmumu[LO]))(7) /alemu  + (*(allcoeffmumu[NLO]))(7) * alsmu/alemu 
                //                    + (*(allcoeffmumu[NNLO]))(7) * alsmu * alsmu/alemu + (*(allcoeffmumu[LO_QED ]))(7) /alsmu
                //                    + (*(allcoeffmumu[NLO_QED11]))(7) + (*(allcoeffmumu[NLO_QED02]))(7) * alemu /alsmu /alsmu 
                //                    + (*(allcoeffmumu[NLO_QED21]))(7) * alsmu 
                //                    + (*(allcoeffmumu[NLO_QED12]))(7) * alemu /alsmu+ (*(allcoeffmumu[NLO_QED22]))(7) * alemu;

                /* Temporary usage of MVll result */
                //            std::cout << " C10_SM " << C10_SM << std::endl;
                //            gslpp::complex C10_SM_plus_NP = CKM_factor * sw * sw * ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9));
                //            std::cout << " C10_SM_plus_NP " << C10_SM_plus_NP / sw / sw / CKM_factor << std::endl;

                //            gslpp::complex CC_P = C10_SM; /* Double GF for including EW corrections*/
                //            gslpp::complex CC_P = CKM_factor * sw * sw * ( C_10 - C_10p + mBs*mBs / ( 2.*mlep*(mb+ms) ) * (C_P - C_Pp) ); /* Spurious sw */
                gslpp::complex CC_P = CKM_factor * (C_10 - C_10p + mBs * mBs / (2. * mlep * (mb + ms)) * (C_P - C_Pp)); /* Single GF for excluding EW corrections*/

                absP = CC_P.abs(); //contains only SM contributions (P, P', S, S' not added)
                argP = CC_P.arg();

                //            gslpp::complex CC_S = CKM_factor * sw * sw * ( beta * mBs*mBs / ( 2.*mlep*(mb+ms) ) * (C_S - C_Sp) ); /* Spurious sw */
                gslpp::complex CC_S = CKM_factor * (beta * mBs * mBs / (2. * mlep * (mb + ms)) * (C_S - C_Sp)); /* Single GF for excluding EW corrections*/

                absS = CC_S.abs();
                argS = CC_S.arg();

                phiNP = 0.;

                ampSq = absP * absP + absS * absS;

            }
                break;
            default:
                std::stringstream out;
                out << order;
                throw std::runtime_error("Bsmumu::computeAmpSq(): order " + out.str() + " not implemented");
        }
    }
        
}


