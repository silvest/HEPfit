/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cmath>
#include <iostream>
#include "EWSM.h"

// include the imaginary part of O(alpha alpha_s) contribution (for a test)
//#define WITHIMTWOLOOPQCD

const double EWSM::Mw_error = 0.00001; /* 0.01 MeV */


EWSM::EWSM(const StandardModel& SM_i) 
: SM(SM_i) 
{
    flag_order[EW1] = true;
    flag_order[EW1QCD1] = true;
    flag_order[EW1QCD2] = true;
    flag_order[EW2] = true;
    flag_order[EW2QCD1] = true;
    flag_order[EW3] = true;

    bUseCacheEWSM = true;// use caches in the current class
    //bUseCacheEWSM = false;// do not use caches in the current class (for test)
    
    //std::string Model = SM.ModelName();
    //std::cout << "Model in EWSM: " << Model << std::endl;
 
    myCache = new EWSMcache(SM);
    myOneLoopEW = new EWSMOneLoopEW(*myCache);
    myTwoLoopQCD = new EWSMTwoLoopQCD(*myCache);
    myThreeLoopQCD = new EWSMThreeLoopQCD(*myCache);
    myTwoLoopEW = new EWSMTwoLoopEW(*myCache);
    myThreeLoopEW2QCD = new EWSMThreeLoopEW2QCD(*myCache);
    myThreeLoopEW = new EWSMThreeLoopEW(*myCache);
    myApproximateFormulae = new EWSMApproximateFormulae(SM);   

    myTwoFermionsLEP2 = new EWSMTwoFermionsLEP2(SM, *myCache);

    /* Default flags (see also StandardModel::setEWSMflags(), etc.) */
    schemeMw = APPROXIMATEFORMULA;
    schemeRhoZ = NORESUM;
    schemeKappaZ = APPROXIMATEFORMULA;
    
    // Initializations of the caches
    DeltaAlphaLepton_cache = 0.0;
    DeltaAlpha_cache = 0.0;
    Mw_cache = 0.0;
    GammaW_cache = 0.0;
    int i, j;
    for (j=0; j<NumSMParams; ++j) {
        DeltaAlphaLepton_params_cache[j] = 0.0;
        DeltaAlpha_params_cache[j] = 0.0;
        Mw_params_cache[j] = 0.0;
        GammaW_params_cache[j] = 0.0;
    }
    for (i=0; i<6; ++i) {
        rhoZ_l_cache[i] = complex(0.0, 0.0, false);
        rhoZ_q_cache[i] = complex(0.0, 0.0, false);
        kappaZ_l_cache[i] = complex(0.0, 0.0, false);
        kappaZ_q_cache[i] = complex(0.0, 0.0, false);
        for (j=0; j<NumSMParams; ++j) {
            rhoZ_l_params_cache[i][j] = 0.0;
            rhoZ_q_params_cache[i][j] = 0.0;
            kappaZ_l_params_cache[i][j] = 0.0;
            kappaZ_q_params_cache[i][j] = 0.0;
        }
    }
    schemeMw_cache = schemes_EW_size;
    schemeRhoZ_cache = schemes_EW_size;
    schemeKappaZ_cache = schemes_EW_size;
}


////////////////////////////////////////////////////////////////////////

bool EWSM::checkSMparams(double Params_cache[], const bool bUpdate) const
{
    // 11 parameters in QCD:
    // "AlsMz", "Mz", "mup", "mdown", "mcharm", "mstrange", "mtop", "mbottom",
    // "mut", "mub", "muc"
    // 13 parameters in StandardModel
    // "GF", "ale", "dAle5Mz", "mHl", 
    // "mneutrino_1", "mneutrino_2", "mneutrino_3", "melectron", "mmu", "mtau",
    // "delMw", "delSin2th_l", "delGammaZ"

    // Note: When modifying the array below, the constant NumSMParams has also
    // to be modified accordingly.
    double SMparams[NumSMParams] = { 
        SM.getAlsMz(), SM.getMz(), SM.getGF(), SM.getAle(), SM.getDAle5Mz(),
        SM.getMHl(), SM.getMtpole(), 
        SM.getLeptons(SM.NEUTRINO_1).getMass(), 
        SM.getLeptons(SM.NEUTRINO_2).getMass(),
        SM.getLeptons(SM.NEUTRINO_3).getMass(),
        SM.getLeptons(SM.ELECTRON).getMass(),
        SM.getLeptons(SM.MU).getMass(),
        SM.getLeptons(SM.TAU).getMass(),
        SM.getQuarks(SM.UP).getMass(),
        SM.getQuarks(SM.DOWN).getMass(),
        SM.getQuarks(SM.CHARM).getMass(),
        SM.getQuarks(SM.STRANGE).getMass(),
        SM.getQuarks(SM.BOTTOM).getMass(),
        SM.getMut(), SM.getMub(), SM.getMuc(),
        SM.getDelMw(), SM.getDelSin2th_l(), SM.getDelGammaZ() 
    };
        
    // check updated parameters
    bool bCache = true;
    for(int i=0; i<NumSMParams; ++i) {
        if (Params_cache[i] != SMparams[i]) { 
            if (bUpdate) Params_cache[i] = SMparams[i];
            bCache &= false;
        }
    }
    
    return bCache;
}


bool EWSM::checkScheme(schemes_EW& scheme_cache, const schemes_EW scheme_current,
                       const bool bUpdate) const
{
    bool bCache = true;
    if (scheme_cache != scheme_current) {
        if (bUpdate) scheme_cache = scheme_current;
        bCache = false;
    }

    return bCache;
}


////////////////////////////////////////////////////////////////////////

double EWSM::DeltaAlphaLepton(const double s) const 
{
    if (s==SM.getMz()*SM.getMz())
        if (bUseCacheEWSM)
            if (checkSMparams(DeltaAlphaLepton_params_cache))
                return DeltaAlphaLepton_cache;
    
    double DeltaAlphaL = 0.0;
    if (flag_order[EW1]) 
        DeltaAlphaL += myOneLoopEW->DeltaAlpha_l(s);
    if (flag_order[EW1QCD1]) 
        DeltaAlphaL += myTwoLoopQCD->DeltaAlpha_l(s);
    if (flag_order[EW1QCD2]) 
        DeltaAlphaL += myThreeLoopQCD->DeltaAlpha_l(s);
    if (flag_order[EW2]) 
        DeltaAlphaL += myTwoLoopEW->DeltaAlpha_l(s);
    if (flag_order[EW2QCD1]) 
        DeltaAlphaL += myThreeLoopEW2QCD->DeltaAlpha_l(s);
    if (flag_order[EW3]) 
        DeltaAlphaL += myThreeLoopEW->DeltaAlpha_l(s);

    if (s==SM.getMz()*SM.getMz())
        DeltaAlphaLepton_cache = DeltaAlphaL;
    return DeltaAlphaL;
}


double EWSM::DeltaAlphaL5q() const 
{
    double Mz2 = SM.getMz()*SM.getMz();    
    return (DeltaAlphaLepton(Mz2) + SM.getDAle5Mz());
}


double EWSM::DeltaAlphaTop(const double s) const 
{
    double DeltaAlpha = 0.0;
    if (flag_order[EW1]) 
        DeltaAlpha += myOneLoopEW->DeltaAlpha_t(s);
    if (flag_order[EW1QCD1]) 
        DeltaAlpha += myTwoLoopQCD->DeltaAlpha_t(s);
    if (flag_order[EW1QCD2]) 
        DeltaAlpha += myThreeLoopQCD->DeltaAlpha_t(s);
    if (flag_order[EW2]) 
        DeltaAlpha += myTwoLoopEW->DeltaAlpha_t(s);
    if (flag_order[EW2QCD1]) 
        DeltaAlpha += myThreeLoopEW2QCD->DeltaAlpha_t(s);
    if (flag_order[EW3]) 
        DeltaAlpha += myThreeLoopEW->DeltaAlpha_t(s);

    return DeltaAlpha; 
}


double EWSM::DeltaAlpha() const 
{
    if (bUseCacheEWSM)
        if (checkSMparams(DeltaAlpha_params_cache))
            return DeltaAlpha_cache;
    
    double Mz2 = SM.getMz()*SM.getMz();
    DeltaAlpha_cache = DeltaAlphaL5q() + DeltaAlphaTop(Mz2);
    return DeltaAlpha_cache; 
}


double EWSM::alphaMz() const 
{
    return (SM.getAle()/(1.0 - DeltaAlpha()));
}


////////////////////////////////////////////////////////////////////////

double EWSM::Mw0() const
{
    return ( sqrt(c02())*SM.getMz() );
}

double EWSM::s02() const
{
    double tmp = 1.0 - 4.0*M_PI*alphaMz()/sqrt(2.0)/SM.getGF()/SM.getMz()/SM.getMz();
    if (tmp < 0.0)
        throw std::runtime_error("Error in EWSM::s02()");

    return ( ( 1.0 - sqrt(tmp) )/2.0 );
}

double EWSM::c02() const
{
    return ( 1.0 - s02() );
}


////////////////////////////////////////////////////////////////////////

double EWSM::Mw_SM() const 
{
    /* Debug */
    //std::cout << std::boolalpha
    //          << checkScheme(schemeMw_cache,schemeMw,false)
    //          << " [cache:" << schemeMw_cache
    //          << " current:" << schemeMw << "]" << std::endl;

    if (bUseCacheEWSM)
        if (checkSMparams(Mw_params_cache)
                && checkScheme(schemeMw_cache,schemeMw))
            return Mw_cache;

    double Mw;
    if (schemeMw==APPROXIMATEFORMULA)
        Mw = myApproximateFormulae->Mw();
    else if (schemeMw==FIXED) {
        Mw = 80.385;
        std::cout << "Mw is fixed to " << Mw << " GeV for test!" << std::endl;
    } else {
        //std::cout << std::setprecision(12) 
        //          << "TEST: Mw_tree = " << SM.Mw_tree() << std::endl;
        
        double DeltaRho[orders_EW_size], DeltaR_rem[orders_EW_size];
        ComputeDeltaRho(SM.Mw_tree(), DeltaRho);
        ComputeDeltaR_rem(SM.Mw_tree(), DeltaR_rem);
        Mw = resumMw(SM.Mw_tree(), DeltaRho, DeltaR_rem);
        
        /* Mw from iterations */
        double Mw_org = SM.Mw_tree();
        while (fabs(Mw - Mw_org) > Mw_error) {
            Mw_org = Mw;
            ComputeDeltaRho(Mw, DeltaRho);
            ComputeDeltaR_rem(Mw, DeltaR_rem);
            Mw = resumMw(Mw, DeltaRho, DeltaR_rem);
            /* TEST */
            //int prec_def = std::cout.precision();
            //std::cout << std::setprecision(12) << "TEST: Mw_org = " << Mw_org 
            //        << "  Mw_new = " << Mw << std::endl;
            //std::cout.precision(prec_def);
        }
    } 

    Mw_cache = Mw;
    return Mw;
}


double EWSM::DeltaR_SM() const 
{
    /* in the experimental/running-width scheme */
    double Mw = Mw_SM();
    double sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
    double tmp = sqrt(2.0)*SM.getGF()*sW2*Mw*Mw/M_PI/SM.getAle();
    if (schemeMw==NORESUM || schemeMw==APPROXIMATEFORMULA) {
        return (tmp - 1.0);
    } else {
        return (1.0 - 1.0/tmp);
    }
}


double EWSM::cW2_SM() const 
{
    double Mw = Mw_SM();
    return ( Mw*Mw/SM.getMz()/SM.getMz() );
}


double EWSM::sW2_SM() const 
{
    return ( 1.0 - cW2_SM() );
}


////////////////////////////////////////////////////////////////////////

double EWSM::Mzbar() const
{
    double Gz = 2.4952; // experimental data
    return ( SM.getMz() - Gz*Gz/2.0/SM.getMz() );
}


double EWSM::MwbarFromMw(const double Mw) const
{
    double AlsMw = SM.Als(Mw, FULLNNLO);
    double Gw_SM = 3.0*SM.getGF()*pow(Mw, 3.0)/2.0/sqrt(2.0)/M_PI
                   *(1.0 + 2.0*AlsMw/3.0/M_PI);

    return ( Mw - Gw_SM*Gw_SM/2.0/Mw );
}


double EWSM::MwFromMwbar(const double Mwbar) const
{
    double AlsMw = SM.Als(Mwbar, FULLNNLO);
    double Gw_SM = 3.0*SM.getGF()*pow(Mwbar, 3.0)/2.0/sqrt(2.0)/M_PI
                   *(1.0 + 2.0*AlsMw/3.0/M_PI);

    return (Mwbar + Gw_SM*Gw_SM/2.0/Mwbar);
}


double EWSM::DeltaRbar_SM() const
{
    double Mwbar_SM = MwbarFromMw(Mw_SM());
    double sW2bar = 1.0 - Mwbar_SM*Mwbar_SM/Mzbar()/Mzbar();
    double tmp = sqrt(2.0)*SM.getGF()*sW2bar*Mwbar_SM*Mwbar_SM/M_PI/SM.getAle();

    return (tmp - 1.0);
}


////////////////////////////////////////////////////////////////////////

complex EWSM::rhoZ_l_SM(const StandardModel::lepton l) const 
{    
    if (schemeRhoZ==APPROXIMATEFORMULA)
        throw std::runtime_error("No approximate formula is available for rhoZ^f"); 
    else {
        
        if (bUseCacheEWSM)        
            if (checkSMparams(rhoZ_l_params_cache[(int)l])
                     && checkScheme(schemeRhoZ_cache,schemeRhoZ))
                return rhoZ_l_cache[(int)l];
        
        double Mw = Mw_SM();
        
        /* compute Delta rho */
        double DeltaRho[orders_EW_size];
        ComputeDeltaRho(Mw, DeltaRho);
        
        /* compute delta rho_rem^f */
        complex deltaRho_rem_f[orders_EW_size];
        deltaRho_rem_f[EW1] = complex(0.0, 0.0, false);        
        deltaRho_rem_f[EW1QCD1] = complex(0.0, 0.0, false);
        deltaRho_rem_f[EW1QCD2] = complex(0.0, 0.0, false);        
        deltaRho_rem_f[EW2] = complex(0.0, 0.0, false);
        deltaRho_rem_f[EW2QCD1] = complex(0.0, 0.0, false);
        deltaRho_rem_f[EW3] = complex(0.0, 0.0, false);    
        if (flag_order[EW1]) 
            deltaRho_rem_f[EW1] = myOneLoopEW->deltaRho_rem_l(l,Mw);
        if (flag_order[EW1QCD1]) 
        #ifdef WITHIMTWOLOOPQCD
            deltaRho_rem_f[EW1QCD1] = complex(myTwoLoopQCD->deltaRho_rem_l(l,Mw).real(), 
                                              myTwoLoopQCD->deltaRho_rem_l(l,Mw).imag(), false);
        #else
            deltaRho_rem_f[EW1QCD1] = complex(myTwoLoopQCD->deltaRho_rem_l(l,Mw).real(), 0.0, false);
        #endif
        if (flag_order[EW1QCD2]) 
            deltaRho_rem_f[EW1QCD2] = complex(myThreeLoopQCD->deltaRho_rem_l(l,Mw).real(), 0.0, false);
        if (flag_order[EW2])
            deltaRho_rem_f[EW2] = complex(myTwoLoopEW->deltaRho_rem_l(l,Mw).real(), 0.0, false);
        if (flag_order[EW2QCD1]) 
            deltaRho_rem_f[EW2QCD1] = complex(myThreeLoopEW2QCD->deltaRho_rem_l(l,Mw).real(), 0.0, false);
        if (flag_order[EW3]) 
            deltaRho_rem_f[EW3] = complex(myThreeLoopEW->deltaRho_rem_l(l,Mw).real(), 0.0, false);  
            
        /* compute Delta rbar_rem */
        double DeltaRbar_rem = 0.0;
        if (flag_order[EW1])
            DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(Mw);    

        /* Re[rho_Z^f] with or without resummation */
        double deltaRho_rem_f_real[orders_EW_size];
        for (int j=0; j<orders_EW_size; ++j)
            deltaRho_rem_f_real[j] = deltaRho_rem_f[j].real();
        double ReRhoZf = resumRhoZ(DeltaRho, deltaRho_rem_f_real, 
                                   DeltaRbar_rem, false);
        
        /* Im[rho_Z^f] without resummation */
        double ImRhoZf = 0.0;
        for (int j=0; j<orders_EW_size; ++j)
            ImRhoZf += deltaRho_rem_f[j].imag();    

        rhoZ_l_cache[(int)l] = complex(ReRhoZf, ImRhoZf, false);
        return (complex(ReRhoZf, ImRhoZf, false));    
    }
}


complex EWSM::rhoZ_q_SM(const StandardModel::quark q) const 
{    
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));
    if (schemeRhoZ==APPROXIMATEFORMULA)
        throw std::runtime_error("No approximate formula is available for rhoZ^f"); 
    else {

        if (bUseCacheEWSM)        
            if (checkSMparams(rhoZ_q_params_cache[(int)q])
                    && checkScheme(schemeRhoZ_cache,schemeRhoZ))
                return rhoZ_q_cache[(int)q];
        
        double Mw = Mw_SM();
        
        /* compute Delta rho */
        double DeltaRho[orders_EW_size];
        ComputeDeltaRho(Mw, DeltaRho);
        
        /* compute delta rho_rem^f */
        complex deltaRho_rem_f[orders_EW_size];
        deltaRho_rem_f[EW1] = complex(0.0, 0.0, false); 
        deltaRho_rem_f[EW1QCD1] = complex(0.0, 0.0, false); 
        deltaRho_rem_f[EW1QCD2] = complex(0.0, 0.0, false); 
        deltaRho_rem_f[EW2] = complex(0.0, 0.0, false); 
        deltaRho_rem_f[EW2QCD1] = complex(0.0, 0.0, false); 
        deltaRho_rem_f[EW3] = complex(0.0, 0.0, false); 
        if (flag_order[EW1]) 
            deltaRho_rem_f[EW1] = myOneLoopEW->deltaRho_rem_q(q,Mw);
        if (flag_order[EW1QCD1]) 
        #ifdef WITHIMTWOLOOPQCD
            deltaRho_rem_f[EW1QCD1] = complex(myTwoLoopQCD->deltaRho_rem_q(q,Mw).real(), 
                                              myTwoLoopQCD->deltaRho_rem_q(q,Mw).imag(), false);
        #else
            deltaRho_rem_f[EW1QCD1] = complex(myTwoLoopQCD->deltaRho_rem_q(q,Mw).real(), 0.0, false);
        #endif
        if (flag_order[EW1QCD2]) 
            deltaRho_rem_f[EW1QCD2] = complex(myThreeLoopQCD->deltaRho_rem_q(q,Mw).real(), 0.0, false);
        if (flag_order[EW2]) 
            deltaRho_rem_f[EW2] = complex(myTwoLoopEW->deltaRho_rem_q(q,Mw).real(), 0.0, false);
        if (flag_order[EW2QCD1]) 
            deltaRho_rem_f[EW2QCD1] = complex(myThreeLoopEW2QCD->deltaRho_rem_q(q,Mw).real(), 0.0, false);
        if (flag_order[EW3]) 
            deltaRho_rem_f[EW3] = complex(myThreeLoopEW->deltaRho_rem_q(q,Mw).real(), 0.0, false);  
        
        /* compute Delta rbar_rem */
        double DeltaRbar_rem = 0.0;
        if (flag_order[EW1])
            DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(Mw);    
        
        /* Re[rho_Z^f] with or without resummation */
        double deltaRho_rem_f_real[orders_EW_size];
        for (int j=0; j<orders_EW_size; ++j)
            deltaRho_rem_f_real[j] = deltaRho_rem_f[j].real();
        bool bool_Zbb = false;
        if (q==StandardModel::BOTTOM) bool_Zbb = true;
        double ReRhoZf = resumRhoZ(DeltaRho, deltaRho_rem_f_real, 
                                   DeltaRbar_rem, bool_Zbb);
        
        /* Im[rho_Z^f] without resummation */
        double ImRhoZf = 0.0;
        for (int j=0; j<orders_EW_size; ++j)
            ImRhoZf += deltaRho_rem_f[j].imag();    

        rhoZ_q_cache[(int)q] = complex(ReRhoZf, ImRhoZf, false);
        return (complex(ReRhoZf, ImRhoZf, false));    
    }
}


complex EWSM::kappaZ_l_SM(const StandardModel::lepton l) const 
{

    if (bUseCacheEWSM)    
        if (checkSMparams(kappaZ_l_params_cache[(int)l])
                && checkScheme(schemeKappaZ_cache,schemeKappaZ))
            return kappaZ_l_cache[(int)l];

    double Mw = Mw_SM();
    
    double ReKappaZf = 0.0, ImKappaZf = 0.0;
    if (schemeKappaZ==APPROXIMATEFORMULA) {
        ReKappaZf = myApproximateFormulae->sin2thetaEff_l(l)/sW2_SM(); 
        ImKappaZf = myOneLoopEW->deltaKappa_rem_l(l,Mw).imag();
        #ifdef WITHIMTWOLOOPQCD
        ImKappaZf += myTwoLoopQCD->deltaKappa_rem_l(l,Mw).imag();

        /* TEST */
        //ImKappaZf -= SM.getAle()*SM.getAlsMz()/24.0/M_PI*(cW2_SM() - sW2_SM())/sW2_SM()/sW2_SM();        
        #endif         
    } else {
        /* compute Delta rho */
        double DeltaRho[orders_EW_size];
        ComputeDeltaRho(Mw, DeltaRho);
        
        /* compute delta kappa_rem^f */
        complex deltaKappa_rem_f[orders_EW_size];
        deltaKappa_rem_f[EW1] = complex(0.0, 0.0, false); 
        deltaKappa_rem_f[EW1QCD1] = complex(0.0, 0.0, false); 
        deltaKappa_rem_f[EW1QCD2] = complex(0.0, 0.0, false); 
        deltaKappa_rem_f[EW2] = complex(0.0, 0.0, false); 
        deltaKappa_rem_f[EW2QCD1] = complex(0.0, 0.0, false);         
        deltaKappa_rem_f[EW3] = complex(0.0, 0.0, false); 
        if (flag_order[EW1]) 
            deltaKappa_rem_f[EW1] = myOneLoopEW->deltaKappa_rem_l(l,Mw);
        if (flag_order[EW1QCD1]) 
        #ifdef WITHIMTWOLOOPQCD
            deltaKappa_rem_f[EW1QCD1] = complex(myTwoLoopQCD->deltaKappa_rem_l(l,Mw).real(), 
                                                myTwoLoopQCD->deltaKappa_rem_l(l,Mw).imag(), false);
        #else
            deltaKappa_rem_f[EW1QCD1] = complex(myTwoLoopQCD->deltaKappa_rem_l(l,Mw).real(), 0.0, false);
        #endif
        if (flag_order[EW1QCD2]) 
            deltaKappa_rem_f[EW1QCD2] = complex(myThreeLoopQCD->deltaKappa_rem_l(l,Mw).real(), 0.0, false);
        if (flag_order[EW2]) 
            deltaKappa_rem_f[EW2] = complex(myTwoLoopEW->deltaKappa_rem_l(l,Mw).real(), 0.0, false);
        if (flag_order[EW2QCD1]) 
            deltaKappa_rem_f[EW2QCD1] = complex(myThreeLoopEW2QCD->deltaKappa_rem_l(l,Mw).real(), 0.0, false);
        if (flag_order[EW3]) 
            deltaKappa_rem_f[EW3] = complex(myThreeLoopEW->deltaKappa_rem_l(l,Mw).real(), 0.0, false);   
    
        /* compute Delta rbar_rem */
        double DeltaRbar_rem = 0.0;
        if (flag_order[EW1])
            DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(Mw);    

        /* Re[kappa_Z^f] with or without resummation */
        double deltaKappa_rem_f_real[orders_EW_size];
        for (int j=0; j<orders_EW_size; ++j)
            deltaKappa_rem_f_real[j] = deltaKappa_rem_f[j].real();
        ReKappaZf = resumKappaZ(DeltaRho, deltaKappa_rem_f_real, 
                                DeltaRbar_rem, false);
        
        /* O(alpha^2) correction to Re[kappa_Z^f] from the Z-gamma mixing */
        ReKappaZf += 35.0*alphaMz()*alphaMz()/18.0/sW2_SM()
                     *(1.0 - 8.0/3.0*ReKappaZf*sW2_SM());
    
        /* Im[kappa_Z^f] without resummation */
        for (int j=0; j<orders_EW_size; ++j)
            ImKappaZf += deltaKappa_rem_f[j].imag();   
    }
        
    kappaZ_l_cache[(int)l] = complex(ReKappaZf, ImKappaZf, false);
    return (complex(ReKappaZf, ImKappaZf, false));       
}


complex EWSM::kappaZ_q_SM(const StandardModel::quark q) const 
{
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));
    
    if (bUseCacheEWSM)     
        if (checkSMparams(kappaZ_q_params_cache[(int)q])
                && checkScheme(schemeKappaZ_cache,schemeKappaZ))
            return kappaZ_q_cache[(int)q];

    double Mw = Mw_SM();
    
    double ReKappaZf = 0.0,  ImKappaZf = 0.0;    
    if (schemeKappaZ==APPROXIMATEFORMULA) {
        ReKappaZf = myApproximateFormulae->sin2thetaEff_q(q)/sW2_SM(); 
        ImKappaZf = myOneLoopEW->deltaKappa_rem_q(q,Mw).imag();
        #ifdef WITHIMTWOLOOPQCD
        ImKappaZf += myTwoLoopQCD->deltaKappa_rem_q(q,Mw).imag();
        
        /* TEST */
        //ImKappaZf -= SM.getAle()*SM.getAlsMz()/24.0/M_PI*(cW2_SM() - sW2_SM())/sW2_SM()/sW2_SM();        
        #endif                   
    } else { 
        /* compute Delta rho */
        double DeltaRho[orders_EW_size];
        ComputeDeltaRho(Mw, DeltaRho); 

        /* compute delta kappa_rem^f */        
        complex deltaKappa_rem_f[orders_EW_size];
        deltaKappa_rem_f[EW1] = complex(0.0, 0.0, false);
        deltaKappa_rem_f[EW1QCD1] = complex(0.0, 0.0, false);
        deltaKappa_rem_f[EW1QCD2] = complex(0.0, 0.0, false);
        deltaKappa_rem_f[EW2] = complex(0.0, 0.0, false);
        deltaKappa_rem_f[EW2QCD1] = complex(0.0, 0.0, false);
        deltaKappa_rem_f[EW3] = complex(0.0, 0.0, false);
        if (flag_order[EW1]) 
            deltaKappa_rem_f[EW1] = myOneLoopEW->deltaKappa_rem_q(q,Mw);
        if (flag_order[EW1QCD1]) 
        #ifdef WITHIMTWOLOOPQCD
            deltaKappa_rem_f[EW1QCD1] = complex(myTwoLoopQCD->deltaKappa_rem_q(q,Mw).real(), 
                                                myTwoLoopQCD->deltaKappa_rem_q(q,Mw).imag(), false);
        #else
            deltaKappa_rem_f[EW1QCD1] = complex(myTwoLoopQCD->deltaKappa_rem_q(q,Mw).real(), 0.0, false);
        #endif
        if (flag_order[EW1QCD2]) 
            deltaKappa_rem_f[EW1QCD2] = complex(myThreeLoopQCD->deltaKappa_rem_q(q,Mw).real(), 0.0, false);
        if (flag_order[EW2]) 
            deltaKappa_rem_f[EW2] = complex(myTwoLoopEW->deltaKappa_rem_q(q,Mw).real(), 0.0, false);
        if (flag_order[EW2QCD1]) 
            deltaKappa_rem_f[EW2QCD1] = complex(myThreeLoopEW2QCD->deltaKappa_rem_q(q,Mw).real(), 0.0, false);
        if (flag_order[EW3]) 
            deltaKappa_rem_f[EW3] = complex(myThreeLoopEW->deltaKappa_rem_q(q,Mw).real(), 0.0, false);   
    
        /* compute Delta rbar_rem */
        double DeltaRbar_rem = 0.0;
        if (flag_order[EW1])
            DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(Mw);    
        
        /* Re[kappa_Z^f] with or without resummation */
        double deltaKappa_rem_f_real[orders_EW_size];
        for (int j=0; j<orders_EW_size; ++j)
            deltaKappa_rem_f_real[j] = deltaKappa_rem_f[j].real();
        bool bool_Zbb = false;
        if (q==StandardModel::BOTTOM) bool_Zbb = true;
        ReKappaZf = resumKappaZ(DeltaRho, deltaKappa_rem_f_real, 
                                DeltaRbar_rem, bool_Zbb);
        
        /* O(alpha^2) correction to Re[kappa_Z^f] from the Z-gamma mixing */
        ReKappaZf += 35.0*alphaMz()*alphaMz()/18.0/sW2_SM()
                     *(1.0 - 8.0/3.0*ReKappaZf*sW2_SM());

        /* Im[kappa_Z^f] without resummation */
        for (int j=0; j<orders_EW_size; ++j)
            ImKappaZf += deltaKappa_rem_f[j].imag();    
    }
    
    kappaZ_q_cache[(int)q] = complex(ReKappaZf, ImKappaZf, false);
    return (complex(ReKappaZf, ImKappaZf, false));       
}


complex EWSM::gVl_SM(const StandardModel::lepton l) const
{
    double Ql = SM.getLeptons(l).getCharge();
    return ( gAl_SM(l)
             *(1.0 - 4.0*fabs(Ql)*(kappaZ_l_SM(l))*sW2_SM()) );
}


complex EWSM::gVq_SM(const StandardModel::quark q) const
{
    double Qq = SM.getQuarks(q).getCharge();
    return ( gAq_SM(q)
             *(1.0 - 4.0*fabs(Qq)*(kappaZ_q_SM(q))*sW2_SM()) );
}


complex EWSM::gAl_SM(const StandardModel::lepton l) const
{
    double I3l = SM.getLeptons(l).getIsospin();
    return ( sqrt(rhoZ_l_SM(l))*I3l );
}


complex EWSM::gAq_SM(const StandardModel::quark q) const
{
    double I3q = SM.getQuarks(q).getIsospin();
    return ( sqrt(rhoZ_q_SM(q))*I3q );
}


////////////////////////////////////////////////////////////////////////

complex EWSM::rhoZ_l(const StandardModel::lepton l) const
{
    return rhoZ_l_SM(l);
}


complex EWSM::rhoZ_q(const StandardModel::quark q) const
{
    return rhoZ_q_SM(q);
}


complex EWSM::kappaZ_l(const StandardModel::lepton l) const
{
    return kappaZ_l_SM(l);
}


complex EWSM::kappaZ_q(const StandardModel::quark q) const
{
    return kappaZ_q_SM(q);
}


complex EWSM::gVl(const StandardModel::lepton l) const
{
    return gVl_SM(l);
}


complex EWSM::gVq(const StandardModel::quark q) const
{
    return gVq_SM(q);
}


complex EWSM::gAl(const StandardModel::lepton l) const
{
    return gAl_SM(l);
}


complex EWSM::gAq(const StandardModel::quark q) const
{
    return gAq_SM(q);
}


////////////////////////////////////////////////////////////////////////  

double EWSM::taub() const 
{
    double taub_tmp = 0.0;
    double Xt = myCache->Xt_GF();
    if (flag_order[EW1]) 
        taub_tmp += -2.0*Xt; 
    if (flag_order[EW1QCD1]) 
        taub_tmp += 2.0/3.0*M_PI*Xt*myCache->alsMt(); 
    if (flag_order[EW1QCD2]) 
        taub_tmp += 0.0;
    if (flag_order[EW2]) 
        taub_tmp += -2.0*Xt*Xt*myTwoLoopEW->tau_2();
    if (flag_order[EW2QCD1]) 
        taub_tmp += 0.0;
    if (flag_order[EW3]) 
        taub_tmp += 0.0;
    
    return taub_tmp;
}


////////////////////////////////////////////////////////////////////////  

complex EWSM::rhoZ_l_SM_FlavorDep(const StandardModel::lepton l) const 
{
    double Mz = SM.getMz(); 
    double Mw = Mw_SM();
    double cW2 = Mw*Mw/Mz/Mz, sW2 = 1.0 - cW2;
    StandardModel::lepton ELE = SM.ELECTRON;
    complex ul = ( 3.0*myCache->vl(ELE,Mw)*myCache->vl(ELE,Mw) 
                   + myCache->al(ELE)*myCache->al(ELE) )/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                 + myOneLoopEW->FW_l(Mz*Mz,ELE,Mw);
    complex uf = ( 3.0*myCache->vl(l,Mw)*myCache->vl(l,Mw) 
                   + myCache->al(l)*myCache->al(l) )/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                 + myOneLoopEW->FW_l(Mz*Mz,l,Mw);
    
    complex dRho = 2.0*(uf - ul);
    dRho *= SM.getAle()/4.0/M_PI/sW2;
    return dRho; 
}


complex EWSM::rhoZ_q_SM_FlavorDep(StandardModel::quark q) const 
{
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));

    /* In the case of BOTTOM, the top contribution has to be subtracted.
     * The remaining contribution is the same as that for DOWN and STRANGE. */
    if (q==StandardModel::BOTTOM) q=StandardModel::DOWN;

    double Mz = SM.getMz(); 
    double Mw = Mw_SM();
    double cW2 = Mw*Mw/Mz/Mz, sW2 = 1.0 - cW2;
    StandardModel::lepton ELE = SM.ELECTRON;
    complex ul = ( 3.0*myCache->vl(ELE,Mw)*myCache->vl(ELE,Mw) 
                   + myCache->al(ELE)*myCache->al(ELE) )/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                 + myOneLoopEW->FW_l(Mz*Mz,ELE,Mw);
    complex uf = ( 3.0*myCache->vq(q,Mw)*myCache->vq(q,Mw) 
                   + myCache->aq(q)*myCache->aq(q) )/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                 + myOneLoopEW->FW_q(Mz*Mz,q,Mw);
    
    complex dRho = 2.0*(uf - ul);
    dRho *= SM.getAle()/4.0/M_PI/sW2;
    return dRho; 
}


complex EWSM::kappaZ_l_SM_FlavorDep(const StandardModel::lepton l) const 
{
    double Mz = SM.getMz(); 
    double Mw = Mw_SM();
    double cW2 = Mw*Mw/Mz/Mz, sW2 = 1.0 - cW2;
    StandardModel::lepton ELE = SM.ELECTRON;
    complex ul = ( 3.0*myCache->vl(ELE,Mw)*myCache->vl(ELE,Mw) 
                   + myCache->al(ELE)*myCache->al(ELE) )/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                 + myOneLoopEW->FW_l(Mz*Mz,ELE,Mw);
    double deltal = myCache->deltal(ELE, Mw);
    complex uf = ( 3.0*myCache->vl(l,Mw)*myCache->vl(l,Mw) 
                   + myCache->al(l)*myCache->al(l) )/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                 + myOneLoopEW->FW_l(Mz*Mz,l,Mw);
    double deltaf = myCache->deltal(l, Mw);
    
    complex dKappa = (deltaf*deltaf - deltal*deltal)/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                     - uf + ul;
    dKappa *= SM.getAle()/4.0/M_PI/sW2;
    return dKappa;
}


complex EWSM::kappaZ_q_SM_FlavorDep(StandardModel::quark q) const 
{
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));

    /* In the case of BOTTOM, the top contribution has to be subtracted.
     * The remaining contribution is the same as that for DOWN and STRANGE. */
    if (q==StandardModel::BOTTOM) q=StandardModel::DOWN;
    
    double Mz = SM.getMz(); 
    double Mw = Mw_SM();
    double cW2 = Mw*Mw/Mz/Mz, sW2 = 1.0 - cW2;
    StandardModel::lepton ELE = SM.ELECTRON;
    complex ul = ( 3.0*myCache->vl(ELE,Mw)*myCache->vl(ELE,Mw) 
                   + myCache->al(ELE)*myCache->al(ELE) )/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                 + myOneLoopEW->FW_l(Mz*Mz,ELE,Mw);
    double deltal = myCache->deltal(ELE, Mw);
    complex uf = ( 3.0*myCache->vq(q,Mw)*myCache->vq(q,Mw) 
                   + myCache->aq(q)*myCache->aq(q) )/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                 + myOneLoopEW->FW_q(Mz*Mz,q,Mw);
    double deltaf = myCache->deltaq(q, Mw);
    
    complex dKappa = (deltaf*deltaf - deltal*deltal)/4.0/cW2*myOneLoopEW->FZ(Mz*Mz,Mw) 
                     - uf + ul;
    dKappa *= SM.getAle()/4.0/M_PI/sW2;
    return dKappa; 
}
    
    
////////////////////////////////////////////////////////////////////////

double EWSM::epsilon1_SM() const
{
    double rhoZe = rhoZ_l_SM(SM.ELECTRON).real();
    double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );

    return DeltaRhoPrime;
}


double EWSM::epsilon2_SM() const
{
    double s_W2 = sW2_SM(), c_W2 = cW2_SM();
    double rhoZe = rhoZ_l_SM(SM.ELECTRON).real();
    double sin2thetaEff = kappaZ_l_SM(SM.ELECTRON).real()*s_W2;
    double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    double DeltaKappaPrime = sin2thetaEff/s02() - 1.0;
    double DeltaRW = 1.0 - M_PI*alphaMz()/(sqrt(2.0)*SM.getGF()*SM.getMz()*SM.getMz()*s_W2*c_W2);

    return ( c02()*DeltaRhoPrime + s02()*DeltaRW/(c02() - s02())
             - 2.0*s02()*DeltaKappaPrime );
}


double EWSM::epsilon3_SM() const
{
    double rhoZe = rhoZ_l_SM(SM.ELECTRON).real();
    double sin2thetaEff = kappaZ_l_SM(SM.ELECTRON).real()*sW2_SM();
    double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    double DeltaKappaPrime = sin2thetaEff/s02() - 1.0;

    return ( c02()*DeltaRhoPrime + (c02() - s02())*DeltaKappaPrime );
}


double EWSM::epsilonb_SM() const
{
    /* epsilon_b from g_A^b
     * see Eq.(13) of IJMP A7, 1031 (1998) by Altarelli et al. */
    //double rhoZe = rhoZ_l_SM(SM.ELECTRON).real();
    //double rhoZb = rhoZ_q_SM(SM.BOTTOM).real();
    //double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    //double eps1 = DeltaRhoPrime;
    //return ( - 1.0 + sqrt(rhoZb)/(1.0 + eps1/2.0) );

    /* epsilon_b from Re(g_V^b/g_A^b), i.e. Re(kappaZ_b)
     * see Eq.(13) of IJMP A7, 1031 (1998) by Altarelli et al. */
    complex kappaZe = kappaZ_l_SM(SM.ELECTRON);
    complex kappaZb = kappaZ_q_SM(SM.BOTTOM);
    if (SM.IsFlagWithoutNonUniversalVC())
        return ( kappaZe.real()/kappaZb.real() - 1.0 );
    else
        return ( (kappaZe.real() + kappaZ_q_SM_FlavorDep(SM.BOTTOM).real())
                 /kappaZb.real() - 1.0 );

    /* epsilon_b from Gamma_b via Eqs.(11), (12) and (16) of IJMP A7,
     * 1031 (1998) by Altarelli et al.
     * Note: mb has to be mb=4.7, since Eq.(16) were derived with this value.
     */
    //double als_Mz = SM.Als(SM.getMz(), FULLNNLO);
    //double delta_als = (als_Mz - 0.119)/M_PI;
    //double delta_alpha = (alphaMz() - 1.0/128.90)/SM.getAle();
    //double Gamma_b_Born = 0.3798*( 1.0 + delta_als - 0.42*delta_alpha);
    //double a = als_Mz/M_PI;
    //double RQCD = 1.0 + 1.2*a - 1.1*a*a - 13.0*a*a*a;
    //double mb = SM.Mrun(SM.getMz(), SM.getQuarks(SM.BOTTOM).getMass(), FULLNNLO);// This is wrong!
    //double mb = 4.7;
    //std::cout << "mb = " << mb << std::endl;
    //double beta = sqrt(1.0 - 4.0*mb*mb/SM.getMz()/SM.getMz());
    //double Nc = 3.0;
    //double factor = SM.getGF()*SM.getMz()*SM.getMz()*SM.getMz()/6.0/M_PI/sqrt(2.0);
    //double Gamma_b = factor*beta*((3.0 - beta*beta)/2.0*gVq_SM(SM.BOTTOM).abs2()
    //                              + beta*beta*gAq_SM(SM.BOTTOM).abs2())
    //                 *Nc*RQCD*(1.0 + alphaMz()/12.0/M_PI);
    //return ( (Gamma_b/Gamma_b_Born - 1.0 - 1.42*epsilon1_SM()
    //          + 0.54*epsilon3_SM() )/2.29 );
}


////////////////////////////////////////////////////////////////////////     

double EWSM::rho_GammaW_l_SM(const StandardModel::lepton li, 
                             const StandardModel::lepton lj) const 
{
    double Mw = Mw_SM();
    double rhoW = 0.0;
    if (flag_order[EW1])
        rhoW = myOneLoopEW->rho_GammaW_l(li,lj,Mw);    
    return rhoW;
}


double EWSM::rho_GammaW_q_SM(const StandardModel::quark qi, 
                             const StandardModel::quark qj) const 
{
    double Mw = Mw_SM();
    double rhoW = 0.0;
    if (flag_order[EW1])
        rhoW = myOneLoopEW->rho_GammaW_q(qi,qj,Mw);    
    return rhoW;
}


double EWSM::GammaW_l_SM(const StandardModel::lepton li, 
                         const StandardModel::lepton lj) const 
{
    if ( ((int)li+2)%2 || ((int)lj+3)%2 ) 
        throw std::runtime_error("Error in EWSM::GammaW_l_SM()"); 
    
    double G0 = SM.getGF()*pow(Mw_SM(),3.0)/6.0/sqrt(2.0)/M_PI;    
    complex V(0.0, 0.0, false);
    if ( (li==StandardModel::NEUTRINO_1 && lj==StandardModel::ELECTRON) ||
         (li==StandardModel::NEUTRINO_2 && lj==StandardModel::MU) ||
         (li==StandardModel::NEUTRINO_3 && lj==StandardModel::TAU) )        
        V.real() = 1.0;
    return ( V.abs2()*G0*rho_GammaW_l_SM(li,lj) ); 
}
    
    
double EWSM::GammaW_q_SM(const StandardModel::quark qi, 
                         const StandardModel::quark qj) const 
{
    if ( ((int)qi+2)%2 || ((int)qj+3)%2 ) 
        throw std::runtime_error("Error in EWSM::GammaW_q_SM()"); 
    
    double G0 = SM.getGF()*pow(Mw_SM(),3.0)/6.0/sqrt(2.0)/M_PI;    
    complex V(0.0, 0.0, false);

    if ( qi==StandardModel::UP && qj==StandardModel::DOWN )
        //V = SM.getCKM().V_ud();
        V = complex(1.0, 0.0, false);
    else if ( qi==StandardModel::UP && qj==StandardModel::STRANGE )
        //V = SM.getCKM().V_us();
        V = complex(0.0, 0.0, false);
    else if ( qi==StandardModel::UP && qj==StandardModel::BOTTOM )
        //V = SM.getCKM().V_ub();
        V = complex(0.0, 0.0, false);
    else if ( qi==StandardModel::CHARM && qj==StandardModel::DOWN )
        //V = SM.getCKM().V_cd();
        V = complex(0.0, 0.0, false);
    else if ( qi==StandardModel::CHARM && qj==StandardModel::STRANGE )
        //V = SM.getCKM().V_cs();
        V = complex(1.0, 0.0, false);
    else if ( qi==StandardModel::CHARM && qj==StandardModel::BOTTOM )
        //V = SM.getCKM().V_cb();
        V = complex(0.0, 0.0, false);
    else if ( qi==StandardModel::TOP || qj==StandardModel::TOP )
        return (0.0);
    double AlsMw = SM.AlsWithInit(Mw_SM(), SM.getAlsMz(), SM.getMz(), FULLNLO);
    //double AlsMw = SM.Als(Mw_SM(), FULLNNLO);
    return ( 3.0*V.abs2()*G0*rho_GammaW_q_SM(qi,qj)*(1.0 + AlsMw/M_PI) );
}


double EWSM::GammaW_SM() const 
{
    if (bUseCacheEWSM)      
        if (checkSMparams(GammaW_params_cache))
            return GammaW_cache;
    
    double GammaW = GammaW_l_SM(SM.NEUTRINO_1, SM.ELECTRON) 
                    + GammaW_l_SM(SM.NEUTRINO_2, SM.MU) 
                    + GammaW_l_SM(SM.NEUTRINO_3, SM.TAU)             
                    + GammaW_q_SM(SM.UP, SM.DOWN) 
                    + GammaW_q_SM(SM.UP, SM.STRANGE) 
                    + GammaW_q_SM(SM.UP, SM.BOTTOM) 
                    + GammaW_q_SM(SM.CHARM, SM.DOWN)
                    + GammaW_q_SM(SM.CHARM, SM.STRANGE) 
                    + GammaW_q_SM(SM.CHARM, SM.BOTTOM);
    GammaW_cache = GammaW;
    return GammaW;
}


////////////////////////////////////////////////////////////////////////

void EWSM::ComputeDeltaRho(const double Mw_i,
                           double DeltaRho[orders_EW_size]) const
{
    if (flag_order[EW1])
        DeltaRho[EW1] = myOneLoopEW->DeltaRho(Mw_i);
    else
        DeltaRho[EW1] = 0.0;
    if (flag_order[EW1QCD1])
        DeltaRho[EW1QCD1] = myTwoLoopQCD->DeltaRho(Mw_i);
    else
        DeltaRho[EW1QCD1] = 0.0;
    if (flag_order[EW1QCD2])
        DeltaRho[EW1QCD2] = myThreeLoopQCD->DeltaRho(Mw_i);
    else
        DeltaRho[EW1QCD2] = 0.0;
    if (flag_order[EW2])
        DeltaRho[EW2] = myTwoLoopEW->DeltaRho(Mw_i);
    else
        DeltaRho[EW2] = 0.0;
    if (flag_order[EW2QCD1])
        DeltaRho[EW2QCD1] = myThreeLoopEW2QCD->DeltaRho(Mw_i);
    else
        DeltaRho[EW2QCD1] = 0.0;
    if (flag_order[EW3])
        DeltaRho[EW3] = myThreeLoopEW->DeltaRho(Mw_i);
    else
        DeltaRho[EW3] = 0.0;
}


void EWSM::ComputeDeltaR_rem(const double Mw_i,
                             double DeltaR_rem[orders_EW_size]) const
{
    if (flag_order[EW1])
        DeltaR_rem[EW1] = myOneLoopEW->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW1] = 0.0;
    if (flag_order[EW1QCD1])
        DeltaR_rem[EW1QCD1] = myTwoLoopQCD->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW1QCD1] = 0.0;
    if (flag_order[EW1QCD2])
        DeltaR_rem[EW1QCD2] = myThreeLoopQCD->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW1QCD2] = 0.0;
    if (flag_order[EW2])
        DeltaR_rem[EW2] = myTwoLoopEW->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW2] = 0.0;
    if (flag_order[EW2QCD1])
        DeltaR_rem[EW2QCD1] = myThreeLoopEW2QCD->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW2QCD1] = 0.0;
    if (flag_order[EW3])
        DeltaR_rem[EW3] = myThreeLoopEW->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW3] = 0.0;
}


////////////////////////////////////////////////////////////////////////

double EWSM::resumMw(const double Mw_i, const double DeltaRho[orders_EW_size],
                     const double DeltaR_rem[orders_EW_size]) const
{
    if ( (schemeMw==APPROXIMATEFORMULA) 
            || (DeltaR_rem[EW2QCD1]!=0.0) 
            || (DeltaR_rem[EW3]!=0.0) )
        throw std::runtime_error("Error in EWSM::resumMw()"); 

    if (!flag_order[EW2] && schemeMw!=NORESUM)
        throw std::runtime_error("Error in EWSM::resumMw()");       
    
    double cW2_TMP = Mw_i*Mw_i/SM.getMz()/SM.getMz();
    double sW2_TMP = 1.0 - cW2_TMP;
    
    double f_AlphaToGF, DeltaRho_sum = 0.0, DeltaRho_G;
    if (schemeMw==NORESUM) {
        for (int j=0; j<orders_EW_size; ++j) {
            //f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)*sW2_TMP*cW2_TMP/M_PI/SM.getAle();
            //if (j==(int)EW1QCD2)
            //    DeltaRho_sum += f_AlphaToGF*DeltaRho[(orders_EW)j];
            //else

            DeltaRho_sum += DeltaRho[(orders_EW)j];
        }
    } else {
        // conversion: alpha(0) --> G_F
        f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)*sW2_TMP*cW2_TMP/M_PI/SM.getAle();
        DeltaRho_sum = f_AlphaToGF*DeltaRho[EW1]
                       + f_AlphaToGF*DeltaRho[EW1QCD1]
                       + f_AlphaToGF*DeltaRho[EW1QCD2]                
                       + pow(f_AlphaToGF,2.0)*DeltaRho[EW2]
                       + pow(f_AlphaToGF,2.0)*DeltaRho[EW2QCD1]
                       + pow(f_AlphaToGF,3.0)*DeltaRho[EW3];
        DeltaRho_G = f_AlphaToGF*DeltaRho[EW1];
    }
        
    double R;
    double DeltaR_rem_sum = 0.0;
    double DeltaR_EW1 = 0.0, DeltaR_EW2_rem = 0.0;
    switch (schemeMw) {
        case NORESUM:
            for (int j=0; j<orders_EW_size; ++j)
                DeltaR_rem_sum += DeltaR_rem[(orders_EW)j];

            // Full EW one-loop contribution (without the full DeltaAlphaL5q)
            DeltaR_EW1 = - cW2_TMP/sW2_TMP*DeltaRho[EW1] + DeltaR_rem[EW1];

            // Full EW two-loop contribution without reducible corrections
            DeltaR_EW2_rem = myApproximateFormulae->DeltaR_TwoLoopEW_rem(Mw_i);

            // subtract the EW two-loop contributions from DeltaRho_sum and DeltaR_rem_sum
            DeltaRho_sum -= DeltaRho[EW2];
            DeltaR_rem_sum -= DeltaR_rem[EW2];

            // R = 1 + Delta r, including the full EW two-loop contribution
            R = 1.0 + DeltaAlphaL5q() - cW2_TMP/sW2_TMP*DeltaRho_sum
                + DeltaR_rem_sum;
            R += DeltaAlphaL5q()*DeltaAlphaL5q() + 2.0*DeltaAlphaL5q()*DeltaR_EW1
                 + DeltaR_EW2_rem;

            break;
        case OMSI:
            // R = 1/(1 - Delta r)
            R = 1.0/(1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum)
                /(1.0 - DeltaAlphaL5q()
                  - DeltaR_rem[EW1] - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2]);
            break;
        case INTERMEDIATE:
            // R = 1/(1 - Delta r)
            R = 1.0/( (1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum)
                      *(1.0 - DeltaAlphaL5q() - DeltaR_rem[EW1]) 
                        - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2] );
            break;        
        case OMSII:
            // R = 1/(1 - Delta r)
            R = 1.0/( (1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum)*(1.0 - DeltaAlphaL5q())
                      - (1.0 + cW2_TMP/sW2_TMP*DeltaRho_G)*DeltaR_rem[EW1]
                      - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2] );
            break;
        default:
            throw std::runtime_error("Error in EWSM::resumMw()");             
    }   

    if (schemeMw == NORESUM) {
        /* Mzbar and Mwbar are defined in the complex-pole scheme. */

        double tmp = 4.0*M_PI*SM.getAle()/sqrt(2.0)/SM.getGF()/Mzbar()/Mzbar();
        if (tmp*R > 1.0) throw std::runtime_error("EWSM::resumMw(): Negative (1-tmp*R)");
        double Mwbar = Mzbar()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp*R));

        return MwFromMwbar(Mwbar);
    } else {
        double tmp = 4.0*M_PI*SM.getAle()/sqrt(2.0)/SM.getGF()/SM.getMz()/SM.getMz();
        if (tmp*R > 1.0) throw std::runtime_error("EWSM::resumMw(): Negative (1-tmp*R)");
    
        return (SM.getMz()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp*R)));
    }
}


double EWSM::resumRhoZ(const double DeltaRho[orders_EW_size],
                       const double deltaRho_rem[orders_EW_size],
                       const double DeltaRbar_rem, const bool bool_Zbb) const 
{
    if ( (schemeRhoZ==APPROXIMATEFORMULA) 
            || (deltaRho_rem[EW1QCD2]!=0.0) 
            || (deltaRho_rem[EW2QCD1]!=0.0) 
            || (deltaRho_rem[EW3]!=0.0) )
        throw std::runtime_error("Error in EWSM::resumRhoZ()");   

    if (!flag_order[EW2] && schemeRhoZ!=NORESUM)
        throw std::runtime_error("Error in EWSM::resumRhoZ()");       
    
    double Mw_TMP = Mw_SM();
    double cW2_TMP = cW2_SM();
    double sW2_TMP = sW2_SM();
    
    double f_AlphaToGF, DeltaRho_sum = 0.0, DeltaRho_G;
    double DeltaRbar_rem_G, deltaRho_rem_G, deltaRho_rem_G2;
    // conversion: alpha(0) --> G_F
    f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)
                  *sW2_TMP*cW2_TMP/M_PI/SM.getAle();
    DeltaRho_sum = f_AlphaToGF*DeltaRho[EW1]
                   + f_AlphaToGF*DeltaRho[EW1QCD1]
                   + f_AlphaToGF*DeltaRho[EW1QCD2]
                   + pow(f_AlphaToGF,2.0)*DeltaRho[EW2]
                   + pow(f_AlphaToGF,2.0)*DeltaRho[EW2QCD1]
                   + pow(f_AlphaToGF,3.0)*DeltaRho[EW3];
    DeltaRho_G = f_AlphaToGF*DeltaRho[EW1];
    DeltaRbar_rem_G = f_AlphaToGF*DeltaRbar_rem;
    deltaRho_rem_G = f_AlphaToGF*(deltaRho_rem[EW1] 
                                  + deltaRho_rem[EW1QCD1]);
    deltaRho_rem_G2 = pow(f_AlphaToGF,2.0)*deltaRho_rem[EW2];

    /* Real parts */
    double rhoZ;
    if (!bool_Zbb) {
        switch (schemeRhoZ) {
            case OMSI:
                rhoZ = (1.0 + deltaRho_rem_G + deltaRho_rem_G2)
                        /(1.0 - DeltaRho_sum*(1.0 - DeltaRbar_rem_G));
                break;
            case INTERMEDIATE:
                rhoZ = (1.0 + deltaRho_rem_G)
                        /(1.0 - DeltaRho_sum*(1.0 - DeltaRbar_rem_G))
                        + deltaRho_rem_G2;            
                break;
            case NORESUM:
            case OMSII:
                rhoZ = 1.0 + DeltaRho_sum - DeltaRho_G*DeltaRbar_rem_G
                        + DeltaRho_G*DeltaRho_G
                        + deltaRho_rem_G*(1.0 + DeltaRho_G) + deltaRho_rem_G2;
                break;
            default:
                throw std::runtime_error("Error in EWSM::resumRhoZ()"); 
        }
    } else { 
        /* Z to bb */
        double OnePlusTaub = 1.0 + taub(); 
        double OnePlusTaub2 = OnePlusTaub*OnePlusTaub;
        double rhoZbL;
        deltaRho_rem_G += f_AlphaToGF*myCache->ale()/4.0/M_PI/sW2_TMP
                          *pow(myCache->Mt()/Mw_TMP, 2.0);        
        switch (schemeRhoZ) {
            case NORESUM: 
                rhoZ = (1.0 + DeltaRho_sum - DeltaRho_G*DeltaRbar_rem_G
                        + DeltaRho_G*DeltaRho_G
                        + deltaRho_rem_G*(1.0 + DeltaRho_G) + deltaRho_rem_G2)
                        *OnePlusTaub2;
                break;
            case OMSI:
                rhoZbL = OnePlusTaub2/(1.0 - DeltaRho_sum);
                rhoZ = rhoZbL/(1.0 - rhoZbL*deltaRho_rem_G);
                break;
            case INTERMEDIATE:
                rhoZbL = OnePlusTaub2/(1.0 - DeltaRho_sum);
                rhoZ = rhoZbL*(1.0 + rhoZbL*deltaRho_rem_G);
                break;
            case OMSII:
                rhoZbL = OnePlusTaub2/(1.0 - DeltaRho_sum);
                rhoZ = rhoZbL*(1.0 + deltaRho_rem_G);
                break;
            default:
                throw std::runtime_error("Error in EWSM::resumRhoZ()");         
        }
    }
    
    return rhoZ;
}


double EWSM::resumKappaZ(const double DeltaRho[orders_EW_size],
                         const double deltaKappa_rem[orders_EW_size],
                         const double DeltaRbar_rem, const bool bool_Zbb) const 
{
    if ( (schemeKappaZ==APPROXIMATEFORMULA)
            || (deltaKappa_rem[EW2QCD1]!=0.0)
            || (deltaKappa_rem[EW3]!=0.0) )
        throw std::runtime_error("Error in EWSM::resumKappaZ()");      

    if (!flag_order[EW2] && schemeKappaZ!=NORESUM)
        throw std::runtime_error("Error in EWSM::resumKappaZ()");       
    
    double Mw_TMP = Mw_SM();
    double cW2_TMP = cW2_SM();
    double sW2_TMP = sW2_SM();
    
    double f_AlphaToGF, DeltaRho_sum = 0.0, DeltaRho_G;
    double DeltaRbar_rem_G, deltaKappa_rem_G, deltaKappa_rem_G2;
    // conversion: alpha(0) --> G_F
    f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)
                  *sW2_TMP*cW2_TMP/M_PI/SM.getAle();
    DeltaRho_sum = f_AlphaToGF*DeltaRho[EW1]
                   + f_AlphaToGF*DeltaRho[EW1QCD1]
                   + f_AlphaToGF*DeltaRho[EW1QCD2]
                   + pow(f_AlphaToGF,2.0)*DeltaRho[EW2]
                   + pow(f_AlphaToGF,2.0)*DeltaRho[EW2QCD1]
                   + pow(f_AlphaToGF,3.0)*DeltaRho[EW3];
    DeltaRho_G = f_AlphaToGF*DeltaRho[EW1];
    DeltaRbar_rem_G = f_AlphaToGF*DeltaRbar_rem;
    deltaKappa_rem_G = f_AlphaToGF*(deltaKappa_rem[EW1] 
                                    + deltaKappa_rem[EW1QCD1]
                                    + deltaKappa_rem[EW1QCD2]);
    deltaKappa_rem_G2 = pow(f_AlphaToGF,2.0)*deltaKappa_rem[EW2];

    /* Real parts */
    double kappaZ;
    if (!bool_Zbb) {
        switch (schemeKappaZ) {
            case OMSI:
                kappaZ = (1.0 + deltaKappa_rem_G + deltaKappa_rem_G2)
                        *(1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum*(1.0 - DeltaRbar_rem_G));
                break;
            case INTERMEDIATE:
                kappaZ = (1.0 + deltaKappa_rem_G)
                        *(1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum*(1.0 - DeltaRbar_rem_G))
                        + deltaKappa_rem_G2;
                break;        
            case NORESUM: 
            case OMSII:
                kappaZ = 1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum
                         - cW2_TMP/sW2_TMP*DeltaRho_G*DeltaRbar_rem_G
                         + deltaKappa_rem_G*(1.0 + cW2_TMP/sW2_TMP*DeltaRho_G)
                         + deltaKappa_rem_G2;
                break;
            default:
                throw std::runtime_error("Error in EWSM::resumKappaZ()"); 
        }
    } else {
        /* Z to bb */
        double OnePlusTaub = 1.0 + taub(); 
        double kappaZbL;
        deltaKappa_rem_G -= f_AlphaToGF*myCache->ale()/8.0/M_PI/sW2_TMP
                            *pow(myCache->Mt()/Mw_TMP, 2.0);
        switch (schemeKappaZ) {
            case NORESUM: 
                kappaZ = (1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum
                          - cW2_TMP/sW2_TMP*DeltaRho_G*DeltaRbar_rem_G
                          + deltaKappa_rem_G*(1.0 + cW2_TMP/sW2_TMP*DeltaRho_G)
                          + deltaKappa_rem_G2)/OnePlusTaub;
                break;
            case OMSI:
                kappaZbL = (1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum)/OnePlusTaub;
                kappaZ = kappaZbL*(1.0 + deltaKappa_rem_G);
                break;
            case INTERMEDIATE:
            case OMSII:
                kappaZbL = (1.0 + cW2_TMP/sW2_TMP*DeltaRho_sum)/OnePlusTaub;
                kappaZ = kappaZbL + deltaKappa_rem_G;
                break;
            default:
                throw std::runtime_error("Error in EWSM::resumKappaZ()"); 
        }
    }

    return kappaZ;
}


