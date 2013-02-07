/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

//#include <bt/assign/list_oost/assign/list_of.hpp> // for 'map_list_of()'
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdexcept>
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>
#include <gsl/gsl_sf_zeta.h>
#include "StandardModel.h"
#include "CKM.h"
#include "EWSM.h"
#include "StandardModelMatching.h"


const std::string StandardModel::SMvars[NSMvars] = {"GF", "mneutrino_1", "mneutrino_2",
    "mneutrino_3", "melectron", "mmu", "mtau", "lambda", "A", "rhob", "etab", "ale",
    "dAle5Mz", "mHl", "muw", "phiEpsK","DeltaMK", "KbarEpsK", "Dmk", "SM_M12D",
    "delMw", "delSin2th_l"};

/**
 * FixedSMparams: true if all the SM parameters are fixed to constants in the fit. 
 * Flags for the EW precision observables (see EW.h for detail):
 *   EWCHMN: use EW_CHMN class
 *   EWABC:  use EW_ABC class based on the formulae in Eqs.(7)-(14) of IJMP, A7, 
 *           1031-1058 (1998) by Altarelli et al.
 *   EWABC2: use use the approximate formulae in Eqs.(16)-(20) of IJMP, A7, 
 *           1031-1058 (1998) by Altarelli et al.
 *   EWBURGESS: use the formulae for STU contributions by Burgess et al.
 *   R0bApproximate: use the two-loop approximate formula for R_b by Freitas and Huang
 */
const std::string StandardModel::SMflags[NSMflags] 
    = {"FixedAllSMparams", "EWCHMN", "EWABC", "EWABC2", "EWBURGESS", "R0bApproximate", 
       "withoutNonUniversalVCinEpsilons"};


StandardModel::StandardModel(const bool bDebug_i) : QCD(), VCKM(3, 3, 0.), 
        UPMNS(3, 3, 0.), Yu(3, 3, 0.), Yd(3, 3, 0.), Yn(3, 3, 0.), Ye(3, 3, 0.), 
        bDebug(bDebug_i) {
    FlagFixedAllSMparams = false;
    FlagEWCHMN = false;
    FlagEWABC = false;
    FlagEWABC2 = false;
    FlagEWBURGESS = false;
    FlagR0bApproximate = false;
    FlagWithoutNonUniversalVC = false;
    
    leptons[NEUTRINO_1].setCharge(0.);
    leptons[NEUTRINO_2].setCharge(0.);    
    leptons[NEUTRINO_3].setCharge(0.);    
    leptons[ELECTRON].setCharge(-1.);
    leptons[MU].setCharge(-1.);    
    leptons[TAU].setCharge(-1.); 
    leptons[NEUTRINO_1].setIsospin(1./2.);
    leptons[NEUTRINO_2].setIsospin(1./2.);    
    leptons[NEUTRINO_3].setIsospin(1./2.);    
    leptons[ELECTRON].setIsospin(-1./2.);
    leptons[MU].setIsospin(-1./2.);   
    leptons[TAU].setIsospin(-1./2.);
}

bool StandardModel::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return(CheckParameters(DPars));
}

bool StandardModel::PreUpdate(){
    computeCKM = false;
    computeYe = false;
    computeYn = false;
    
    if(!QCD::PreUpdate())  return (false);
    
    return (true);
}

bool StandardModel::Update(const std::map<std::string, double>& DPars) {
    if(!PreUpdate()) return (false);

    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate())  return (false);
    
    return (true);
}

bool StandardModel::PostUpdate(){
    if(!QCD::PostUpdate()) return (false);
    
    if (computeCKM) {
        myCKM.setWolfenstein(lambda, A, rhob, etab);
        myCKM.getCKM(VCKM);
    }
    if (computeYu || computeCKM) {
        Yu = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yu.assign(i, i, this->quarks[UP + 2 * i].getMass() / v() * sqrt(2.));
        Yu = VCKM.transpose()*Yu;
    }
    if (computeYd) {
        for (int i = 0; i < 3; i++)
            Yd.assign(i, i, this->QCD::quarks[DOWN + 2 * i].getMass() / v() * sqrt(2.));
    }
    if (computeYe) {
        for (int i = 0; i < 3; i++)
            Ye.assign(i, i, this->leptons[ELECTRON + 2 * i].getMass() / v() * sqrt(2.));
    }
    if (computeYn) {
        Yn = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yn.assign(i, i, this->leptons[NEUTRINO_1 + 2 * i].getMass() / v() * sqrt(2.));
        Yn = Yn * UPMNS.hconjugate();
    }
    
     return (true);
}

void StandardModel::SetParameter(const std::string name, const double& value) {
    if (name.compare("GF") == 0)
        GF = value;
    else if (name.compare("ale") == 0)
        ale = value;
    else if (name.compare("dAle5Mz") == 0)
        dAle5Mz = value;
    else if (name.compare("mHl") == 0)
        mHl = value;
    else if (name.compare("muw") == 0)
        muw = value;
    else if (name.compare("SM_M12D") == 0)
        SM_M12D = value;
    else if (name.compare("phiEpsK") == 0)
        phiEpsK = value;
    else if (name.compare("KbarEpsK") == 0)
        KbarEpsK = value;
    else if (name.compare("Dmk") == 0)
        Dmk = value;
    else if (name.compare("DeltaMK") == 0)
        DeltaMK = value;
    else if (name.compare("mneutrino_1") == 0) {
        leptons[NEUTRINO_1].setMass(value);
        computeYn = true;
    } else if (name.compare("mneutrino_2") == 0) {
        leptons[NEUTRINO_2].setMass(value);
        computeYn = true;
    } else if (name.compare("mneutrino_3") == 0) {
        leptons[NEUTRINO_3].setMass(value);
        computeYn = true;
    } else if (name.compare("melectron") == 0) {
        leptons[ELECTRON].setMass(value);
        computeYe = true;
    } else if (name.compare("mmu") == 0) {
        leptons[MU].setMass(value);
        computeYe = true;
    } else if (name.compare("mtau") == 0) {
        leptons[TAU].setMass(value);
        computeYe = true;
    } else if (name.compare("lambda") == 0) {
        lambda = value;
        computeCKM = true;
    } else if (name.compare("A") == 0) {
        A = value;
        computeCKM = true;
    } else if (name.compare("rhob") == 0) {
        rhob = value;
        computeCKM = true;
    } else if (name.compare("etab") == 0) {
        etab = value;
        computeCKM = true;
    } else if (name.compare("delMw") == 0)
        delMw = value;
    else if (name.compare("delSin2th_l") == 0)
        delSin2th_l = value;
    else
        QCD::SetParameter(name, value);
}

bool StandardModel::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NSMvars; i++) {
        if (DPars.find(SMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory SM parameter " << SMvars[i] << std::endl;
            return false;
        }
    }
    return(QCD::CheckParameters(DPars));
}


///////////////////////////////////////////////////////////////////////////
// Flags

bool StandardModel::SetFlag(const std::string name, const bool& value) {  
    bool res = false;
    if (name.compare("FixedAllSMparams") == 0) {
        FlagFixedAllSMparams = value;
        res = true;
    } else if (name.compare("EWCHMN") == 0) {
        FlagEWCHMN = value;
        res = true;
    } else if (name.compare("EWABC") == 0) {
        FlagEWABC = value;
        res = true;
    } else if (name.compare("EWABC2") == 0) {
        FlagEWABC2 = value;
        res = true;
    } else if (name.compare("EWBURGESS") == 0) {
        FlagEWBURGESS = value;
        res = true;
    } else if (name.compare("R0bApproximate") == 0) {
        FlagR0bApproximate = value;
        res = true;
    } else if (name.compare("withoutNonUniversalVCinEpsilons") == 0) {
        FlagWithoutNonUniversalVC = value;
        res = true;
    }
    return(res);
}


///////////////////////////////////////////////////////////////////////////
// Initialization and Matching

bool StandardModel::InitializeModel() {
    myStandardModelMatching = new StandardModelMatching(*this);
    SetModelInitialized(true);
    myEWSM = new EWSM(*this);
    this->SetEWSMflags(*myEWSM);
    return(true);
}

void StandardModel::SetEWSMflags(EWSM& myEWSM) {
    myEWSM.setSchemeMw(EWSM::APPROXIMATEFORMULA);
    myEWSM.setSchemeRhoZ(EWSM::OMSI);
    myEWSM.setSchemeKappaZ(EWSM::APPROXIMATEFORMULA);
}


///////////////////////////////////////////////////////////////////////////

double StandardModel::v() const {
    return ( 1. / sqrt(sqrt(2.) * GF) );
}

double StandardModel::Mw_tree() const {
    double tmp = 4.0*M_PI*ale/sqrt(2.0)/GF/Mz/Mz;
    return ( Mz/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp)) );
}


////////////////////////////////////////////////////////////////////////

double StandardModel::ale_OS(const double mu, orders order) const {
    if (mu < 50.0) 
        throw std::runtime_error("out of range in StandardModel::Als_OS()"); 
    
    double N = 20.0/3.0;
    double beta1 = N/3.0;
    double beta2 = N/4.0;
    double alpha_ini = alphaMz();
    double v = 1.0 + 2.0*beta1*alpha_ini/M_PI*log(Mz/mu);

    switch (order) {
        case LO:
            return ( alpha_ini/v );
        case FULLNLO:
            return ( alpha_ini/v*(1.0 - beta2/beta1*alpha_ini/M_PI*log(v)/v) );
        default:
            throw std::runtime_error("Error in StandardModel::Als_OS()"); 
    }
}


////////////////////////////////////////////////////////////////////////

double StandardModel::Mw0() const {
    return ( sqrt(c02())*Mz );
}

double StandardModel::s02() const {
    double tmp = 1.0 - 4.0*M_PI*alphaMz()/sqrt(2.0)/GF/Mz/Mz;
    if (tmp < 0.0)
        throw std::runtime_error("Error in StandardModel::s02()");
    
    return ( ( 1.0 - sqrt(tmp) )/2.0 );
}

double StandardModel::c02() const {
    return ( 1.0 - s02() );
}


////////////////////////////////////////////////////////////////////////

double StandardModel::DeltaAlphaLepton(const double s) const {
    return myEWSM->DeltaAlphaLepton(s);
}

double StandardModel::DeltaAlphaL5q() const {
    return myEWSM->DeltaAlphaL5q();
}
    
double StandardModel::DeltaAlpha() const {
    return myEWSM->DeltaAlpha();
}
    
double StandardModel::alphaMz() const {
    return myEWSM->alphaMz();
}
    
double StandardModel::Mw() const {
    return myEWSM->Mw_SM();
}

double StandardModel::cW2() const {
    return myEWSM->cW2_SM();
}

double StandardModel::sW2() const {
    return myEWSM->sW2_SM();
}
    
complex StandardModel::rhoZ_l(const lepton l) const {
    return myEWSM->rhoZ_l_SM(l);
}
    
complex StandardModel::rhoZ_q(const quark q) const {
    //if (q!=BOTTOM) // TEST
    if (q!=BOTTOM || (q==BOTTOM && !FlagR0bApproximate) ) 
        return myEWSM->rhoZ_q_SM(q);
        
    // Gamma_u + Gamma_d + Gamma_c + Gamma_s + 4*N_c*Gamma_0*RVh
    double G0 = GF*pow(Mz,3.0)/24.0/sqrt(2.0)/M_PI;    
    double Gamma = 0.0;
    Gamma += 3.0*G0*rhoZ_q(UP).abs()
             * ( (myEWSM->gVq_SM(UP)/myEWSM->gAq_SM(UP)).abs2()*RVq(UP) + RAq(UP) ) + Delta_EWQCD(UP); 
    Gamma += 3.0*G0*rhoZ_q(DOWN).abs()
             * ( (myEWSM->gVq_SM(DOWN)/myEWSM->gAq_SM(DOWN)).abs2()*RVq(DOWN) + RAq(DOWN) ) + Delta_EWQCD(DOWN); 
    Gamma += 3.0*G0*rhoZ_q(CHARM).abs()
             * ( (myEWSM->gVq_SM(CHARM)/myEWSM->gAq_SM(CHARM)).abs2()*RVq(CHARM) + RAq(CHARM) ) + Delta_EWQCD(CHARM); 
    Gamma += 3.0*G0*rhoZ_q(STRANGE).abs()
             * ( (myEWSM->gVq_SM(STRANGE)/myEWSM->gAq_SM(STRANGE)).abs2()*RVq(STRANGE) + RAq(STRANGE) ) + Delta_EWQCD(STRANGE); 
    //Gamma += 4.0*3.0*G0*RVh(); /* RVh depends on rho_Z_q(BOTTOM). */
    
    // |kappaZ_b| from R_0^b
    double R0b = myEWSM->R0_bottom_SM();
    double Qb = getQuarks(BOTTOM).getCharge();  
    double gVb_over_gAb_abs2 = (1.0 - 4.0*fabs(Qb)*myEWSM->kappaZ_q_SM(BOTTOM)*myEWSM->sW2_SM()).abs2();
    double absRhoZb = ( Gamma*R0b/(1.0-R0b) - Delta_EWQCD(BOTTOM) )
                      /( 3.0*G0*( gVb_over_gAb_abs2*RVq(BOTTOM) + RAq(BOTTOM) ) );
    
    // Im(kappaZ_b)
    double ImRhoZb = myEWSM->rhoZ_q_SM(BOTTOM).imag();
    if (absRhoZb < ImRhoZb)
        throw std::runtime_error("Error in StandardModel::rhoZ_q"); 

    return complex(sqrt(absRhoZb*absRhoZb - ImRhoZb*ImRhoZb), ImRhoZb, false);
}
    
complex StandardModel::kappaZ_l(const lepton l) const {
    return myEWSM->kappaZ_l_SM(l);
}
    
complex StandardModel::kappaZ_q(const quark q) const {
    return myEWSM->kappaZ_q_SM(q);
}

complex StandardModel::gVl(const lepton l) const {
    return myEWSM->gVl_SM(l);
}

complex StandardModel::gVq(const quark q) const {
    if (q!=BOTTOM || (q==BOTTOM && !FlagR0bApproximate) ) 
        return myEWSM->gVq_SM(q);

    double Qb = getQuarks(BOTTOM).getCharge();  
    return ( StandardModel::gAq(BOTTOM)
             *(1.0 - 4.0*fabs(Qb)*myEWSM->kappaZ_q_SM(BOTTOM)*myEWSM->sW2_SM()) );
}
    
complex StandardModel::gAl(const lepton l) const {
    return myEWSM->gAl_SM(l);
}

complex StandardModel::gAq(const quark q) const {
    if (q!=BOTTOM || (q==BOTTOM && !FlagR0bApproximate) ) 
        return myEWSM->gAq_SM(q);

    return ( sqrt(StandardModel::rhoZ_q(BOTTOM))*getQuarks(BOTTOM).getIsospin() );
}

double StandardModel::Delta_EWQCD(const StandardModel::quark q) const 
{
    switch(q) {
        case StandardModel::UP:
        case StandardModel::CHARM:
            return ( -0.000113 );
        case StandardModel::TOP:
            return ( 0.0 );
        case StandardModel::DOWN:
        case StandardModel::STRANGE:
            return ( -0.000160 );
        case StandardModel::BOTTOM:
            return ( -0.000040 );
        default:
            throw std::runtime_error("Error in StandardModel::Delta_EWQCD");   
    }
}

double StandardModel::RVq(const StandardModel::quark q) const
{
    if (q==StandardModel::TOP) return 0.0;
    
    double mcMz, mbMz;
    if (!bDebug) {
        mcMz = Mrun(Mz, getQuarks(CHARM).getMass(), FULLNNLO);
        mbMz = Mrun(Mz, getQuarks(BOTTOM).getMass(), FULLNNLO);  
    } else {
        mcMz = 0.56381685; 
        mbMz = 2.8194352;
    }

    /* electric charge squared */
    double Qf2 = pow(getQuarks(q).getCharge(),2.0);

    /* s = Mz^2 */
    double s = Mz*Mz;

    /* products of the charm and bottom masses at Mz */
    double mcMz2 = mcMz*mcMz;
    double mbMz2 = mbMz*mbMz;
    double mqMz2, mqdash4;
    switch(q) {
        case StandardModel::CHARM:
            mqMz2 = mcMz*mcMz;
            mqdash4 = mbMz2*mbMz2;
            break;
        case StandardModel::BOTTOM:
            mqMz2 = mbMz*mbMz;
            mqdash4 = mcMz2*mcMz2;
            break;
        default:
            mqMz2 = 0.0;
            mqdash4 = 0.0;
            break;
    }

    /* Logarithms */
    //double log_t = log(pow(getQuarks(TOP).getMass(),2.0)/s);
    double log_t = log(mtpole*mtpole/s); // the pole mass
    double log_c = log(mcMz2/s);
    double log_b = log(mbMz2/s);
    double log_q;
    switch(q) {
        case StandardModel::CHARM:
        case StandardModel::BOTTOM:
            log_q = log(mqMz2/s);
            break;
        default:
            log_q = 0.0;
            break;
    }    
    
    /* the active number of flavour */
    double nf = 5.0;

    /* zeta functions */
    double zeta2 = gsl_sf_zeta_int(2);
    double zeta3 = gsl_sf_zeta_int(3);
    //double zeta4 = gsl_sf_zeta_int(4);
    double zeta5 = gsl_sf_zeta_int(5);

    /* massless non-singlet corrections */
    double C02 = 365.0/24.0 - 11.0*zeta3 + (-11.0/12.0 + 2.0/3.0*zeta3)*nf;
    double C03 = 87029.0/288.0 - 121.0/8.0*zeta2 - 1103.0/4.0*zeta3
                 + 275.0/6.0*zeta5 
                 + (-7847.0/216.0 + 11.0/6.0*zeta2 + 262.0/9.0*zeta3
                    - 25.0/9.0*zeta5)*nf
                 + (151.0/162.0 - zeta2/18.0 - 19.0/27.0*zeta3)*nf*nf;
    double C04 = -156.61 + 18.77*nf - 0.7974*nf*nf + 0.0215*nf*nf*nf;
    //std::cout << "TEST: C02 = " << C02 << std::endl;// TEST (should be 1.40923)
    //std::cout << "TEST: C03 = " << C03 << std::endl;// TEST (should be -12.7671)
    //std::cout << "TEST: C04 = " << C04 << std::endl;// TEST (should be -80.0075)

    /* quadratic massive corrections */
    double C23  = -80.0 + 60.0*zeta3 + (32.0/9.0 - 8.0/3.0*zeta3)*nf;
    double C21V = 12.0;
    double C22V = 253.0/2.0 - 13.0/3.0*nf;
    double C23V = 2522.0 - 855.0/2.0*zeta2 + 310.0/3.0*zeta3 - 5225.0/6.0*zeta5
                  + (-4942.0/27.0 + 34.0*zeta2 - 394.0/27.0*zeta3
                     + 1045.0/27.0*zeta5)*nf
                  + (125.0/54.0 - 2.0/3.0*zeta2)*nf*nf;

    /* quartic massive corrections */
    double C42  = 13.0/3.0 - 4.0*zeta3;
    //double C40V = -6.0; /* not used */
    double C41V = -22.0;
    double C42V = -3029.0/12.0 + 162.0*zeta2 + 112.0*zeta3
                  + (143.0/18.0 - 4.0*zeta2 - 8.0/3.0*zeta3)*nf;
    double C42VL= -11.0/2.0 + nf/3.0;

    /* power suppressed top-mass correction */
    //double xt = s/pow(getQuarks(TOP).getMass(),2.0);
    double xt = s/mtpole/mtpole; // the pole mass
    double C2t = xt*(44.0/675.0 - 2.0/135.0*(-log_t));

    /* rescaled strong coupling constant */
    double AlsMzPi  = AlsMz/M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;
    
    /* radiator function to the vector current */
    double RVf;
    RVf = 1.0 + 3.0/4.0*Qf2*alphaMz()/M_PI + AlsMzPi - Qf2/4.0*alphaMz()/M_PI*AlsMzPi
            + (C02 + C2t)*AlsMzPi2 + C03*AlsMzPi3 + C04*AlsMzPi4
            + (mcMz2 + mbMz2)/s*C23*AlsMzPi3
            + mqMz2/s*(C21V*AlsMzPi + C22V*AlsMzPi2 + C23V*AlsMzPi3)
            + mcMz2*mcMz2/s/s*(C42 - log_c)*AlsMzPi2
            + mbMz2*mbMz2/s/s*(C42 - log_b)*AlsMzPi2
            + mqMz2*mqMz2/s/s*(C41V*AlsMzPi + (C42V + C42VL*log_q)*AlsMzPi2)
            + 12.0*mqdash4/s/s*AlsMzPi2
            - mqMz2*mqMz2*mqMz2/s/s/s
              *(8.0+16.0/27.0*(155.0 + 6.0*log_q)*AlsMzPi);    
    return RVf;    
}

double StandardModel::RAq(const StandardModel::quark q) const
{
    if (q==StandardModel::TOP) return 0.0;

    double mcMz, mbMz;
    if (!bDebug) {
        mcMz = Mrun(Mz, getQuarks(CHARM).getMass(), FULLNNLO);
        mbMz = Mrun(Mz, getQuarks(BOTTOM).getMass(), FULLNNLO);  
    } else {
        mcMz = 0.56381685; 
        mbMz = 2.8194352;
    }

    /* z-component of isospin */
    double I3q = getQuarks(q).getIsospin();
    /* electric charge squared */
    double Qf2 = pow(getQuarks(q).getCharge(),2.0);

    /* s = Mz^2 */
    double s = Mz*Mz;

    /* products of the charm and bottom masses at Mz */
    double mcMz2 = mcMz*mcMz;
    double mbMz2 = mbMz*mbMz;
    double mqMz2, mqdash4;
    switch(q) {
        case StandardModel::CHARM:
            mqMz2 = mcMz*mcMz;
            mqdash4 = mbMz2*mbMz2;
            break;
        case StandardModel::BOTTOM:
            mqMz2 = mbMz*mbMz;
            mqdash4 = mcMz2*mcMz2;
            break;
        default:
            mqMz2 = 0.0;
            mqdash4 = 0.0;
            break;
    }

    /* Logarithms */
    //double log_t = log(pow(getQuarks(TOP).getMass(),2.0)/s);
    double log_t = log(mtpole*mtpole/s); // the pole mass
    double log_c = log(mcMz2/s);
    double log_b = log(mbMz2/s);
    double log_q;
    switch(q) {
        case StandardModel::CHARM:
        case StandardModel::BOTTOM:
            log_q = log(mqMz2/s);
            break;
        default:
            log_q = 0.0;
            break;
    }    
    
    /* the active number of flavour */
    double nf = 5.0;

    /* zeta functions */
    double zeta2 = gsl_sf_zeta_int(2);
    double zeta3 = gsl_sf_zeta_int(3);
    double zeta4 = gsl_sf_zeta_int(4);
    double zeta5 = gsl_sf_zeta_int(5);

    /* massless non-singlet corrections */
    double C02 = 365.0/24.0 - 11.0*zeta3 + (-11.0/12.0 + 2.0/3.0*zeta3)*nf;
    double C03 = 87029.0/288.0 - 121.0/8.0*zeta2 - 1103.0/4.0*zeta3
                 + 275.0/6.0*zeta5 
                 + (-7847.0/216.0 + 11.0/6.0*zeta2 + 262.0/9.0*zeta3
                    - 25.0/9.0*zeta5)*nf
                 + (151.0/162.0 - zeta2/18.0 - 19.0/27.0*zeta3)*nf*nf;
    double C04 = -156.61 + 18.77*nf - 0.7974*nf*nf + 0.0215*nf*nf*nf;
    //std::cout << "TEST: C02 = " << C02 << std::endl;// TEST (should be 1.40923)
    //std::cout << "TEST: C03 = " << C03 << std::endl;// TEST (should be -12.7671)
    //std::cout << "TEST: C04 = " << C04 << std::endl;// TEST (should be -80.0075)

    /* quadratic massive corrections */
    double C23  = -80.0 + 60.0*zeta3 + (32.0/9.0 - 8.0/3.0*zeta3)*nf;
    double C20A = -6.0;
    double C21A = -22.0;
    double C22A = -8221.0/24.0 + 57.0*zeta2 + 117.0*zeta3
                  + (151.0/12.0 - 2.0*zeta2 - 4.0*zeta3)*nf;
    double C23A = -4544045.0/864.0 + 1340.0*zeta2 + 118915.0/36.0*zeta3
                  - 127.0*zeta5
                  + (71621.0/162.0 - 209.0/2.0*zeta2 - 216.0*zeta3
                     + 5.0*zeta4 + 55.0*zeta5)*nf
                  + (-13171.0/1944.0 + 16.0/9.0*zeta2 + 26.0/9.0*zeta3)*nf*nf;

    /* quartic massive corrections */
    double C42  = 13.0/3.0 - 4.0*zeta3;
    double C40A = 6.0;
    double C41A = 10.0;
    double C42A = 3389.0/12.0 - 162.0*zeta2 - 220.0*zeta3
                  + (-41.0/6.0 + 4.0*zeta2 + 16.0/3.0*zeta3)*nf;
    double C42AL= 77.0/2.0 - 7.0/3.0*nf;

    /* power suppressed top-mass correction */
    //double xt = s/pow(getQuarks(TOP).getMass(),2.0);
    double xt = s/mtpole/mtpole; // the pole mass
    double C2t = xt*(44.0/675.0 - 2.0/135.0*(-log_t));

    /* singlet axial-vector corrections */
    double I2 = -37.0/12.0 + (-log_t) + 7.0/81.0*xt + 0.0132*xt*xt;
    double I3 = -5075.0/216.0 + 23.0/6.0*zeta2 + zeta3 + 67.0/18.0*(-log_t)
                + 23.0/12.0*log_t*log_t;
    double I4 = 49.0309 - 17.6637*(-log_t) + 14.6597*log_t*log_t 
                + 3.6736*(-log_t*log_t*log_t);
    
    /* rescaled strong coupling constant */
    double AlsMzPi  = AlsMz/M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;    
    
    /* radiator function to the axial-vector current */
    double RAf;
    RAf = 1.0 + 3.0/4.0*Qf2*alphaMz()/M_PI + AlsMzPi - Qf2/4.0*alphaMz()/M_PI*AlsMzPi
            + (C02 + C2t - 2.0*I3q*I2)*AlsMzPi2
            + (C03 - 2.0*I3q*I3)*AlsMzPi3
            + (C04 - 2.0*I3q*I4)*AlsMzPi4
            + (mcMz2 + mbMz2)/s*C23*AlsMzPi3
            + mqMz2/s*(C20A + C21A*AlsMzPi + C22A*AlsMzPi2
                       + 6.0*(3.0 + log_t)*AlsMzPi2 + C23A*AlsMzPi3)
            //- 10.0*mqMz2/pow(getQuarks(TOP).getMass(),2.0)
            - 10.0*mqMz2/mtpole/mtpole // the pole mass
              *(8.0/81.0 + log_t/54.0)*AlsMzPi2
            + mcMz2*mcMz2/s/s*(C42 - log_c)*AlsMzPi2
            + mbMz2*mbMz2/s/s*(C42 - log_b)*AlsMzPi2
            + mqMz2*mqMz2/s/s*(C40A + C41A*AlsMzPi
                               + (C42A + C42AL*log_q)*AlsMzPi2)
            - 12.0*mqdash4/s/s*AlsMzPi2 ;  
    return RAf;
}

double StandardModel::RVh() const
{
    /* rescaled strong coupling constant */
    double AlsMzPi  = AlsMz/M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;

    complex gV_sum(0.0, 0.0); 
    complex gV_q;
    for (int q=0; q<6; q++) {
        gV_q = gVq((StandardModel::quark)q);
        if (q==(int)(StandardModel::TOP)) 
            gV_q = 0.0;
        gV_sum += gV_q;
    }
    
    // singlet vector corrections
    return ( gV_sum.abs2()*(-0.4132*AlsMzPi3 - 4.9841*AlsMzPi4) );
}

double StandardModel::GammaW() const {
    return myEWSM->GammaW_SM();
}

double StandardModel::epsilon1_SM() const {
    double rhoZe = myEWSM->rhoZ_l_SM(ELECTRON).real();
    double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );

    return DeltaRhoPrime;
}

double StandardModel::epsilon2_SM() const {
    double Qe = getLeptons(ELECTRON).getCharge();
    double s_W2 = myEWSM->sW2_SM(), c_W2 = myEWSM->cW2_SM();
    double rhoZe = myEWSM->rhoZ_l_SM(ELECTRON).real();
    double sin2thetaEff = myEWSM->kappaZ_l_SM(ELECTRON).real()*s_W2;
    double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    double DeltaKappaPrime = sin2thetaEff/s02() - 1.0;
    double DeltaRW = 1.0 - M_PI*alphaMz()/(sqrt(2.0)*GF*Mz*Mz*s_W2*c_W2);
    
    return ( c02()*DeltaRhoPrime + s02()*DeltaRW/(c02() - s02()) 
             - 2.0*s02()*DeltaKappaPrime );
}

double StandardModel::epsilon3_SM() const {
    double Qe = getLeptons(ELECTRON).getCharge();
    double rhoZe = myEWSM->rhoZ_l_SM(ELECTRON).real();
    double sin2thetaEff = myEWSM->kappaZ_l_SM(ELECTRON).real()*myEWSM->sW2_SM();
    double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    double DeltaKappaPrime = sin2thetaEff/s02() - 1.0;

    return ( c02()*DeltaRhoPrime + (c02() - s02())*DeltaKappaPrime );
}

double StandardModel::epsilonb_SM() const {
    /* epsilon_b from g_A^b
     * see Eq.(13) of IJMP A7, 1031 (1998) by Altarelli et al. */
    //double rhoZe = myEWSM->rhoZ_l_SM(ELECTRON).real();
    //double rhoZb = myEWSM->rhoZ_q_SM(BOTTOM).real();
    //double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    //double eps1 = DeltaRhoPrime;
    //return ( - 1.0 + sqrt(rhoZb)/(1.0 + eps1/2.0) );

    /* epsilon_b from Re(g_V^b/g_A^b), i.e. Re(kappaZ_b)
     * see Eq.(13) of IJMP A7, 1031 (1998) by Altarelli et al. */
    complex kappaZe = myEWSM->kappaZ_l_SM(ELECTRON);
    complex kappaZb = myEWSM->kappaZ_q_SM(BOTTOM);
    if (FlagWithoutNonUniversalVC) 
        return ( kappaZe.real()/kappaZb.real() - 1.0 ); 
    else 
        return ( (kappaZe.real() + myEWSM->kappaZ_q_SM_FlavorDep(BOTTOM).real())
                 /kappaZb.real() - 1.0 );     
    
    /* epsilon_b from Gamma_b via Eqs.(11), (12) and (16) of IJMP A7, 
     * 1031 (1998) by Altarelli et al. 
     * Note: mb has to be mb=4.7, since Eq.(16) were derived with this value. 
     */
    //double als_Mz = Als(Mz, FULLNNLO);
    //double delta_als = (als_Mz - 0.119)/M_PI;
    //double delta_alpha = (alphaMz() - 1.0/128.90)/ale;
    //double Gamma_b_Born = 0.3798*( 1.0 + delta_als - 0.42*delta_alpha);
    //double a = als_Mz/M_PI;
    //double RQCD = 1.0 + 1.2*a - 1.1*a*a - 13.0*a*a*a;
    //double mb = Mrun(Mz, getQuarks(BOTTOM).getMass(), FULLNNLO);// This is wrong!
    //double mb = 4.7;
    //std::cout << "mb = " << mb << std::endl;
    //double beta = sqrt(1.0 - 4.0*mb*mb/Mz/Mz);
    //double Nc = 3.0; 
    //double factor = GF*Mz*Mz*Mz/6.0/M_PI/sqrt(2.0);
    //double Gamma_b = factor*beta*((3.0 - beta*beta)/2.0*myEWSM->gVq_SM(BOTTOM).abs2()
    //                              + beta*beta*myEWSM->gAq_SM(BOTTOM).abs2())
    //                 *Nc*RQCD*(1.0 + alphaMz()/12.0/M_PI);
    //return ( (Gamma_b/Gamma_b_Born - 1.0 - 1.42*epsilon1_SM() 
    //          + 0.54*epsilon3_SM() )/2.29 );
}

double StandardModel::epsilon1() const{ 
    return epsilon1_SM();
}

double StandardModel::epsilon2() const {
    return epsilon2_SM();    
}
    
double StandardModel::epsilon3() const {
    return epsilon3_SM();
}

double StandardModel::epsilonb() const {
    return epsilonb_SM();    
}

double StandardModel::deltaGVb() const 
{
    return 0.0;
}

double StandardModel::deltaGAb() const
{
    return 0.0;
}


////////////////////////////////////////////////////////////////////////
// CKM parameters

// Angles

double StandardModel::getBeta() const {
    return (-VCKM(1, 0) * VCKM(1, 2).conjugate() / (VCKM(2, 0) * VCKM(2, 2).conjugate())).arg();
}

double StandardModel::getGamma() const {
    return (-VCKM(0, 0) * VCKM(0, 2).conjugate() / (VCKM(1, 0) * VCKM(1, 2).conjugate())).arg();
}

double StandardModel::getAlpha() const {
    return (-VCKM(2, 0) * VCKM(2, 2).conjugate() / (VCKM(0, 0) * VCKM(0, 2).conjugate())).arg();
}

double StandardModel::getBetas() const {
    return (-VCKM(2, 1) * VCKM(2, 2).conjugate() / (VCKM(1, 1) * VCKM(1, 2).conjugate())).arg();
}

// Lambda_q

complex StandardModel::getlamt() const {
    return VCKM(2, 0) * VCKM(2, 1).conjugate();
}

complex StandardModel::getlamc() const {
    return VCKM(1, 0) * VCKM(1, 1).conjugate();
}

complex StandardModel::getlamu() const {
    return VCKM(0, 0) * VCKM(0, 1).conjugate();
}

complex StandardModel::getlamt_d() const {
    return VCKM(2, 0) * VCKM(2, 2).conjugate();
}

complex StandardModel::getlamc_d() const {
    return VCKM(1, 0) * VCKM(1, 2).conjugate();
}

complex StandardModel::getlamu_d() const {
    return VCKM(0, 0) * VCKM(0, 2).conjugate();
}

complex StandardModel::getlamt_s() const {
    return VCKM(2, 1) * VCKM(2, 2).conjugate();
}

complex StandardModel::getlamc_s() const {
    return VCKM(1, 1) * VCKM(1, 2).conjugate();
}

complex StandardModel::getlamu_s() const {
    return VCKM(0, 1) * VCKM(0, 2).conjugate();
}

