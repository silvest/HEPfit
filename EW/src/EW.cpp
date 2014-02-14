/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <EWSM.h>
#include <NPEpsilons.h>
#include <NPSTU.h>
#include <NPZbbbar.h>
#include <ThObservable.h> // for ThObservable::GeVminus2_to_nb
#include "EW.h"


EW::EW(const StandardModel& SM_i) 
: ThObsType(SM_i), SM(SM_i), myEW_NPZff(SM_i)
{
}


////////////////////////////////////////////////////////////////////////

bool EW::checkNPZff_linearized() const
{
    std::string Model = SM.ModelName();
    if (Model.compare("NPZbbbar") == 0) {
        if ((static_cast<const NPZbbbar*> (&SM))->IsFlagNotLinearizedNP())
            return false;
        else
            return true;
    } else if (Model.compare("StandardModel") == 0
            || Model.compare("NPHiggsST") == 0
            || Model.compare("NPSTU") == 0
            || Model.compare("NPSTUVWXY") == 0
            || Model.compare("NPEpsilons_pureNP") == 0
            || Model.compare("NPEffective1") == 0
            || Model.compare("NPEffective2") == 0)
        return true;
    else if (Model.compare("NPEpsilons") == 0)
        return false;
    else 
        throw std::runtime_error("EW::checkNPZff_linearized(): " + Model
                                 + " cannot be dealt with this function");
}


////////////////////////////////////////////////////////////////////////

double EW::Delta_EWQCD(const QCD::quark q) const
{
    switch(q) {
        case QCD::UP:
        case QCD::CHARM:
            return ( -0.000113 );
        case QCD::TOP:
            return ( 0.0 );
        case QCD::DOWN:
        case QCD::STRANGE:
            return ( -0.000160 );
        case QCD::BOTTOM:
            return ( -0.000040 );
        default:
            throw std::runtime_error("Error in EWSM::Delta_EWQCD");
    }
}


double EW::RVq(const QCD::quark q) const
{
    if (q==QCD::TOP) return 0.0;

    double mcMz, mbMz;
    mcMz = SM.Mrun(SM.getMz(), SM.getQuarks(SM.CHARM).getMass(), FULLNNLO);
    mbMz = SM.Mrun(SM.getMz(), SM.getQuarks(SM.BOTTOM).getMass(), FULLNNLO);  
    //mcMz = 0.56381685; /* for debug */
    //mbMz = 2.8194352; /* for debug */

    double MtPole = SM.getMtpole();

    /* electric charge squared */
    double Qf2 = pow(SM.getQuarks(q).getCharge(),2.0);

    /* s = Mz^2 */
    double s = SM.getMz()*SM.getMz();

    /* products of the charm and bottom masses at Mz */
    double mcMz2 = mcMz*mcMz;
    double mbMz2 = mbMz*mbMz;
    double mqMz2, mqdash4;
    switch(q) {
        case QCD::CHARM:
            mqMz2 = mcMz*mcMz;
            mqdash4 = mbMz2*mbMz2;
            break;
        case QCD::BOTTOM:
            mqMz2 = mbMz*mbMz;
            mqdash4 = mcMz2*mcMz2;
            break;
        default:
            mqMz2 = 0.0;
            mqdash4 = 0.0;
            break;
    }

    /* Logarithms */
    //double log_t = log(pow(SM.getQuarks(TOP).getMass(),2.0)/s);
    double log_t = log(MtPole*MtPole/s); // the pole mass
    double log_c = log(mcMz2/s);
    double log_b = log(mbMz2/s);
    double log_q;
    switch(q) {
        case QCD::CHARM:
        case QCD::BOTTOM:
            log_q = log(mqMz2/s);
            break;
        default:
            log_q = 0.0;
            break;
    }

    /* the active number of flavour */
    double nf = 5.0;

    /* zeta functions */
    double zeta2 = SM.getEWSM()->getMyCache()->getZeta2();
    double zeta3 = SM.getEWSM()->getMyCache()->getZeta3();
    //double zeta4 = SM.getEWSM()->getMyCache()->GetZeta4();
    double zeta5 = SM.getEWSM()->getMyCache()->getZeta5();

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
    double C40V = -6.0;
    double C41V = -22.0;
    double C42V = -3029.0/12.0 + 162.0*zeta2 + 112.0*zeta3
                  + (143.0/18.0 - 4.0*zeta2 - 8.0/3.0*zeta3)*nf;
    double C42VL= -11.0/2.0 + nf/3.0;

    /* power suppressed top-mass correction */
    //double xt = s/pow(getQuarks(TOP).getMass(),2.0);
    double xt = s/MtPole/MtPole; // the pole mass
    double C2t = xt*(44.0/675.0 - 2.0/135.0*(-log_t));

    /* rescaled strong coupling constant */
    double AlsMzPi  = SM.getAlsMz()/M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;

    /* electromagnetic coupling at Mz */
    double alpMz = SM.alphaMz();

    /* radiator function to the vector current */
    double RVf;
    RVf = 1.0 + 3.0/4.0*Qf2*alpMz/M_PI + AlsMzPi - Qf2/4.0*alpMz/M_PI*AlsMzPi
            + (C02 + C2t)*AlsMzPi2 + C03*AlsMzPi3 + C04*AlsMzPi4
            + (mcMz2 + mbMz2)/s*C23*AlsMzPi3
            + mqMz2/s*(C21V*AlsMzPi + C22V*AlsMzPi2 + C23V*AlsMzPi3)
            + mcMz2*mcMz2/s/s*(C42 - log_c)*AlsMzPi2
            + mbMz2*mbMz2/s/s*(C42 - log_b)*AlsMzPi2
            + mqMz2*mqMz2/s/s*(C40V + C41V*AlsMzPi + (C42V + C42VL*log_q)*AlsMzPi2)
            + 12.0*mqdash4/s/s*AlsMzPi2
            - mqMz2*mqMz2*mqMz2/s/s/s
              *(8.0+16.0/27.0*(155.0 + 6.0*log_q)*AlsMzPi);
    return RVf;
}


double EW::RAq(const QCD::quark q) const
{
    if (q==QCD::TOP) return 0.0;

    double mcMz, mbMz;
    mcMz = SM.Mrun(SM.getMz(), SM.getQuarks(SM.CHARM).getMass(), FULLNNLO);
    mbMz = SM.Mrun(SM.getMz(), SM.getQuarks(SM.BOTTOM).getMass(), FULLNNLO);
    //mcMz = 0.56381685; /* for debug */
    //mbMz = 2.8194352; /* for debug */

    double MtPole = SM.getMtpole();

    /* z-component of isospin */
    double I3q = SM.getQuarks(q).getIsospin();
    /* electric charge squared */
    double Qf2 = pow(SM.getQuarks(q).getCharge(),2.0);

    /* s = Mz^2 */
    double s = SM.getMz()*SM.getMz();

    /* products of the charm and bottom masses at Mz */
    double mcMz2 = mcMz*mcMz;
    double mbMz2 = mbMz*mbMz;
    double mqMz2, mqdash4;
    switch(q) {
        case QCD::CHARM:
            mqMz2 = mcMz*mcMz;
            mqdash4 = mbMz2*mbMz2;
            break;
        case QCD::BOTTOM:
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
    double log_t = log(MtPole*MtPole/s); // the pole mass
    double log_c = log(mcMz2/s);
    double log_b = log(mbMz2/s);
    double log_q;
    switch(q) {
        case QCD::CHARM:
        case QCD::BOTTOM:
            log_q = log(mqMz2/s);
            break;
        default:
            log_q = 0.0;
            break;
    }

    /* the active number of flavour */
    double nf = 5.0;

    /* zeta functions */
    double zeta2 = SM.getEWSM()->getMyCache()->getZeta2();
    double zeta3 = SM.getEWSM()->getMyCache()->getZeta3();
    double zeta4 = SM.getEWSM()->getMyCache()->getZeta4();
    double zeta5 = SM.getEWSM()->getMyCache()->getZeta5();

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
    double xt = s/MtPole/MtPole; // the pole mass
    double C2t = xt*(44.0/675.0 - 2.0/135.0*(-log_t));

    /* singlet axial-vector corrections */
    double I2 = -37.0/12.0 + (-log_t) + 7.0/81.0*xt + 0.0132*xt*xt;
    double I3 = -5075.0/216.0 + 23.0/6.0*zeta2 + zeta3 + 67.0/18.0*(-log_t)
                + 23.0/12.0*log_t*log_t;
    double I4 = 49.0309 - 17.6637*(-log_t) + 14.6597*log_t*log_t
                + 3.6736*(-log_t*log_t*log_t);

    /* rescaled strong coupling constant */
    double AlsMzPi  = SM.getAlsMz()/M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;

    /* electromagnetic coupling at Mz */
    double alpMz = SM.alphaMz();

    /* radiator function to the axial-vector current */
    double RAf;
    RAf = 1.0 + 3.0/4.0*Qf2*alpMz/M_PI + AlsMzPi - Qf2/4.0*alpMz/M_PI*AlsMzPi
            + (C02 + C2t - 2.0*I3q*I2)*AlsMzPi2
            + (C03 - 2.0*I3q*I3)*AlsMzPi3
            + (C04 - 2.0*I3q*I4)*AlsMzPi4
            + (mcMz2 + mbMz2)/s*C23*AlsMzPi3
            + mqMz2/s*(C20A + C21A*AlsMzPi + C22A*AlsMzPi2
                       + 6.0*(3.0 + log_t)*AlsMzPi2 + C23A*AlsMzPi3)
            //- 10.0*mqMz2/pow(getQuarks(TOP).getMass(),2.0)
            - 10.0*mqMz2/MtPole/MtPole // the pole mass
              *(8.0/81.0 + log_t/54.0)*AlsMzPi2
            + mcMz2*mcMz2/s/s*(C42 - log_c)*AlsMzPi2
            + mbMz2*mbMz2/s/s*(C42 - log_b)*AlsMzPi2
            + mqMz2*mqMz2/s/s*(C40A + C41A*AlsMzPi
                               + (C42A + C42AL*log_q)*AlsMzPi2)
            - 12.0*mqdash4/s/s*AlsMzPi2 ;
    return RAf;
}


double EW::RVh() const
{
    /* rescaled strong coupling constant */
    double AlsMzPi  = SM.getAlsMz()/M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;

    complex gV_sum(0.0, 0.0);
    complex gV_q;
    for (int q=0; q<6; q++) {
        gV_q = SM.getEWSM()->gVq_SM((QCD::quark)q);
        if (q==(int)(QCD::TOP))
            gV_q = 0.0;
        gV_sum += gV_q;
    }

    // singlet vector corrections
    return ( gV_sum.abs2()*(-0.4132*AlsMzPi3 - 4.9841*AlsMzPi4) );
}


////////////////////////////////////////////////////////////////////////

double EW::A_l(const StandardModel::lepton l) const
{
    double Re_gV_over_gA;
    if (checkNPZff_linearized() && SM.getFlagKappaZ().compare("APPROXIMATEFORMULA") == 0) {
        /* SM contribution with the approximate formula */
        double sin2thEff = SM.getEWSM()->getMyApproximateFormulae()->sin2thetaEff_l(l);
        Re_gV_over_gA = 1.0 - 4.0*fabs(SM.getLeptons(l).getCharge())*sin2thEff;
    } else
        Re_gV_over_gA = (SM.getEWSM()->gVl(l)/SM.getEWSM()->gAl(l)).real();
    return ( 2.0*Re_gV_over_gA/(1.0+pow(Re_gV_over_gA,2.0)) );
}


double EW::A_q(const QCD::quark q) const
{
    double Re_gV_over_gA;
    if (checkNPZff_linearized() && SM.getFlagKappaZ().compare("APPROXIMATEFORMULA") == 0) {
        /* SM contribution with the approximate formula */
        double sin2thEff = SM.getEWSM()->getMyApproximateFormulae()->sin2thetaEff_q(q);
        Re_gV_over_gA = 1.0 - 4.0*fabs(SM.getQuarks(q).getCharge())*sin2thEff;
    } else
        Re_gV_over_gA = (SM.getEWSM()->gVq(q)/SM.getEWSM()->gAq(q)).real();
    return ( 2.0*Re_gV_over_gA/(1.0+pow(Re_gV_over_gA,2.0)) );
}


double EW::sin2thetaEff(const StandardModel::lepton l) const 
{
    if (checkNPZff_linearized() && SM.getFlagKappaZ().compare("APPROXIMATEFORMULA") == 0)
        /* SM contribution with the approximate formula */
        return SM.getEWSM()->getMyApproximateFormulae()->sin2thetaEff_l(l);
    else {
        double Re_kappa = SM.getEWSM()->kappaZ_l(l).real();
        return ( Re_kappa*SM.sW2() );
    }
}


double EW::sin2thetaEff(const QCD::quark q) const
{
    if (checkNPZff_linearized() && SM.getFlagKappaZ().compare("APPROXIMATEFORMULA") == 0)
        /* SM contribution with the approximate formula */
        return SM.getEWSM()->getMyApproximateFormulae()->sin2thetaEff_q(q);
    else {
        double Re_kappa = SM.getEWSM()->kappaZ_q(q).real();
        return ( Re_kappa*SM.sW2() );
    }
}


double EW::Gamma_l(const StandardModel::lepton l) const 
{
    double Gamma;
    if (checkNPZff_linearized() && !SM.IsFlagNoApproximateGammaZ()) {
        /* SM contribution with the approximate formula */
        switch (l) {
            case StandardModel::NEUTRINO_1:
            case StandardModel::NEUTRINO_2:
            case StandardModel::NEUTRINO_3:
                Gamma = SM.getEWSM()->getMyApproximateFormulae()->X_extended("Gamma_nu");
                break;
            case StandardModel::ELECTRON:
            case StandardModel::MU:
                Gamma = SM.getEWSM()->getMyApproximateFormulae()->X_extended("Gamma_e_mu");
                break;
            case StandardModel::TAU:
                Gamma = SM.getEWSM()->getMyApproximateFormulae()->X_extended("Gamma_tau");
                break;
            default:
                throw std::runtime_error("Error in EWSM::Gamma_l()");
        }
    } else {
        complex rhoZ_l = SM.getEWSM()->rhoZ_l(l);
        complex gV_over_gA = SM.getEWSM()->gVl(l)/SM.getEWSM()->gAl(l);
        double alphaMz = SM.alphaMz();
        double Q = SM.getLeptons(l).getCharge();
        double xl = pow(SM.getLeptons(l).getMass()/SM.getMz(), 2.0);
        double G0 = SM.getGF()*pow(SM.getMz(),3.0)/24.0/sqrt(2.0)/M_PI;
        Gamma = G0*rhoZ_l.abs()*sqrt(1.0 - 4.0*xl)
                       * ( (1.0 + 2.0*xl)*(gV_over_gA.abs2() + 1.0) - 6.0*xl )
                       * ( 1.0 + 3.0/4.0*alphaMz/M_PI*pow(Q,2.0) );
    }

    return Gamma;
}


double EW::Gamma_q(const QCD::quark q) const 
{
    if (q==QCD::TOP) return 0.0;

    double Gamma;
    if (checkNPZff_linearized() && !SM.IsFlagNoApproximateGammaZ()) {
        /* SM contribution with the approximate formula */
        switch (q) {
            case QCD::UP:
                Gamma = SM.getEWSM()->getMyApproximateFormulae()->X_extended("Gamma_u");
                break;
            case QCD::CHARM:
                Gamma = SM.getEWSM()->getMyApproximateFormulae()->X_extended("Gamma_c");
                break;
            case QCD::DOWN:
            case QCD::STRANGE:
                Gamma = SM.getEWSM()->getMyApproximateFormulae()->X_extended("Gamma_d_s");
                break;
            case QCD::BOTTOM:
                Gamma = SM.getEWSM()->getMyApproximateFormulae()->X_extended("Gamma_b");
                break;
            default:
                throw std::runtime_error("Error in EWSM::Gamma_q()");
        }
    } else {
        complex rhoZ_q = SM.getEWSM()->rhoZ_q(q);
        complex gV_over_gA = SM.getEWSM()->gVq(q)/SM.getEWSM()->gAq(q);
        double G0 = SM.getGF()*pow(SM.getMz(),3.0)/24.0/sqrt(2.0)/M_PI;
        Gamma = 3.0*G0*rhoZ_q.abs()*( gV_over_gA.abs2()*RVq(q) + RAq(q) );

        /* Nonfactorizable EW-QCD corrections */
        Gamma += Delta_EWQCD(q);
    }

    return Gamma;
}


double EW::Gamma_inv() const 
{
    return ( Gamma_l(SM.NEUTRINO_1) + Gamma_l(SM.NEUTRINO_2) 
             + Gamma_l(SM.NEUTRINO_3) );
}


double EW::Gamma_had() const 
{
    double Gamma_had_tmp = Gamma_q(SM.UP) + Gamma_q(SM.DOWN) + Gamma_q(SM.CHARM)
                           + Gamma_q(SM.STRANGE) + Gamma_q(SM.BOTTOM);

    /* Singlet vector contribution (not included in the approximate formula) */
    double G0 = SM.getGF()*pow(SM.getMz(),3.0)/24.0/sqrt(2.0)/M_PI; 
    Gamma_had_tmp += 4.0*3.0*G0*RVh();

    return Gamma_had_tmp;    
}


double EW::Gamma_Z() const 
{
    if (checkNPZff_linearized() && !SM.IsFlagNoApproximateGammaZ())
        /* SM contribution with the approximate formula */
        return SM.getEWSM()->getMyApproximateFormulae()->X_extended("GammaZ");
    else
        return ( Gamma_l(SM.ELECTRON) + Gamma_l(SM.MU) + Gamma_l(SM.TAU)
                 + Gamma_inv() + Gamma_had() );
}


double EW::sigma0_had() const 
{
    if (checkNPZff_linearized() && !SM.IsFlagNoApproximateGammaZ())
        /* SM contribution with the approximate formula */
        return (SM.getEWSM()->getMyApproximateFormulae()->X_extended("sigmaHadron")
                /ThObservable::GeVminus2_to_nb);
    else
        return (12.0*M_PI*Gamma_l(SM.ELECTRON)*Gamma_had()
                /SM.getMz()/SM.getMz()/Gamma_Z()/Gamma_Z());
}


double EW::R0_l() const
{
    if (checkNPZff_linearized() && !SM.IsFlagNoApproximateGammaZ())
        /* SM contribution with the approximate formula */
        return (SM.getEWSM()->getMyApproximateFormulae()->X_extended("R0_lepton"));
    else
        return (Gamma_had()/Gamma_l(SM.ELECTRON));
}


double EW::R0_c() const
{
    if (checkNPZff_linearized() && !SM.IsFlagNoApproximateGammaZ())
        /* SM contribution with the approximate formula */
        return (SM.getEWSM()->getMyApproximateFormulae()->X_extended("R0_charm"));
    else
        return (Gamma_q(SM.CHARM)/Gamma_had());
}


double EW::R0_b() const
{
    if (checkNPZff_linearized() && !SM.IsFlagNoApproximateGammaZ())
        /* SM contribution with the approximate formula */
        return (SM.getEWSM()->getMyApproximateFormulae()->X_extended("R0_bottom"));
    else
        return (Gamma_q(SM.BOTTOM)/Gamma_had());
}


