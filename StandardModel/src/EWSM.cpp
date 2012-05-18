/* 
 * File:   EWSM.cpp
 * Author: mishima
 */

#include <cmath>
#include "EWSM.h"
#include "StandardModel.h"


const double EWSM::Mw_error = 0.00001; /* 0.01 MeV */ 


EWSM::EWSM(const StandardModel& SM_i) : SM(SM_i) {
    flag_order[EW1] = true;
    flag_order[EW1QCD1] = true;
    flag_order[EW1QCD2] = true;
    flag_order[EW2] = true;
    flag_order[EW2QCD1] = true;
    flag_order[EW3] = true;
            
    //schemeMw = OMSI;// for test
    schemeMw = APPROXIMATEFORMULA;    
    schemeRhoZ = OMSI;
    //schemeKappaZ = OMSI;// for test
    schemeKappaZ = APPROXIMATEFORMULA;
    
    myCache = new EWSMcache(SM);
    myOneLoopEW = new EWSMOneLoopEW(*myCache);
    myTwoLoopQCD = new EWSMTwoLoopQCD(*myCache);
    myThreeLoopQCD = new EWSMThreeLoopQCD(*myCache);
    myTwoLoopEW = new EWSMTwoLoopEW(*myCache);
    myThreeLoopEW2QCD = new EWSMThreeLoopEW2QCD(*myCache);
    myThreeLoopEW = new EWSMThreeLoopEW(*myCache);
    myApproximateFormulae = new EWSMApproximateFormulae(SM);    
}


////////////////////////////////////////////////////////////////////////

double EWSM::DeltaAlphaLepton() const {
    double DeltaAlphaL = 0.0;
    if (flag_order[EW1]) 
        DeltaAlphaL += myOneLoopEW->DeltaAlpha_l();
    if (flag_order[EW1QCD1]) 
        DeltaAlphaL += myTwoLoopQCD->DeltaAlpha_l();
    if (flag_order[EW1QCD2]) 
        DeltaAlphaL += myThreeLoopQCD->DeltaAlpha_l();
    if (flag_order[EW2]) 
        DeltaAlphaL += myTwoLoopEW->DeltaAlpha_l();
    if (flag_order[EW2QCD1]) 
        DeltaAlphaL += myThreeLoopEW2QCD->DeltaAlpha_l();
    if (flag_order[EW3]) 
        DeltaAlphaL += myThreeLoopEW->DeltaAlpha_l();
    return DeltaAlphaL;
}


double EWSM::DeltaAlphaL5q() const {
    return (DeltaAlphaLepton() + SM.getDAle5Mz());
}


double EWSM::DeltaAlpha() const {
    // leptonic + hadronic contributions
    double DeltaAlpha = DeltaAlphaL5q(); 

    // Top-quark contribution
    if (flag_order[EW1]) 
        DeltaAlpha += myOneLoopEW->DeltaAlpha_t();
    if (flag_order[EW1QCD1]) 
        DeltaAlpha += myTwoLoopQCD->DeltaAlpha_t();
    if (flag_order[EW1QCD2]) 
        DeltaAlpha += myThreeLoopQCD->DeltaAlpha_t();
    if (flag_order[EW2]) 
        DeltaAlpha += myTwoLoopEW->DeltaAlpha_t();
    if (flag_order[EW2QCD1]) 
        DeltaAlpha += myThreeLoopEW2QCD->DeltaAlpha_t();
    if (flag_order[EW3]) 
        DeltaAlpha += myThreeLoopEW->DeltaAlpha_t();
    return DeltaAlpha; 
}


double EWSM::alphaMz() const {
    return (SM.getAle()/(1.0 - DeltaAlpha()));
}


////////////////////////////////////////////////////////////////////////

double EWSM::Mw() const {
    return Mw_SM();
}


double EWSM::cW2() const {
    return cW2_SM();
}


double EWSM::sW2() const {
    return sW2_SM();    
}

    
complex EWSM::rhoZ_l(const StandardModel::lepton l) const {
    return rhoZ_l_SM(l);
}


complex EWSM::rhoZ_q(const StandardModel::quark q) const {
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));
    return rhoZ_q_SM(q);
}
  

complex EWSM::kappaZ_l(const StandardModel::lepton l) const {
    return kappaZ_l_SM(l);
}


complex EWSM::kappaZ_q(const StandardModel::quark q) const {
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));
    return kappaZ_q_SM(q);
}


complex EWSM::gVl(const StandardModel::lepton l) const {
    return gVl_SM(l);
}


complex EWSM::gVq(const StandardModel::quark q) const {
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));
    return gVq_SM(q);
}


complex EWSM::gAl(const StandardModel::lepton l) const {
    return gAl_SM(l);
}


complex EWSM::gAq(const StandardModel::quark q) const {
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));
    return gAq_SM(q);
}


double EWSM::GammaW() const {
    return GammaW_SM();
}


////////////////////////////////////////////////////////////////////////

double EWSM::Mw_SM() const {
    double Mw;
    if (schemeMw==APPROXIMATEFORMULA) {    
        int NumParMw = 12;
        double params[] = { 
            // Note: Delta alpha depends on the lepton and top-quark masses. 
            // Does Mt depend on threshold scales? 
            SM.getAlsMz(), SM.getMz(), SM.getAle(), SM.getDAle5Mz(),
            SM.getMHl(), 
            SM.getLeptons(SM.NEUTRINO_1).getMass(), 
            SM.getLeptons(SM.NEUTRINO_2).getMass(),
            SM.getLeptons(SM.NEUTRINO_3).getMass(),
            SM.getLeptons(SM.ELECTRON).getMass(),
            SM.getLeptons(SM.MU).getMass(),
            SM.getLeptons(SM.TAU).getMass(),
            SM.getQuarks(SM.TOP).getMass()
        };
        int i = myCache->CacheCheck(Mw_cache, NumParMw, params);
        if (i>=0) {
            return ( Mw_cache[NumParMw][i] );
        } else {
            Mw = myApproximateFormulae->Mw(DeltaAlpha());
            myCache->CacheShift(Mw_cache, NumParMw, params, Mw);
            return Mw;
        }
    } else {
        int NumParMw = 22;
        double params[] = { 
            SM.getAlsMz(), SM.getMz(), SM.getGF(), SM.getAle(), SM.getDAle5Mz(),
            SM.getMHl(), 
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
            SM.getQuarks(SM.TOP).getMass(),
            SM.getQuarks(SM.BOTTOM).getMass(),
            SM.getMut(), SM.getMub(), SM.getMuc(), SM.getMuw()
            //SM.GetLambda(), SM.GetA(), SM.GetRhob(), SM.GetEtab()
        };
        int i = myCache->CacheCheck(Mw_cache, NumParMw, params);
        if (i>=0) {
            return ( Mw_cache[NumParMw][i] );
        } else {
            //std::cout << std::setprecision(12) 
            //          << "TEST: Mw_tree = " << SM.Mw_tree() << std::endl;
            
            double DeltaRho[orders_EW_size], DeltaR_rem[orders_EW_size];
            ComputeDeltaRho(SM.Mw_tree(), DeltaRho);
            ComputeDeltaR_rem(SM.Mw_tree(), DeltaR_rem);
            Mw = resumMw(DeltaRho, DeltaR_rem);
            
            /* Mw from iterations */
            double Mw_org = SM.Mw_tree();
            while (fabs(Mw - Mw_org) > Mw_error) {
                Mw_org = Mw;
                ComputeDeltaRho(Mw, DeltaRho);
                ComputeDeltaR_rem(Mw, DeltaR_rem);
                Mw = resumMw(DeltaRho, DeltaR_rem);
                /* TEST */
                //int prec_def = std::cout.precision();
                //std::cout << std::setprecision(12) << "TEST: Mw_org = " << Mw_org 
                //        << "  Mw_new = " << Mw << std::endl;
                //std::cout.precision(prec_def);
            }
            myCache->CacheShift(Mw_cache, NumParMw, params, Mw);
            return Mw;
        }
    }
}


double EWSM::DeltaR_SM() const {
    double tmp = sqrt(2.0)*SM.getGF()*sW2_SM()*Mw_SM()*Mw_SM()/M_PI/SM.getAle();
    if (schemeMw==NORESUM || schemeMw==APPROXIMATEFORMULA) {
        return (tmp - 1.0);
    } else {
        return (1.0 - 1.0/tmp);
    }
}


double EWSM::cW2_SM() const {
    return ( Mw_SM()*Mw_SM()/SM.getMz()/SM.getMz() );
}


double EWSM::sW2_SM() const {
    return ( 1.0 - cW2_SM() );
}


complex EWSM::rhoZ_l_SM(const StandardModel::lepton l) const {    
    if (schemeRhoZ==APPROXIMATEFORMULA) {
        throw "No approximate formula is available for rhoZ^f";
    } else {
        double Mw = Mw_SM();
        
        /* compute Delta rho */
        double DeltaRho[orders_EW_size];
        ComputeDeltaRho(Mw, DeltaRho);
        
        /* compute delta rho_rem^f */
        complex deltaRho_rem_f[orders_EW_size];
        if (flag_order[EW1]) 
            deltaRho_rem_f[EW1] = myOneLoopEW->deltaRho_rem_l(l,Mw);
        if (flag_order[EW1QCD1]) 
            deltaRho_rem_f[EW1QCD1] = myTwoLoopQCD->deltaRho_rem_l(l,Mw);
        if (flag_order[EW1QCD2]) 
            deltaRho_rem_f[EW1QCD2] = myThreeLoopQCD->deltaRho_rem_l(l,Mw);
        if (flag_order[EW2]) 
            deltaRho_rem_f[EW2] = myTwoLoopEW->deltaRho_rem_l(l,Mw);
        if (flag_order[EW2QCD1]) 
            deltaRho_rem_f[EW2QCD1] = myThreeLoopEW2QCD->deltaRho_rem_l(l,Mw);
        if (flag_order[EW3]) 
            deltaRho_rem_f[EW3] = myThreeLoopEW->deltaRho_rem_l(l,Mw);    
        
        /* compute Delta rbar_rem */
        double DeltaRbar_rem = 0.0;
        if (flag_order[EW1])
            DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(Mw);    
        
        /* Re[rho_Z^f] with or without resummation */
        double deltaRho_rem_f_real[orders_EW_size];
        for (int j=0; j<orders_EW_size; j++)
            deltaRho_rem_f_real[j] = deltaRho_rem_f[j].real();
        double ReRhoZf = resumRhoZ(DeltaRho, deltaRho_rem_f_real, DeltaRbar_rem);
        
        /* Im[rho_Z^f] without resummation */
        double ImRhoZf = 0.0;
        for (int j=0; j<orders_EW_size; j++)
            ImRhoZf += deltaRho_rem_f[j].imag();    
        
        return (complex(ReRhoZf, ImRhoZf, false));    
    }
}


complex EWSM::rhoZ_q_SM(const StandardModel::quark q) const {    
    if (schemeRhoZ==APPROXIMATEFORMULA) {
        throw "No approximate formula is available for rhoZ^f";
    } else {
        double Mw = Mw_SM();
        
        /* compute Delta rho */
        double DeltaRho[orders_EW_size];
        ComputeDeltaRho(Mw, DeltaRho);
        
        /* compute delta rho_rem^f */
        complex deltaRho_rem_f[orders_EW_size];
        if (flag_order[EW1]) 
            deltaRho_rem_f[EW1] = myOneLoopEW->deltaRho_rem_q(q,Mw);
        if (flag_order[EW1QCD1]) 
            deltaRho_rem_f[EW1QCD1] = myTwoLoopQCD->deltaRho_rem_q(q,Mw);
        if (flag_order[EW1QCD2]) 
            deltaRho_rem_f[EW1QCD2] = myThreeLoopQCD->deltaRho_rem_q(q,Mw);
        if (flag_order[EW2]) 
            deltaRho_rem_f[EW2] = myTwoLoopEW->deltaRho_rem_q(q,Mw);
        if (flag_order[EW2QCD1]) 
            deltaRho_rem_f[EW2QCD1] = myThreeLoopEW2QCD->deltaRho_rem_q(q,Mw);
        if (flag_order[EW3]) 
            deltaRho_rem_f[EW3] = myThreeLoopEW->deltaRho_rem_q(q,Mw);    
        
        /* compute Delta rbar_rem */
        double DeltaRbar_rem = 0.0;
        if (flag_order[EW1])
            DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(Mw);    
        
        /* Re[rho_Z^f] with or without resummation */
        double deltaRho_rem_f_real[orders_EW_size];
        for (int j=0; j<orders_EW_size; j++)
            deltaRho_rem_f_real[j] = deltaRho_rem_f[j].real();
        double ReRhoZf = resumRhoZ(DeltaRho, deltaRho_rem_f_real, DeltaRbar_rem);
        
        /* Im[rho_Z^f] without resummation */
        double ImRhoZf = 0.0;
        for (int j=0; j<orders_EW_size; j++)
            ImRhoZf += deltaRho_rem_f[j].imag();    
        
        /* Corrections to the Z-b-bbar vertex */    
        if (q==StandardModel::BOTTOM) {
            throw "Write codes for Zbb in EWSM::rhoZ_f_SM()!";
            
            // Write codes!!    
            
        }
        
        return (complex(ReRhoZf, ImRhoZf, false));    
    }
}


complex EWSM::kappaZ_l_SM(const StandardModel::lepton l) const {
    double Mw = Mw_SM();
    
    /* compute Delta rho */
    double DeltaRho[orders_EW_size];
    ComputeDeltaRho(Mw, DeltaRho);
    
    /* compute delta kappa_rem^f */
    complex deltaKappa_rem_f[orders_EW_size];
    if (flag_order[EW1]) 
        deltaKappa_rem_f[EW1] = myOneLoopEW->deltaKappa_rem_l(l,Mw);
    if (flag_order[EW1QCD1]) 
        deltaKappa_rem_f[EW1QCD1] = myTwoLoopQCD->deltaKappa_rem_l(l,Mw);
    if (flag_order[EW1QCD2]) 
        deltaKappa_rem_f[EW1QCD2] = myThreeLoopQCD->deltaKappa_rem_l(l,Mw);
    if (flag_order[EW2]) 
        deltaKappa_rem_f[EW2] = myTwoLoopEW->deltaKappa_rem_l(l,Mw);
    if (flag_order[EW2QCD1]) 
        deltaKappa_rem_f[EW2QCD1] = myThreeLoopEW2QCD->deltaKappa_rem_l(l,Mw);
    if (flag_order[EW3]) 
        deltaKappa_rem_f[EW3] = myThreeLoopEW->deltaKappa_rem_l(l,Mw);    
    
    /* compute Delta rbar_rem */
    double DeltaRbar_rem = 0.0;
    if (flag_order[EW1])
        DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(Mw);    
    
    double ReKappaZf = 0.0;
    if (schemeKappaZ==APPROXIMATEFORMULA) {
        ReKappaZf = myApproximateFormulae->sin2thetaEff_l(l, DeltaAlpha())/sW2_SM(); 
    } else {    
        /* Re[kappa_Z^f] with or without resummation */
        double deltaKappa_rem_f_real[orders_EW_size];
        for (int j=0; j<orders_EW_size; j++)
            deltaKappa_rem_f_real[j] = deltaKappa_rem_f[j].real();
        ReKappaZf = resumKappaZ(DeltaRho, deltaKappa_rem_f_real, DeltaRbar_rem);
        
        /* O(alpha^2) correction to Re[kappa_Z^f] from the Z-gamma mixing */
        ReKappaZf += 35.0*alphaMz()*alphaMz()/18.0/sW2_SM()
                *(1.0 - 8.0/3.0*ReKappaZf*sW2_SM());
    }
    
    /* Im[kappa_Z^f] without resummation */
    double ImKappaZf = 0.0;
    for (int j=0; j<orders_EW_size; j++)
        ImKappaZf += deltaKappa_rem_f[j].imag();    

    return (complex(ReKappaZf, ImKappaZf, false));       
}


complex EWSM::kappaZ_q_SM(const StandardModel::quark q) const {
    double Mw = Mw_SM();
    
    /* compute Delta rho */
    double DeltaRho[orders_EW_size];
    ComputeDeltaRho(Mw, DeltaRho);
    
    /* compute delta kappa_rem^f */
    complex deltaKappa_rem_f[orders_EW_size];
    if (flag_order[EW1]) 
        deltaKappa_rem_f[EW1] = myOneLoopEW->deltaKappa_rem_q(q,Mw);
    if (flag_order[EW1QCD1]) 
        deltaKappa_rem_f[EW1QCD1] = myTwoLoopQCD->deltaKappa_rem_q(q,Mw);
    if (flag_order[EW1QCD2]) 
        deltaKappa_rem_f[EW1QCD2] = myThreeLoopQCD->deltaKappa_rem_q(q,Mw);
    if (flag_order[EW2]) 
        deltaKappa_rem_f[EW2] = myTwoLoopEW->deltaKappa_rem_q(q,Mw);
    if (flag_order[EW2QCD1]) 
        deltaKappa_rem_f[EW2QCD1] = myThreeLoopEW2QCD->deltaKappa_rem_q(q,Mw);
    if (flag_order[EW3]) 
        deltaKappa_rem_f[EW3] = myThreeLoopEW->deltaKappa_rem_q(q,Mw);    
    
    /* compute Delta rbar_rem */
    double DeltaRbar_rem = 0.0;
    if (flag_order[EW1])
        DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(Mw);    
    
    double ReKappaZf = 0.0;
    if (schemeKappaZ==APPROXIMATEFORMULA) {
        ReKappaZf = myApproximateFormulae->sin2thetaEff_q(q, DeltaAlpha())/sW2_SM(); 
    } else {    
        /* Re[kappa_Z^f] with or without resummation */
        double deltaKappa_rem_f_real[orders_EW_size];
        for (int j=0; j<orders_EW_size; j++)
            deltaKappa_rem_f_real[j] = deltaKappa_rem_f[j].real();
        ReKappaZf = resumKappaZ(DeltaRho, deltaKappa_rem_f_real, DeltaRbar_rem);
        
        /* O(alpha^2) correction to Re[kappa_Z^f] from the Z-gamma mixing */
        ReKappaZf += 35.0*alphaMz()*alphaMz()/18.0/sW2_SM()
                *(1.0 - 8.0/3.0*ReKappaZf*sW2_SM());
    }
    
    /* Im[kappa_Z^f] without resummation */
    double ImKappaZf = 0.0;
    for (int j=0; j<orders_EW_size; j++)
        ImKappaZf += deltaKappa_rem_f[j].imag();    

    /* Corrections to the Z-b-bbar vertex */    
    if (q==StandardModel::BOTTOM) {
        throw "Write codes for Zbb in EWSM::kappaZ_f_SM()!";
        
        // Write codes!!    
        
    }
    
    return (complex(ReKappaZf, ImKappaZf, false));       
}


complex EWSM::gVl_SM(const StandardModel::lepton l) const {
    return ( gAl_SM(l)
             *(1.0 - 4.0*fabs(myCache->Ql(l))*kappaZ_l_SM(l)*sW2_SM()) );
}


complex EWSM::gVq_SM(const StandardModel::quark q) const {
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));
    return ( gAq_SM(q)
             *(1.0 - 4.0*fabs(myCache->Qq(q))*kappaZ_q_SM(q)*sW2_SM()) );
}


complex EWSM::gAl_SM(const StandardModel::lepton l) const {
    double I3f;
    if ( l==StandardModel::NEUTRINO_1 || l==StandardModel::NEUTRINO_2 
            || l==StandardModel::NEUTRINO_3 )
        I3f = 1.0/2.0;
    else if (l==StandardModel::ELECTRON || l==StandardModel::MU
               || l==StandardModel::TAU )
        I3f = - 1.0/2.0;
    else 
        throw "Error in EWSM::gAl_SM()";
    return ( sqrt(rhoZ_l(l))*I3f );
}


complex EWSM::gAq_SM(const StandardModel::quark q) const {
    if (q==StandardModel::TOP) return (complex(0.0, 0.0, false));
    double I3f;
    if ( q==StandardModel::UP || q==StandardModel::CHARM)
        I3f = 1.0/2.0;
    else if (q==StandardModel::DOWN || q==StandardModel::STRANGE
               || q==StandardModel::BOTTOM )
        I3f = - 1.0/2.0;
    else 
        throw "Error in EWSM::gAq_SM()";
    return ( sqrt(rhoZ_q(q))*I3f );
}


double EWSM::rho_GammaW_l_SM(const StandardModel::lepton li, 
                             const StandardModel::lepton lj) const {
    double Mw = Mw_SM();
    double rhoW = 0.0;
    if (flag_order[EW1])
        rhoW = myOneLoopEW->rho_GammaW_l(li,lj,Mw);    
    return rhoW;
}


double EWSM::rho_GammaW_q_SM(const StandardModel::quark qi, 
                             const StandardModel::quark qj) const {
    double Mw = Mw_SM();
    double rhoW = 0.0;
    if (flag_order[EW1])
        rhoW = myOneLoopEW->rho_GammaW_q(qi,qj,Mw);    
    return rhoW;
}


double EWSM::GammaW_l_SM(const StandardModel::lepton li, 
                         const StandardModel::lepton lj) const {
    if ( ((int)li+2)%2 || ((int)lj+3)%2 ) 
        throw "Error in EWSM::GammaW_l_SM()";
    
    double G0 = SM.getGF()*pow(Mw_SM(),3.0)/6.0/sqrt(2.0)/M_PI;    
    complex V(0.0,0.0,false);
    if ( li==StandardModel::NEUTRINO_1 && lj==StandardModel::ELECTRON ||
         li==StandardModel::NEUTRINO_2 && lj==StandardModel::MU ||
         li==StandardModel::NEUTRINO_3 && lj==StandardModel::TAU )        
        V.real() = 1.0;
    return ( V.abs2()*G0*rho_GammaW_l_SM(li,lj) ); 
}
    
    
double EWSM::GammaW_q_SM(const StandardModel::quark qi, 
                         const StandardModel::quark qj) const {
    if ( ((int)qi+2)%2 || ((int)qj+3)%2 ) 
        throw "Error in EWSM::GammaW_q_SM()";
    
    double G0 = SM.getGF()*pow(Mw_SM(),3.0)/6.0/sqrt(2.0)/M_PI;    
    complex V(0.0,0.0,false);
    
    if ( qi==StandardModel::UP && qj==StandardModel::DOWN )
        V = SM.getCKM().V_ud();
    else if ( qi==StandardModel::UP && qj==StandardModel::STRANGE )
        V = SM.getCKM().V_us();
    else if ( qi==StandardModel::UP && qj==StandardModel::BOTTOM )
        V = SM.getCKM().V_ub();
    else if ( qi==StandardModel::CHARM && qj==StandardModel::DOWN )
        V = SM.getCKM().V_cd();
    else if ( qi==StandardModel::CHARM && qj==StandardModel::STRANGE )
        V = SM.getCKM().V_cs();
    else if ( qi==StandardModel::CHARM && qj==StandardModel::BOTTOM )
        V = SM.getCKM().V_cb();
    else if ( qi==StandardModel::TOP || qj==StandardModel::TOP) 
        return (0.0);
    double AlsMw = SM.Als(Mw_SM(), 5.0, SM.getAlsMz(), SM.getMz(), FULLNLO); 
    return ( 3.0*V.abs2()*G0*rho_GammaW_q_SM(qi,qj)*(1.0 + AlsMw/M_PI) );    
}


double EWSM::GammaW_SM() const {
    return ( GammaW_l_SM(SM.NEUTRINO_1, SM.ELECTRON) 
             + GammaW_l_SM(SM.NEUTRINO_2, SM.MU) 
             + GammaW_l_SM(SM.NEUTRINO_3, SM.TAU)             
             + GammaW_q_SM(SM.UP, SM.DOWN) 
             + GammaW_q_SM(SM.UP, SM.STRANGE) 
             + GammaW_q_SM(SM.UP, SM.BOTTOM) 
             + GammaW_q_SM(SM.CHARM, SM.DOWN)
             + GammaW_q_SM(SM.CHARM, SM.STRANGE) 
             + GammaW_q_SM(SM.CHARM, SM.BOTTOM) );
}


////////////////////////////////////////////////////////////////////////     
//
//double EWSM::Mw_NP_fromSTU(const double S, const double T, const double U) const {
//    double CW2minusSW2 = cW2_SM() - sW2_SM();
//    return ( Mw_SM()*Mw_SM()
//             *(1.0 - SM.getAle()*S/2.0/CW2minusSW2 
//               + cW2_SM()*SM.getAle()*T/CW2minusSW2 
//               + SM.getAle()*U/4.0/sW2_SM()) );
//}
// 
//
//complex EWSM::rhoZ_l_NP_fromSTU(const StandardModel::lepton l, const double T) const {
//    return ( (1.0 + SM.getAle()*T)*rhoZ_l_SM(l) );
//}
//
//
//complex EWSM::rhoZ_q_NP_fromSTU(const StandardModel::quark q, const double T) const {
//    return ( (1.0 + SM.getAle()*T)*rhoZ_q_SM(q) );
//}
//
//
//complex EWSM::gZl_over_gAl_NP_fromSTU(const StandardModel::lepton l,
//                                      const double S, const double T) const {
//    double CW2minusSW2 = cW2_SM() - sW2_SM();
//    return ( 1.0 - 4.0*fabs(myCache->Ql(l))
//                   *( 1.0 + SM.getAle()/4.0/sW2_SM()/CW2minusSW2
//                            *(S - 4.0*cW2_SM()*sW2_SM()*T) )
//                   *kappaZ_l_SM(l)*sW2_SM() );    
//}
//
//
//complex EWSM::gZq_over_gAq_NP_fromSTU(const StandardModel::quark q,
//                                      const double S, const double T) const {
//    double CW2minusSW2 = cW2_SM() - sW2_SM();
//    return ( 1.0 - 4.0*fabs(myCache->Qq(q))
//                   *( 1.0 + SM.getAle()/4.0/sW2_SM()/CW2minusSW2
//                            *(S - 4.0*cW2_SM()*sW2_SM()*T) )
//                   *kappaZ_q_SM(q)*sW2_SM() );    
//}
//
//
//double EWSM::GammaW_NP_fromSTU(const double S, const double T, const double U) const {
//    double CW2minusSW2 = cW2_SM() - sW2_SM();
//    return ( GammaW_SM()
//             *(1.0 - SM.getAle()*S/2.0/CW2minusSW2 
//               + cW2_SM()*SM.getAle()*T/CW2minusSW2 
//               + SM.getAle()*U/4.0/sW2_SM()) );
//}
//    
//
////////////////////////////////////////////////////////////////////////

void EWSM::ComputeDeltaRho(const double Mw_i,
                           double DeltaRho[orders_EW_size]) const {
    if (flag_order[EW1]) 
        DeltaRho[EW1] = myOneLoopEW->DeltaRho(Mw_i);
    if (flag_order[EW1QCD1]) 
        DeltaRho[EW1QCD1] = myTwoLoopQCD->DeltaRho(Mw_i);
    if (flag_order[EW1QCD2]) 
        DeltaRho[EW1QCD2] = myThreeLoopQCD->DeltaRho(Mw_i);
    if (flag_order[EW2]) 
        DeltaRho[EW2] = myTwoLoopEW->DeltaRho(Mw_i);
    if (flag_order[EW2QCD1]) 
        DeltaRho[EW2QCD1] = myThreeLoopEW2QCD->DeltaRho(Mw_i);
    if (flag_order[EW3]) 
        DeltaRho[EW3] = myThreeLoopEW->DeltaRho(Mw_i);  
}    


void EWSM::ComputeDeltaR_rem(const double Mw_i,
                             double DeltaR_rem[orders_EW_size]) const {    
    if (flag_order[EW1]) 
        DeltaR_rem[EW1] = myOneLoopEW->DeltaR_rem(Mw_i);
    if (flag_order[EW1QCD1]) 
        DeltaR_rem[EW1QCD1] = myTwoLoopQCD->DeltaR_rem(Mw_i);
    if (flag_order[EW1QCD2]) 
        DeltaR_rem[EW1QCD2] = myThreeLoopQCD->DeltaR_rem(Mw_i);
    if (flag_order[EW2]) 
        DeltaR_rem[EW2] = myTwoLoopEW->DeltaR_rem(Mw_i);
    if (flag_order[EW2QCD1]) 
        DeltaR_rem[EW2QCD1] = myThreeLoopEW2QCD->DeltaR_rem(Mw_i);
    if (flag_order[EW3]) 
        DeltaR_rem[EW3] = myThreeLoopEW->DeltaR_rem(Mw_i);    
}   


double EWSM::resumMw(const double DeltaRho[orders_EW_size],
                     const double DeltaR_rem[orders_EW_size]) const {
    if ( (schemeMw==APPROXIMATEFORMULA) 
            || (DeltaR_rem[EW1QCD2]!=0.0) 
            || (DeltaR_rem[EW2QCD1]!=0.0) 
            || (DeltaR_rem[EW3]!=0.0) )
        throw "Error in EWSM::resumMw()";

    double sW2 = sW2_SM();
    double cW2 = cW2_SM();
    
    double f_AlphaToGF, DeltaRho_sum = 0.0, DeltaRho_G;
    if (schemeMw==NORESUM) {
        for (int j=0; j<orders_EW_size; j++)
            DeltaRho_sum += DeltaRho[(orders_EW)j];
    } else {
        // conversion: alpha(0) --> G_F
        f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)*sW2*cW2/M_PI/SM.getAle();
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
    switch (schemeMw) {
        case NORESUM:
            for (int j=0; j<orders_EW_size; j++)
                DeltaR_rem_sum += DeltaR_rem[(orders_EW)j];
            
            // R = 1 + Delta r
            R = 1.0 + DeltaAlphaL5q() - cW2/sW2*DeltaRho_sum
                + DeltaR_rem_sum;
            break;
        case OMSI:
            // R = 1/(1 - Delta r)
            R = 1.0/(1.0 + cW2/sW2*DeltaRho_sum)
                /(1.0 - DeltaAlphaL5q()
                  - DeltaR_rem[EW1] - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2]);
            break;
        case INTERMEDIATE:
            // R = 1/(1 - Delta r)
            R = 1.0/( (1.0 + cW2/sW2*DeltaRho_sum)
                      *(1.0 - DeltaAlphaL5q() - DeltaR_rem[EW1]) 
                        - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2] );
            break;        
        case OMSII:
            // R = 1/(1 - Delta r)
            R = 1.0/( (1.0 + cW2/sW2*DeltaRho_sum)*(1.0 - DeltaAlphaL5q())
                      - (1.0 + cW2/sW2*DeltaRho_G)*DeltaR_rem[EW1]
                      - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2] );
            break;
        default:
            throw "Error in EWSM::resumMw()";            
    }   

    double tmp = 4.0*M_PI*SM.getAle()/sqrt(2.0)/SM.getGF()/SM.getMz()/SM.getMz();
    if (tmp*R > 1.0) throw "Negative (1-tmp*R) in EWSM::resumMw()";
    
    return (SM.getMz()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp*R)));
}


double EWSM::resumRhoZ(const double DeltaRho[orders_EW_size],
                       const double deltaRho_rem[orders_EW_size],
                       const double DeltaRbar_rem) const {
    if ( (schemeRhoZ==APPROXIMATEFORMULA) 
            || (deltaRho_rem[EW1QCD2]!=0.0) 
            || (deltaRho_rem[EW2QCD1]!=0.0) 
            || (deltaRho_rem[EW3]!=0.0) )
        throw "Error in EWSM::resumRhoZ()";  
    
    double f_AlphaToGF, DeltaRho_sum = 0.0, DeltaRho_G, deltaRho_rem_sum=0.0;
    double DeltaRbar_rem_G, deltaRho_rem_G, deltaRho_rem_G2;
    if (schemeRhoZ==NORESUM) {
        for (int j=0; j<EWSM::orders_EW_size; j++) { 
            DeltaRho_sum += DeltaRho[j];
            deltaRho_rem_sum += deltaRho_rem[j];
        }
    } else {
        // conversion: alpha(0) --> G_F
        f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)
                      *sW2_SM()*cW2_SM()/M_PI/SM.getAle();
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
    }
    
    /* Real parts */
    double rhoZ;
    switch (schemeRhoZ) {
        case NORESUM: 
            rhoZ = 1.0 + DeltaRho_sum + deltaRho_rem_sum;
            break;
        case OMSI:
            rhoZ = (1.0 + deltaRho_rem_G + deltaRho_rem_G2)
                   /(1.0 - DeltaRho_sum*(1.0 - DeltaRbar_rem_G));
            break;
        case INTERMEDIATE:
            rhoZ = (1.0 + deltaRho_rem_G)
                   /(1.0 - DeltaRho_sum*(1.0 - DeltaRbar_rem_G))
                   + deltaRho_rem_G2;            
            break;        
        case OMSII:
            rhoZ = 1.0 + DeltaRho_sum + pow(DeltaRho_G, 2.0) 
                   - DeltaRho_G*DeltaRbar_rem_G
                   + deltaRho_rem_G*(1.0 + DeltaRho_G) + deltaRho_rem_G2;  
            break;
        default:
            throw "Error in EWSM::resumRhoZ()";
    }
    return rhoZ;
}


double EWSM::resumKappaZ(const double DeltaRho[orders_EW_size],
                         const double deltaKappa_rem[orders_EW_size],
                         const double DeltaRbar_rem) const {
    if ( (schemeKappaZ==APPROXIMATEFORMULA)
            || (deltaKappa_rem[EW2QCD1]!=0.0)
            || (deltaKappa_rem[EW3]!=0.0) )
        throw "Error in EWSM::resumKappaZ()";     
    
    double f_AlphaToGF, DeltaRho_sum = 0.0, DeltaRho_G, deltaKappa_rem_sum=0.0;
    double DeltaRbar_rem_G, deltaKappa_rem_G, deltaKappa_rem_G2;
    if (schemeKappaZ==NORESUM) {
        for (int j=0; j<EWSM::orders_EW_size; j++) { 
            DeltaRho_sum += DeltaRho[j];
            deltaKappa_rem_sum += deltaKappa_rem[j];
        }
    } else {
        // conversion: alpha(0) --> G_F
        f_AlphaToGF = sqrt(2.0)*SM.getGF()*pow(SM.getMz(),2.0)
                      *sW2_SM()*cW2_SM()/M_PI/SM.getAle();
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
    }    
    
    /* Real parts */
    double kappaZ;
    switch (schemeKappaZ) {
        case NORESUM: 
            kappaZ = 1.0 + cW2_SM()/sW2_SM()*DeltaRho_sum + deltaKappa_rem_sum;
            break;
        case OMSI:
            kappaZ = (1.0 + deltaKappa_rem_G + deltaKappa_rem_G2)
                     *(1.0 + cW2_SM()/sW2_SM()*DeltaRho_sum*(1.0 - DeltaRbar_rem_G));
            break;
        case INTERMEDIATE:
            kappaZ = (1.0 + deltaKappa_rem_G)
                     *(1.0 + cW2_SM()/sW2_SM()*DeltaRho_sum*(1.0 - DeltaRbar_rem_G))
                     + deltaKappa_rem_G2;
            break;        
        case OMSII:
            kappaZ = 1.0 + cW2_SM()/sW2_SM()*DeltaRho_sum
                     - cW2_SM()/sW2_SM()*DeltaRho_G*DeltaRbar_rem_G
                     + deltaKappa_rem_G*(1.0 + cW2_SM()/sW2_SM()*DeltaRho_G)
                     + deltaKappa_rem_G2;
            break;
        default:
            throw "Error in EWSM::resumKappaZ()";
    }
    return kappaZ;
}



