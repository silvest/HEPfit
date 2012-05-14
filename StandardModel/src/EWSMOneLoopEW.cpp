/* 
 * File:   EWSMOneLoopEW.cpp
 * Author: mishima
 */

#include "EWSMOneLoopEW.h"


EWSMOneLoopEW::EWSMOneLoopEW(const EWSMcache& cache_i) : cache(cache_i) {
}


////////////////////////////////////////////////////////////////////////

double EWSMOneLoopEW::DeltaAlpha_l() const {  
    double Mz = cache.Mz();

    double oneLoop[3];
    oneLoop[0] = - PiGammaGamma_fer(Mz, Mz*Mz, StandardModel::ELECTRON).real() 
                 + PiGammaGamma_fer(Mz, 0.0, StandardModel::ELECTRON).real();
    oneLoop[1] = - PiGammaGamma_fer(Mz, Mz*Mz, StandardModel::MU).real() 
                 + PiGammaGamma_fer(Mz, 0.0, StandardModel::MU).real();
    oneLoop[2] = - PiGammaGamma_fer(Mz, Mz*Mz, StandardModel::TAU).real() 
                 + PiGammaGamma_fer(Mz, 0.0, StandardModel::TAU).real();
    
    return( cache.ale()/4.0/M_PI
            *(oneLoop[0] + oneLoop[1] + oneLoop[2]) );
}


double EWSMOneLoopEW::DeltaAlpha_t() const {   
    double xt = pow(cache.Mz()/cache.Mt(), 2.0);
    double tmp = 1.0 + xt*0.1071;
    tmp *= -4.0/45.0*cache.ale()/M_PI*xt;
    return tmp;
}


double EWSMOneLoopEW::DeltaRho(const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    return ( - cache.ale()/4.0/M_PI/cache.sW2(Mw)*DeltaRhobar(cache.Mz(),Mw) );
}


double EWSMOneLoopEW::DeltaR_rem(const double Mw_i) const {
    double Mz = cache.Mz();
    double Mz2 = Mz*Mz;
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);
    
    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw);
    
    double PiGammaGamma_t_0 = PiGammaGamma_fer(Mz,0.0,StandardModel::TOP).real();
    double PiGammaGamma_l5q_Mz2 = PiGammaGamma_fer(Mz,Mz2).real() 
                                  - PiGammaGamma_fer(Mz,Mz2,StandardModel::TOP).real();

    double DR_rem = - 2.0/3.0*sW2 + sW2*PiGammaGamma_t_0
                    + sW2*PiGammaGamma_l5q_Mz2 + DeltaRhobarW(Mz,Mw)
                    + (4.0 - 25.0/4.0*cW2 + 3.0/4.0*cW2*cW2 + 9.0*cW2/4.0/sW2)
                      *log_cW2
                    + 11.0/2.0 - 5.0/8.0*cW2*(1.0 + cW2);
    DR_rem *= cache.ale()/4.0/M_PI/sW2;
    return DR_rem;    
}


double EWSMOneLoopEW::DeltaRbar_rem(const double Mw_i) const {
    double Mz = cache.Mz();
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);

    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw);    
    
    double DRbar_rem = - 2.0/3.0*sW2 + DeltaRhobarW(Mz,Mw)
                       + (4.0 - 25.0/4.0*cW2 + 3.0/4.0*cW2*cW2 + 9.0*cW2/4.0/sW2)
                         *log_cW2
                       + 11.0/2.0 - 5.0/8.0*cW2*(1.0 + cW2);
    DRbar_rem *= cache.ale()/4.0/M_PI/sW2;
    return DRbar_rem;     
}


complex EWSMOneLoopEW::deltaRho_rem_tmp(const complex uf, 
                                        const double Mw_i) const {
    double Mz = cache.Mz();  
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);

    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw); 

    complex dRho_rem(0.0,0.0,false);
    dRho_rem = - ( SigmaPrime_ZZ_bos_Mz2(Mz,Mw).real() 
                   + SigmaPrime_ZZ_fer_Mz2(Mz,Mw).real() )/cW2
               - DeltaRhobarW(Mz,Mw) + 2.0*uf
               - (1.0/6.0/cW2 - 1.0/3.0 + 3.0/4.0*cW2*(1.0 + cW2) + 9.0*cW2/4.0/sW2)
                 *log_cW2
               - 11.0/2.0 + 5.0/8.0*cW2*(1.0 + cW2);
    dRho_rem *= cache.ale()/4.0/M_PI/sW2;
    return dRho_rem;  
}


template<typename T> 
complex EWSMOneLoopEW::deltaRho_rem_f(const T f, const double Mw_i) const {
    cache.checkSMfermion(f, "EWSMOneLoopEW::deltaRho_rem_f");
    if(f==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    
    double Mz = cache.Mz(); 
    double Mw = cache.Mw(Mw_i);
    complex uf = ( 3.0*cache.vf(f,Mw)*cache.vf(f,Mw) + cache.af(f)*cache.af(f) )
                 /4.0/cache.cW2(Mw)*FZ(Mz*Mz,Mw) + FW(Mz*Mz,f,Mw);
    return ( deltaRho_rem_tmp(uf,Mw) );   
}


complex EWSMOneLoopEW::deltaKappa_rem_tmp(const double deltaf, const complex uf,
                                          const double Mw_i) const {
    double Mz = cache.Mz(); 
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);

    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw); 

    complex dKappa_rem(0.0,0.0,false);
    dKappa_rem = - ( PiZgamma_bos(Mz,Mz*Mz,Mw) + PiZgamma_fer(Mz,Mz*Mz,Mw) )
                 + deltaf*deltaf/4.0/cW2*FZ(Mz*Mz,Mw) - uf
                 + (1.0/12.0/cW2 + 4.0/3.0)*log_cW2;
    dKappa_rem *= cache.ale()/4.0/M_PI/sW2;
    return dKappa_rem;    
}


template<typename T> 
complex EWSMOneLoopEW::deltaKappa_rem_f(const T f, const double Mw_i) const {
    cache.checkSMfermion(f, "EWSMOneLoopEW::deltaKappa_rem_f");
    if(f==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    
    double Mz = cache.Mz(); 
    double Mw = cache.Mw(Mw_i);
    complex uf = ( 3.0*cache.vf(f,Mw)*cache.vf(f,Mw) + cache.af(f)*cache.af(f) )
                 /4.0/cache.cW2(Mw)*FZ(Mz*Mz,Mw) + FW(Mz*Mz,f,Mw);
    return ( deltaKappa_rem_tmp(cache.deltaf(f, Mw), uf, Mw) );       
}


double EWSMOneLoopEW::rho_GammaW_tmp(const double Qi, const double Qj, 
                                     const double Mw_i) const { 
    double QiQj = Qi*Qj;
    double Mz = cache.Mz();
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);

    /* Logarithm and one-loop functions */
    double log_cW2 = cache.log_cW2(Mw); 
    complex B0_Mw_Mw2_Mz_Mw = cache.B0_Mw_Mw2_Mz_Mw(Mw);
    complex C0_Mw2_Mw_0_Mz = cache.C0_Mw2_Mw_0_Mz(Mw);
    
    double V1 = FZa_0(Mw*Mw,Mw).real() - 3.0/2.0;
    double V2 = - 2.0*(2.0 + cW2)*Mz*Mz*C0_Mw2_Mw_0_Mz.real()
                - (1.0/12.0/cW2/cW2 + 5.0/3.0/cW2 + 1.0)*B0_Mw_Mw2_Mz_Mw.real()
                + (1.0/12.0/cW2/cW2 + 1.0/cW2 + 1.0)*log_cW2 
                + 1.0/12.0/cW2/cW2 + 13.0/12.0/cW2 + 59.0/18.0;
    
    double deltafij_W, deltafij_QED;
    deltafij_W = - DeltaRhobarW(Mw,Mw) 
                 - SigmaPrime_WW_bos_Mw2(Mw,Mw).real() 
                 - SigmaPrime_WW_fer_Mw2(Mw,Mw).real()
                 + 5.0/8.0*cW2*(1.0 + cW2) - 11.0/2.0 - 9.0*cW2/4.0/sW2*log_cW2
                 + (-1.0 + 1.0/2.0/cW2 + 2.0*sW2*sW2/cW2*QiQj)*(V1 + 3.0/2.0)
                 + 2.0*cW2*(V2 + 3.0/2.0);
    deltafij_W *= cache.ale()/4.0/M_PI/sW2;

    deltafij_QED = 85.0/18.0 - M_PI*M_PI/3.0 + 3.0/4.0*QiQj;
    deltafij_QED *= cache.ale()/M_PI;

    return ( 1.0 + deltafij_W + deltafij_QED );    
}


double EWSMOneLoopEW::rho_GammaW_l(const StandardModel::lepton li, 
                                   const StandardModel::lepton lj, 
                                   const double Mw_i) const {
    if ( ((int)li+(int)lj+3)%2 ) 
        throw "Error in EWSMOneLoopEW::rho_GammaW_l()";
    double Mw = cache.Mw(Mw_i);
    return ( rho_GammaW_tmp(cache.Qf(li), cache.Qf(lj), Mw) );
}


double EWSMOneLoopEW::rho_GammaW_q(const StandardModel::quark qi, 
                                   const StandardModel::quark qj, 
                                   const double Mw_i) const {
    if ( ((int)qi+(int)qj+3)%2 ) 
        throw "Error in EWSMOneLoopEW::rho_GammaW_q()";
    double Mw = cache.Mw(Mw_i);
    return ( rho_GammaW_tmp(cache.Qf(qi), cache.Qf(qj), Mw) );
}


//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW::SigmaWW_bos(const double mu, const double s,
                                   const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();    
    //double Mz2 = Mz*Mz;
    double mh = cache.mh();
    double mh2 = mh*mh;
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);
    double cW4 = cW2*cW2;
    double RW = pow(Mw, 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;
    double rw = pow(mh/Mw, 2.0);
    
    /* Loop functions */
    double A0_Mw, A0_Mz, A0_mh;
    complex B0_s_Mz_Mw, B0_s_0_Mw, B0_s_mh_Mw;
    complex B0p_s_Mz_Mw, B0p_s_mh_Mw; /* for s==0.0 */
    if (mu==Mz && s==0.0) {
        A0_Mw = cache.A0_Mz_Mw(Mw);
        A0_Mz = cache.A0_Mz_Mz();
        A0_mh = cache.A0_Mz_mh();
        B0_s_Mz_Mw = cache.B0_Mz_0_Mz_Mw(Mw);
        B0_s_0_Mw = cache.B0_Mz_0_0_Mw(Mw);
        B0_s_mh_Mw = cache.B0_Mz_0_mh_Mw(Mw);
        B0p_s_Mz_Mw = cache.B0p_Mz_0_Mz_Mw(Mw);
        B0p_s_mh_Mw = cache.B0p_Mz_0_mh_Mw(Mw);         
    } else if (mu==Mz && s==Mw2) {
        A0_Mw = cache.A0_Mz_Mw(Mw);
        A0_Mz = cache.A0_Mz_Mz();
        A0_mh = cache.A0_Mz_mh();
        B0_s_Mz_Mw = cache.B0_Mz_Mw2_Mz_Mw(Mw);
        B0_s_0_Mw = cache.B0_Mz_Mw2_0_Mw(Mw);
        B0_s_mh_Mw = cache.B0_Mz_Mw2_mh_Mw(Mw);
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        A0_Mw = myPV->A0(mu, Mw);
        A0_Mz = myPV->A0(mu, Mz);
        A0_mh = myPV->A0(mu, mh);
        B0_s_Mz_Mw = myPV->B0(mu, s, Mz, Mw);
        B0_s_0_Mw = myPV->B0(mu, s, 0.0, Mw);
        B0_s_mh_Mw = myPV->B0(mu, s, mh, Mw);
        B0p_s_Mz_Mw = myPV->B0p(mu, s, Mz, Mw);
        B0p_s_mh_Mw = myPV->B0p(mu, s, mh, Mw);
        delete myPV;
    }
    
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        Sigma = Mw2*( 2.0/3.0*(1.0/cW2 - 4.0 - 4.0*cW2 + cW4)*B0_s_Mz_Mw
                      + (1.0/12.0/cW4 + 2.0/3.0/cW2 - 3.0/2.0 
                          + 2.0/3.0*cW2 + 1.0/12.0*cW4)*Mw2*B0p_s_Mz_Mw
                      - 17.0*sW2/6.0*B0_s_0_Mw + 5.0/12.0*sW2
                      - 1.0/12.0*(- 10.0 + 2.0*rw)*B0_s_mh_Mw
                      + 1.0/12.0*pow(1.0 - rw, 2.0)*Mw2*B0p_s_mh_Mw
                      - 1.0/12.0*(24.0 - 2.0*cW2 + cW4)*A0_Mw/Mw2
                      - 1.0/6.0*A0_mh/Mw2
                      - 1.0/12.0*(1.0 + 14.0*cW2 + 9.0*cW4)*A0_Mz/Mw2
                      - 1.0/6.0*(1.0/cW2 + 22.0 + cW2 + cW4 + rw) );     
    } else {
        Sigma = Mw2*( ( (1.0/12.0/cW4 + 2.0/3.0/cW2 - 3.0/2.0 + 2.0/3.0*cW2 
                         + 1.0/12.0*cW4)*RW 
                       + 2.0/3.0*(1.0/cW2 - 4.0 - 4.0*cW2 + cW4)
                       - (3.0/2.0 + 8.0/3.0*cW2 + 3.0/2.0*cW4)/RW 
                       + 2.0/3.0*cW2*(1.0 + cW2)/RW2 + 1.0/12.0*cW4/RW3 )*B0_s_Mz_Mw
                     - sW2/6.0*(- 5.0*RW + 17.0 + 17.0/RW - 5.0/RW2)*B0_s_0_Mw
                     - 1.0/12.0*(- pow(1.0 - rw, 2.0)*RW - 10.0 + 2.0*rw - 1.0/RW)
                       *B0_s_mh_Mw
                     - 1.0/12.0*( (1.0/cW2 - 2.0 + cW2 - cW4 + rw)*RW 
                                   + 24.0 - 2.0*cW2 + cW4
                                   + (- 10.0 + cW2 + cW4)/RW - cW4/RW2 )*A0_Mw/Mw2
                     - 1.0/12.0*( - (1.0/cW2 + 9.0 - 9.0*cW2 - cW4)*RW 
                                     + 1.0 + 14.0*cW2 + 9.0*cW4
                                     + cW2/RW*(1.0 - 9.0*cW2) - cW4/RW2 )*A0_Mz/Mw2 
                     + 1.0/12.0*((mh2 - Mw2)/s - 2.0)*A0_mh/Mw2 
                     - 1.0/6.0*(1.0/cW2 + 22.0 + cW2 + cW4 + rw) 
                     + 1.0/9.0*( (6.0 + 3.0*cW2 + 7.0/2.0*cW4)/RW
                                  - (1.0 + 3.0/2.0*cW2 + 5.0/2.0*cW4)/RW2 
                                  + cW4/2.0/RW3) );    
    }
    return Sigma;
}


complex EWSMOneLoopEW::SigmaWW_fer(const double mu, const double s,
                                   const double Mw_i) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = cache.mf((StandardModel::lepton) i);
        mq[i] = cache.mf((StandardModel::quark) i);        
    }
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();        

    /* Loop functions */
    complex B1_s_ml_mlprime[3], B1_s_mq_mqprime[3]; 
    complex B1_s_mlprime_ml[3], B1_s_mqprime_mq[3];    
    complex Bf_s_mlprime_ml[3], Bf_s_mqprime_mq[3];
    if (mu==Mz && s==0.0) {
        for (int gen=0; gen<3; gen++) {
            B1_s_ml_mlprime[gen] = cache.B1_Mz_0_ml_mlprime(gen);
            B1_s_mq_mqprime[gen] = cache.B1_Mz_0_mq_mqprime(gen);
            B1_s_mlprime_ml[gen] = cache.B1_Mz_0_mlprime_ml(gen);
            B1_s_mqprime_mq[gen] = cache.B1_Mz_0_mqprime_mq(gen);
            Bf_s_mlprime_ml[gen] = cache.Bf_Mz_0_mlprime_ml(gen);
            Bf_s_mqprime_mq[gen] = cache.Bf_Mz_0_mqprime_mq(gen);            
        }
    } else if (mu==Mz && s==Mw2) {
        for (int gen=0; gen<3; gen++) {
            B1_s_ml_mlprime[gen] = cache.B1_Mz_Mw2_ml_mlprime(gen,Mw);
            B1_s_mq_mqprime[gen] = cache.B1_Mz_Mw2_mq_mqprime(gen,Mw);   
            B1_s_mlprime_ml[gen] = cache.B1_Mz_Mw2_mlprime_ml(gen,Mw);
            B1_s_mqprime_mq[gen] = cache.B1_Mz_Mw2_mqprime_mq(gen,Mw);           
            Bf_s_mlprime_ml[gen] = cache.Bf_Mz_Mw2_mlprime_ml(gen,Mw);
            Bf_s_mqprime_mq[gen] = cache.Bf_Mz_Mw2_mqprime_mq(gen,Mw);           
        }
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        for (int gen=0; gen<3; gen++) {
            B1_s_ml_mlprime[gen] = myPV->B1(mu,s,ml[2*gen],ml[2*gen+1]);
            B1_s_mq_mqprime[gen] = myPV->B1(mu,s,mq[2*gen],mq[2*gen+1]);
            B1_s_mlprime_ml[gen] = myPV->B1(mu,s,ml[2*gen+1],ml[2*gen]);
            B1_s_mqprime_mq[gen] = myPV->B1(mu,s,mq[2*gen+1],mq[2*gen]);
            Bf_s_mlprime_ml[gen] = myPV->Bf(mu,s,ml[2*gen+1],ml[2*gen]);
            Bf_s_mqprime_mq[gen] = myPV->Bf(mu,s,mq[2*gen+1],mq[2*gen]);            
        }
        delete myPV;
    }
    
    complex Sigma(0.0,0.0,false);
    double ml2, mlprime2, mq2, mqprime2;
    for (int gen=0; gen<3; gen++) {
        ml2 = ml[2*gen]*ml[2*gen];
        mlprime2 = ml[2*gen+1]*ml[2*gen+1];
        if(s!=0.0) Sigma += - s*Bf_s_mlprime_ml[gen];
        Sigma += mlprime2*B1_s_ml_mlprime[gen] + ml2*B1_s_mlprime_ml[gen];
        //
        mq2 = mq[2*gen]*mq[2*gen];
        mqprime2 = mq[2*gen+1]*mq[2*gen+1];
        if(s!=0.0) Sigma += 3.0*( - s*Bf_s_mqprime_mq[gen] );
        Sigma += 3.0*( mqprime2*B1_s_mq_mqprime[gen] + mq2*B1_s_mqprime_mq[gen] );
    }
    return Sigma;
}


complex EWSMOneLoopEW::SigmaZZ_bos(const double mu, const double s,
                                   const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;
    double mh = cache.mh();
    double mh2 = mh*mh;
    double cW2 = cache.cW2(Mw);
    double cW4 = cW2*cW2;
    double RW = pow(Mw, 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;
    double rw = pow(mh/Mw, 2.0);
    
    /* Loop functions */
    double A0_Mw, A0_Mz, A0_mh;
    complex B0_s_Mw_Mw, B0_s_mh_Mz;
    if (mu==Mz && s==Mz2) {
        A0_Mw = cache.A0_Mz_Mw(Mw);
        A0_Mz = cache.A0_Mz_Mz();
        A0_mh = cache.A0_Mz_mh();
        B0_s_Mw_Mw = cache.B0_Mz_Mz2_Mw_Mw(Mw);
        B0_s_mh_Mz = cache.B0_Mz_Mz2_mh_Mz();
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        A0_Mw = myPV->A0(mu, Mw);
        A0_Mz = myPV->A0(mu, Mz);
        A0_mh = myPV->A0(mu, mh);
        B0_s_Mw_Mw = myPV->B0(mu, s, Mw, Mw);
        B0_s_mh_Mz = myPV->B0(mu, s, mh, Mz);
        delete myPV;
    }        

    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw "Missing codes for EWSMOneLoopEW::SigmaZZ_bos(s=0.0)";
    } else {
        Sigma = Mw2*( - cW4*(4.0 + 17.0/3.0/RW - 4.0/3.0/RW2 - 1.0/12.0/RW3 )
                        *B0_s_Mw_Mw
                      + 1.0/12.0*( (1.0/cW4 - 2.0/cW2*rw + rw*rw)*RW
                                   + 10.0/cW2 - 2.0*rw + 1.0/RW )*B0_s_mh_Mz
                      - cW2*(4.0 - 4.0/3.0/RW - 1.0/6.0/RW2)*A0_Mw/Mz2
                      + 1.0/12.0*((Mz2 - mh2)/s + 1.0)*(A0_Mz - A0_mh)/cW2/Mz2
                      - 1.0/12.0*A0_mh/cW2/Mz2
                      - ( 1.0/6.0/cW2 + 4.0*cW4 + 1.0/6.0*rw
                          - (1.0/18.0 + 4.0/3.0*cW4)/RW 
                          + 1.0/9.0*cW4*(5.0 - 1.0/2.0/RW)/RW2 ) );
    }
    return Sigma;
}


complex EWSMOneLoopEW::SigmaZZ_fer(const double mu, const double s, 
                                   const double Mw_i) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = cache.mf((StandardModel::lepton) i);
        mq[i] = cache.mf((StandardModel::quark) i);        
    }
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_ml_ml[6], Bf_s_mq_mq[6];
    complex B0_s_ml_ml[6], B0_s_mq_mq[6];    
    if (mu==Mz && s==Mz2) {
        for (int i=0; i<6; i++) {
            Bf_s_ml_ml[i] = cache.Bf_Mz_Mz2_mf_mf((StandardModel::lepton) i);
            Bf_s_mq_mq[i] = cache.Bf_Mz_Mz2_mf_mf((StandardModel::quark) i);           
            B0_s_ml_ml[i] = cache.B0_Mz_Mz2_mf_mf((StandardModel::lepton) i);
            B0_s_mq_mq[i] = cache.B0_Mz_Mz2_mf_mf((StandardModel::quark) i);           
        }
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        for (int i=0; i<6; i++) {
            Bf_s_ml_ml[i] = myPV->Bf(mu,s,ml[i],ml[i]);
            Bf_s_mq_mq[i] = myPV->Bf(mu,s,mq[i],mq[i]);            
            B0_s_ml_ml[i] = myPV->B0(mu,s,ml[i],ml[i]);
            B0_s_mq_mq[i] = myPV->B0(mu,s,mq[i],mq[i]);            
        }
        delete myPV;
    }
    
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw "Missing codes for EWSMOneLoopEW::SigmaZZ_fer(s=0.0)";
    } else {
        double ml2, vl2, al2, mq2, vq2, aq2;
        for (int i=0; i<6; i++) {
            ml2 = ml[i]*ml[i];
            vl2 = pow(cache.vf((StandardModel::lepton) i, Mw), 2.0);
            al2 = pow(cache.af((StandardModel::lepton) i), 2.0);            
            if(s!=0.0) Sigma += - (vl2 + al2)*s*Bf_s_ml_ml[i];
            Sigma += - 2.0*al2*ml2*B0_s_ml_ml[i];
            //
            mq2 = mq[i]*mq[i];
            vq2 = pow(cache.vf((StandardModel::quark) i, Mw), 2.0);
            aq2 = pow(cache.af((StandardModel::quark) i), 2.0);
            if(s!=0.0) Sigma += - 3.0*(vq2 + aq2)*s*Bf_s_mq_mq[i];
            Sigma += - 3.0*2.0*aq2*mq2*B0_s_mq_mq[i];
        }
    }   
    return Sigma;
}


complex EWSMOneLoopEW::PiGammaGamma_bos(const double mu, const double s,
                                        const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;
    double RW = pow(Mw, 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;    
    
    /* Loop functions */
    double A0_Mw;
    complex B0_s_Mw_Mw;
    if (mu==Mz && s==Mz2) {
        A0_Mw = cache.A0_Mz_Mw(Mw);
        B0_s_Mw_Mw = cache.B0_Mz_Mz2_Mw_Mw(Mw);
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        A0_Mw = myPV->A0(mu, Mw);
        B0_s_Mw_Mw = myPV->B0(mu, s, Mw, Mw);
        delete myPV; 
    }
    
    complex Pi(0.0,0.0,false);
    if (s==0.0) {
        throw "Missing codes for EWSMOneLoopEW::PiGammaGamma_bos(s=0.0)";
    } else {
        Pi = - RW*( (4.0 + 17.0/3.0/RW - 4.0/3.0/RW2 - 1.0/12.0/RW3)*B0_s_Mw_Mw
                    + (4.0 - 4.0/3.0/RW - 1.0/6.0/RW2)*(A0_Mw/Mw2 + 1.0)
                    - 1.0/18.0/RW2*(1.0/RW - 13.0) );
    }
    return Pi;
}


template<typename T> 
complex EWSMOneLoopEW::PiGammaGamma_fer(const double mu, const double s, 
                                        const T f) const {
    double mf;
    if ( (typeid(f)==typeid(StandardModel::quark)) 
         || (typeid(f)==typeid(StandardModel::lepton)) )
        mf = cache.mf(f);
    else 
        throw "Error in EWSMOneLoopEW::PiGammaGamma_fer()";
    
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_mf_mf;
    if (mu==Mz && s==Mz2) {
        if (mf==0.0) { 
            Bf_s_mf_mf = 0.0; 
        } else {
            Bf_s_mf_mf = cache.Bf_Mz_Mz2_mf_mf(f);
        }
    } else if (mu==Mz && s==0.0) {        
        if (mf==0.0) {
            Bf_s_mf_mf = 0.0; 
        } else {
            Bf_s_mf_mf = cache.Bf_Mz_0_mf_mf(f);
        }
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        if (mf==0.0) {
            Bf_s_mf_mf = 0.0; 
        } else {
            Bf_s_mf_mf = myPV->Bf(mu,s,mf,mf);
        }
        delete myPV;
    }
    
    double Qf = cache.Qf(f);
    return ( - 4.0*Qf*Qf*Bf_s_mf_mf);
}


complex EWSMOneLoopEW::PiGammaGamma_fer(const double mu, const double s) const {
    complex Pi(0.0,0.0,false);
    for (int i=0; i<6; i++) {
        Pi += PiGammaGamma_fer(mu, s, (StandardModel::lepton) i);
        Pi += PiGammaGamma_fer(mu, s, (StandardModel::quark) i);        
    }
    return Pi;
}


complex EWSMOneLoopEW::PiZgamma_bos(const double mu, const double s,
                                    const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double cW2 = cache.cW2(Mw);
    return ( - PiGammaGamma_bos(mu,s,Mw)*cW2);
}


complex EWSMOneLoopEW::PiZgamma_fer(const double mu, const double s,
                                    const double Mw_i) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = cache.mf((StandardModel::lepton) i);
        mq[i] = cache.mf((StandardModel::quark) i); 
    }
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    
    /* Loop functions */
    complex Bf_s_ml_ml[6], Bf_s_mq_mq[6];
    if (mu==Mz && s==Mz2) {
        for (int i=0; i<6; i++) {
            Bf_s_ml_ml[i] = cache.Bf_Mz_Mz2_mf_mf((StandardModel::lepton) i);
            Bf_s_mq_mq[i] = cache.Bf_Mz_Mz2_mf_mf((StandardModel::quark) i);           
        }
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        for (int i=0; i<6; i++) {
            Bf_s_ml_ml[i] = myPV->Bf(mu,s,ml[i],ml[i]);
            Bf_s_mq_mq[i] = myPV->Bf(mu,s,mq[i],mq[i]);            
        }
        delete myPV;
    }

    complex Pi(0.0,0.0,false);
    double Ql, Qq;
    for (int i=0; i<6; i++) {
        Ql = cache.Qf((StandardModel::lepton) i);
        Pi += (fabs(Ql) - 4.0*sW2*Ql*Ql)*Bf_s_ml_ml[i];
        //
        Qq = cache.Qf((StandardModel::quark) i);
        Pi += 3.0*(fabs(Qq) - 4.0*sW2*Qq*Qq)*Bf_s_mq_mq[i];
    }   
    return Pi;
}
 

//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW::SigmaPrime_WW_bos_Mw2(const double mu, 
                                             const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();    
    //double Mz2 = Mz*Mz;
    double mh = cache.mh();
    double mh2 = mh*mh;
    double sW2 = cache.sW2(Mw_i);
    double cW2 = cache.cW2(Mw_i);
    double cW4 = cW2*cW2;
    double rw = mh2/Mw2;
    //double rz = mh2/Mz2;
    
    /* Loop functions */
    double A0_Mw, A0_Mz, A0_mh;
    complex B0_Mw2_Mz_Mw, B0_Mw2_0_Mw, B0_Mw2_mh_Mw;
    complex B0p_Mw2_Mz_Mw, B0p_Mw2_0_Mw, B0p_Mw2_mh_Mw;
    if (mu==Mw) {
        A0_Mw = cache.A0_Mw_Mw(Mw);
        A0_Mz = cache.A0_Mw_Mz(Mw);
        A0_mh = cache.A0_Mw_mh(Mw);
        B0_Mw2_Mz_Mw = cache.B0_Mw_Mw2_Mz_Mw(Mw);
        B0_Mw2_0_Mw = cache.B0_Mw_Mw2_0_Mw(Mw);
        B0_Mw2_mh_Mw = cache.B0_Mw_Mw2_mh_Mw(Mw);
        B0p_Mw2_Mz_Mw = cache.B0p_Mw_Mw2_Mz_Mw(Mw);
        B0p_Mw2_0_Mw = cache.B0p_Mw_Mw2_0_Mw(Mw);
        B0p_Mw2_mh_Mw = cache.B0p_Mw_Mw2_mh_Mw(Mw);
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        A0_Mw = myPV->A0(mu, Mw);
        A0_Mz = myPV->A0(mu, Mz);
        A0_mh = myPV->A0(mu, mh);
        B0_Mw2_Mz_Mw = myPV->B0(mu, Mw2, Mz, Mw);
        B0_Mw2_0_Mw = myPV->B0(mu, Mw2, 0.0, Mw);
        B0_Mw2_mh_Mw = myPV->B0(mu, Mw2, mh, Mw);
        B0p_Mw2_Mz_Mw = myPV->B0p(mu, Mw2, Mz, Mw);
        B0p_Mw2_0_Mw = myPV->B0p(mu, Mw2, 0.0, Mw);
        B0p_Mw2_mh_Mw = myPV->B0p(mu, Mw2, mh, Mw);
        delete myPV;
    }
    
    complex Sigma(0.0,0.0,false);
    Sigma = - (1.0/12.0/cW4 + 2.0/3.0/cW2 + 2.0*cW2)*B0_Mw2_Mz_Mw
            + (1.0/12.0/cW4 + 4.0/3.0/cW2 - 17.0/3.0 - 4.0*cW2)*Mw2*B0p_Mw2_Mz_Mw
            - 2.0*sW2*B0_Mw2_0_Mw - 4.0*sW2*Mw2*B0p_Mw2_0_Mw
            + rw/6.0*(1.0 - rw/2.0)*B0_Mw2_mh_Mw
            + (1.0 - rw/3.0 + rw*rw/12.0)*Mw2*B0p_Mw2_mh_Mw
            + (1.0/cW2 + 8.0 + rw)/12.0*A0_Mw/Mw2
            - (1.0/cW2 + 9.0 - 8.0*cW2 - 12.0*cW4)/12.0*A0_Mz/Mw2
            - (rw - 1.0)/12.0*A0_mh/Mw2 + 4.0/9.0;
    return Sigma;    
}


complex EWSMOneLoopEW::SigmaPrime_WW_fer_Mw2(const double mu,
                                             const double Mw_i) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = cache.mf((StandardModel::lepton) i);
        mq[i] = cache.mf((StandardModel::quark) i); 
    }
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;

    /* Loop functions */
    complex Bf_Mw2_mlprime_ml[3], Bf_Mw2_mqprime_mq[3];
    complex Bfp_Mw2_mlprime_ml[3], Bfp_Mw2_mqprime_mq[3];    
    complex B1p_Mw2_mlprime_ml[3], B1p_Mw2_mqprime_mq[3];  
    complex B1p_Mw2_ml_mlprime[3], B1p_Mw2_mq_mqprime[3];    
    if (mu==Mw) {
        for (int gen=0; gen<3; gen++) {
            Bf_Mw2_mlprime_ml[gen] = cache.Bf_Mw_Mw2_mlprime_ml(gen,Mw);
            Bf_Mw2_mqprime_mq[gen] = cache.Bf_Mw_Mw2_mqprime_mq(gen,Mw);              
            Bfp_Mw2_mlprime_ml[gen] = cache.Bfp_Mw_Mw2_mlprime_ml(gen,Mw);
            Bfp_Mw2_mqprime_mq[gen] = cache.Bfp_Mw_Mw2_mqprime_mq(gen,Mw);            
            B1p_Mw2_ml_mlprime[gen] = cache.B1p_Mw_Mw2_ml_mlprime(gen,Mw);
            B1p_Mw2_mq_mqprime[gen] = cache.B1p_Mw_Mw2_mq_mqprime(gen,Mw);
            B1p_Mw2_mlprime_ml[gen] = cache.B1p_Mw_Mw2_mlprime_ml(gen,Mw);
            B1p_Mw2_mqprime_mq[gen] = cache.B1p_Mw_Mw2_mqprime_mq(gen,Mw);
        }
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        for (int gen=0; gen<3; gen++) {
            Bf_Mw2_mlprime_ml[gen] = myPV->Bf(mu,Mw2,ml[2*gen+1],ml[2*gen]);
            Bf_Mw2_mqprime_mq[gen] = myPV->Bf(mu,Mw2,mq[2*gen+1],mq[2*gen]);            
            Bfp_Mw2_mlprime_ml[gen] = myPV->Bfp(mu,Mw2,ml[2*gen+1],ml[2*gen]);
            Bfp_Mw2_mqprime_mq[gen] = myPV->Bfp(mu,Mw2,mq[2*gen+1],mq[2*gen]);            
            B1p_Mw2_ml_mlprime[gen] = myPV->B1p(mu,Mw2,ml[2*gen],ml[2*gen+1]);
            B1p_Mw2_mq_mqprime[gen] = myPV->B1p(mu,Mw2,mq[2*gen],mq[2*gen+1]);
            B1p_Mw2_mlprime_ml[gen] = myPV->B1p(mu,Mw2,ml[2*gen+1],ml[2*gen]);
            B1p_Mw2_mqprime_mq[gen] = myPV->B1p(mu,Mw2,mq[2*gen+1],mq[2*gen]);
        }        
        delete myPV;
    }

    complex Sigma(0.0,0.0,false);
    double ml2, mlprime2, mq2, mqprime2;
    for (int gen=0; gen<3; gen++) {
        ml2 = ml[2*gen]*ml[2*gen];
        mlprime2 = ml[2*gen+1]*ml[2*gen+1];
        Sigma += - (Bf_Mw2_mlprime_ml[gen] + Mw2*Bfp_Mw2_mlprime_ml[gen]);
        Sigma += mlprime2*B1p_Mw2_ml_mlprime[gen] + ml2*B1p_Mw2_mlprime_ml[gen];
        //
        mq2 = mq[2*gen]*mq[2*gen];
        mqprime2 = mq[2*gen+1]*mq[2*gen+1];
        Sigma += - 3.0*(Bf_Mw2_mqprime_mq[gen] + Mw2*Bfp_Mw2_mqprime_mq[gen]);
        Sigma += 3.0*( mqprime2*B1p_Mw2_mq_mqprime[gen] + mq2*B1p_Mw2_mqprime_mq[gen] );
    }
    return Sigma;    
}


complex EWSMOneLoopEW::SigmaPrime_ZZ_bos_Mz2(const double mu,
                                             const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;
    double mh = cache.mh();
    double mh2 = mh*mh;
    double cW2 = cache.cW2(Mw);
    double cW4 = cW2*cW2;
    double rw = mh2/Mw2;
    double rz = mh2/Mz2;
    
    /* Loop functions */
    double A0_Mw, A0_Mz, A0_mh;
    complex B0_Mz2_Mw_Mw, B0_Mz2_mh_Mz;
    complex B0p_Mz2_Mw_Mw, B0p_Mz2_mh_Mz;
    if (mu==Mz) {
        A0_Mw = cache.A0_Mz_Mw(Mw);
        A0_Mz = cache.A0_Mz_Mz();
        A0_mh = cache.A0_Mz_mh();
        B0_Mz2_Mw_Mw = cache.B0_Mz_Mz2_Mw_Mw(Mw);
        B0_Mz2_mh_Mz = cache.B0_Mz_Mz2_mh_Mz();
        B0p_Mz2_Mw_Mw = cache.B0p_Mz_Mz2_Mw_Mw(Mw);
        B0p_Mz2_mh_Mz = cache.B0p_Mz_Mz2_mh_Mz();
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        A0_Mw = myPV->A0(mu, Mw);
        A0_Mz = myPV->A0(mu, Mz);
        A0_mh = myPV->A0(mu, mh);
        B0_Mz2_Mw_Mw = myPV->B0(mu, Mz2, Mw, Mw);
        B0_Mz2_mh_Mz = myPV->B0(mu, Mz2, mh, Mz);
        B0p_Mz2_Mw_Mw = myPV->B0p(mu, Mz2, Mw, Mw);
        B0p_Mz2_mh_Mz = myPV->B0p(mu, Mz2, mh, Mz);
        delete myPV;
    }

    complex Sigma(0.0,0.0,false);
    Sigma = (1.0/4.0/cW2 + 8.0/3.0 - 17.0/3.0*cW2)*B0_Mz2_Mw_Mw
            + (1.0/12.0/cW2 + 4.0/3.0 - 17.0/3.0*cW2 - 4.0*cW4)*Mz2*B0p_Mz2_Mw_Mw
            + rw/6.0*(1.0 - rz/2.0)*B0_Mz2_mh_Mz 
            + (1.0 - rz/3.0 + rz*rz/12.0)*Mz2/cW2*B0p_Mz2_mh_Mz
            + (1.0 + 4.0*cW2)/3.0/cW2*A0_Mw/Mz2
            - (1.0 - rz)/12.0/cW2*(A0_Mz - A0_mh)/Mz2
            + 2.0/9.0/cW2 - 10.0/9.0 + 4.0/3.0*cW2;
    Sigma *= cW2;
    return Sigma;        
}


complex EWSMOneLoopEW::SigmaPrime_ZZ_fer_Mz2(const double mu, const double Mw_i) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = cache.mf((StandardModel::lepton) i);
        mq[i] = cache.mf((StandardModel::quark) i); 
    }
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_Mz2_ml_ml[6], Bf_Mz2_mq_mq[6];
    complex Bfp_Mz2_ml_ml[6], Bfp_Mz2_mq_mq[6];
    complex B0p_Mz2_ml_ml[6], B0p_Mz2_mq_mq[6]; 
    if (mu==Mz) {
         for (int i=0; i<6; i++) {
             Bf_Mz2_ml_ml[i] = cache.Bf_Mz_Mz2_mf_mf((StandardModel::lepton) i);
             Bf_Mz2_mq_mq[i] = cache.Bf_Mz_Mz2_mf_mf((StandardModel::quark) i);           
             Bfp_Mz2_ml_ml[i] = cache.Bfp_Mz_Mz2_mf_mf((StandardModel::lepton) i);
             Bfp_Mz2_mq_mq[i] = cache.Bfp_Mz_Mz2_mf_mf((StandardModel::quark) i);           
             B0p_Mz2_ml_ml[i] = cache.B0p_Mz_Mz2_mf_mf((StandardModel::lepton) i);
             B0p_Mz2_mq_mq[i] = cache.B0p_Mz_Mz2_mf_mf((StandardModel::quark) i);           
         }        
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        for (int i=0; i<6; i++) {
            Bf_Mz2_ml_ml[i] = myPV->Bf(mu,Mz2,ml[i],ml[i]);
            Bf_Mz2_mq_mq[i] = myPV->Bf(mu,Mz2,mq[i],mq[i]);            
            Bfp_Mz2_ml_ml[i] = myPV->Bfp(mu,Mz2,ml[i],ml[i]);
            Bfp_Mz2_mq_mq[i] = myPV->Bfp(mu,Mz2,mq[i],mq[i]);            
            B0p_Mz2_ml_ml[i] = myPV->B0p(mu,Mz2,ml[i],ml[i]);
            B0p_Mz2_mq_mq[i] = myPV->B0p(mu,Mz2,mq[i],mq[i]);            
        }        
        delete myPV;
   }
    
    complex Sigma(0.0,0.0,false);
    double ml2, vl2, al2, mq2, vq2, aq2;
    for (int i=0; i<6; i++) {
        ml2 = ml[i]*ml[i];
        vl2 = pow(cache.vf((StandardModel::lepton) i, Mw), 2.0);
        al2 = pow(cache.af((StandardModel::lepton) i), 2.0);            
        Sigma += - (vl2 + al2)*(Bf_Mz2_ml_ml[i] + Mz2*Bfp_Mz2_ml_ml[i])
                 - 2.0*al2*ml2*B0p_Mz2_ml_ml[i];
        //
        mq2 = mq[i]*mq[i];
        vq2 = pow(cache.vf((StandardModel::quark) i, Mw), 2.0);
        aq2 = pow(cache.af((StandardModel::quark) i), 2.0);
        Sigma += - 3.0*(vq2 + aq2)*(Bf_Mz2_mq_mq[i] + Mz2*Bfp_Mz2_mq_mq[i])
                 - 6.0*aq2*mq2*B0p_Mz2_mq_mq[i];
    }
    return Sigma;    
}


//////////////////////////////////////////////////////////////////////// 

double EWSMOneLoopEW::DeltaRhobar(const double mu, const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();    
    return ( (SigmaWW_bos(mu,Mw*Mw,Mw).real() + SigmaWW_fer(mu,Mw*Mw,Mw).real() 
              - SigmaZZ_bos(mu,Mz*Mz,Mw).real() - SigmaZZ_fer(mu,Mz*Mz,Mw).real())
             /Mw/Mw );
}


double EWSMOneLoopEW::DeltaRhobarW(const double mu, const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    return ( (SigmaWW_bos(mu,0.0,Mw).real() + SigmaWW_fer(mu,0.0,Mw).real() 
              - SigmaWW_bos(mu,Mw*Mw,Mw).real() - SigmaWW_fer(mu,Mw*Mw,Mw).real())
             /Mw/Mw );
}


//////////////////////////////////////////////////////////////////////// 

double EWSMOneLoopEW::TEST_DeltaRhobar_bos(const double Mw_i) const {
    double Mz = cache.Mz();    
    double mh = cache.mh();
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);
    double cW4 = cW2*cW2;
    double rw = pow(mh/Mw, 2.0);
    double rz = pow(mh/Mz, 2.0);
    
    /* Logarithm and one-loop functions */
    double log_cW2 = cache.log_cW2(Mw); 
    
    /* B0 functions for mu=Mw */
    complex B0_Mw_Mz2_Mw_Mw = cache.B0_Mz_Mz2_Mw_Mw(Mw) + log_cW2;
    complex B0_Mw_Mz2_mh_Mz = cache.B0_Mz_Mz2_mh_Mz() + log_cW2;
    complex B0_Mw_Mw2_Mz_Mw = cache.B0_Mz_Mw2_Mz_Mw(Mw) + log_cW2;
    complex B0_Mw_Mw2_mh_Mw = cache.B0_Mz_Mw2_mh_Mw(Mw) + log_cW2;
    
    double DRhobar;    
    DRhobar = - (1.0/12.0/cW2 + 4.0/3.0 - 17.0/3.0*cW2 - 4.0*cW4)
                 *(B0_Mw_Mz2_Mw_Mw.real() - 1.0/cW2*B0_Mw_Mw2_Mz_Mw.real())
              + (1.0 - 1.0/3.0*rw + 1.0/12.0*rw*rw)*B0_Mw_Mw2_mh_Mw.real()
              - (1.0 - 1.0/3.0*rz + 1.0/12.0*rz*rz)/cW2*B0_Mw_Mz2_mh_Mz.real()
              + 1.0/12.0*sW2*rw*rw*(log(rw) - 1.0)
              - (1.0/12.0/cW4 + 1.0/2.0/cW2 - 2.0 + 1.0/12.0*rw)*log_cW2
              - 1.0/12.0/cW4 - 19.0/36.0/cW2 - 133.0/18.0 + 8.0*cW2;
    return DRhobar;
} 


double EWSMOneLoopEW::TEST_DeltaRhobarW_bos(const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);
    double cW4 = cW2*cW2;
    double mh = cache.mh();
    double rw = pow(mh/Mw, 2.0);

    /* Logarithm and one-loop functions */
    double log_cW2 = cache.log_cW2(Mw); 
    
    /* B0 functions for mu=Mw */
    complex B0_Mw2_Mz_Mw = cache.B0_Mz_Mw2_Mz_Mw(Mw) + log_cW2;
    complex B0_Mw2_mh_Mw = cache.B0_Mz_Mw2_mh_Mw(Mw) + log_cW2;
    
    double DRhobarW;    
    DRhobarW = - (1.0/12.0/cW4 + 4.0/3.0/cW2 - 17.0/3.0 - 4.0*cW2)*B0_Mw2_Mz_Mw.real()
               - (1.0 - 1.0/3.0*rw + 1.0/12.0*rw*rw)*B0_Mw2_mh_Mw.real()
               + (3.0/4.0/(1.0-rw) + 1.0/4.0 - 1.0/12.0*rw)*rw*log(rw)
               + (1.0/12.0/cW4 + 17.0/12.0/cW2 - 3.0/sW2 + 1.0/4.0)*log_cW2
               + 1.0/12.0/cW4 + 11.0/8.0/cW2 + 139.0/36.0 - 177.0/24.0*cW2 
               + 5.0/8.0*cW4 - 1.0/12.0*rw*(7.0/2.0 - rw);     
    return DRhobarW;
} 

   
//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW::FZa_0(const double s, const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();   
    double Rz = Mz*Mz/s;

    /* Logarithm and three-point one-loop functions */
    double log_Rz;
    complex C0_s_0_Mz_0;
    if (s==Mz*Mz) {
        log_Rz = 0.0;
        C0_s_0_Mz_0 = cache.C0_Mz2_0_Mz_0();
    } else if (s==Mw*Mw) {
        log_Rz = - cache.log_cW2(Mw);
        C0_s_0_Mz_0 = cache.C0_Mw2_0_Mz_0(Mw);
    } else {
        log_Rz = log(Rz);
        PVfunctions* myPV;
        myPV = new PVfunctions();
        C0_s_0_Mz_0 = myPV->C0(s,0.0,Mz,0.0);
        delete myPV;
    }
    
    complex FZa(0.0,0.0,false);
    FZa = 2.0*pow((Rz + 1.0), 2.0)*s*C0_s_0_Mz_0
          - (2.0*Rz + 3.0)*(log_Rz + M_PI*complex::i()) - 2.0*Rz - 7.0/2.0;
    return FZa;    
}


complex EWSMOneLoopEW::FWa_0(const double s, const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;

    /* Logarithm and three-point one-loop functions */
    double log_Rw;
    complex C0_s_0_Mw_0;
    if (s==Mz*Mz) {    
        log_Rw = cache.log_cW2(Mw);
        C0_s_0_Mw_0 = cache.C0_Mz2_0_Mw_0(Mw); 
    } else {
        log_Rw = log(Rw);
        PVfunctions* myPV;
        myPV = new PVfunctions();
        C0_s_0_Mw_0 = myPV->C0(s,0.0,Mw,0.0);
        delete myPV;    
    }
    
    complex FWa(0.0,0.0,false);
    FWa = 2.0*pow((Rw + 1.0), 2.0)*s*C0_s_0_Mw_0
          - (2.0*Rw + 3.0)*(log_Rw + M_PI*complex::i()) - 2.0*Rw - 7.0/2.0;
    return FWa;     
}


complex EWSMOneLoopEW::FbarWa_0(const double s) const {
    return ( complex(0.0,0.0,false) );
}


complex EWSMOneLoopEW::FWn_0(const double s, const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;

    /* Two- and three-point one-loop functions */
    complex B0_Mw_s_Mw_Mw;
    complex C0_s_Mw_0_Mw;
    if (s==Mz*Mz) {
        B0_Mw_s_Mw_Mw = cache.B0_Mw_Mz2_Mw_Mw(Mw);
        C0_s_Mw_0_Mw = cache.C0_Mz2_Mw_0_Mw(Mw);
    } else {    
        PVfunctions* myPV;
        myPV = new PVfunctions();
        B0_Mw_s_Mw_Mw = myPV->B0(Mw,s,Mw,Mw); 
        C0_s_Mw_0_Mw = myPV->C0(s,Mw,0.0,Mw);
        delete myPV;  
    }
 
    complex FWn(0.0,0.0,false);    
    FWn = - 2.0*(Rw + 2.0)*Mw*Mw*C0_s_Mw_0_Mw 
          - (2.0*Rw + 7.0/3.0 - 3.0/2.0/Rw - 1.0/12.0/Rw/Rw)*B0_Mw_s_Mw_Mw
          + 2.0*Rw + 9.0/2.0 - 11.0/18.0/Rw + 1.0/18.0/Rw/Rw;  
    return FWn;
}
    

complex EWSMOneLoopEW::FWa_t(const double s, const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;
    double Mt = cache.Mt();
    double wt = Mt*Mt/Mw/Mw; 

    /* Logarithm and two- and three-point one-loop functions */
    double log_wt = - 2.0*cache.logMZtoMTOP() - cache.log_cW2(Mw);
    double log_Rw;
    complex B0_Mw_s_Mt_Mt;
    complex C0_s_Mt_Mw_Mt, C0_s_0_Mw_0;
    if (s==Mz*Mz) {
        log_Rw = cache.log_cW2(Mw);
        B0_Mw_s_Mt_Mt = cache.B0_Mw_Mz2_Mt_Mt(Mw);
        C0_s_Mt_Mw_Mt = cache.C0_Mz2_Mt_Mw_Mt(Mw);
        C0_s_0_Mw_0 = cache.C0_Mz2_0_Mw_0(Mw);
    } else {
        log_Rw = log(Rw);
        PVfunctions* myPV;
        myPV = new PVfunctions();
        B0_Mw_s_Mt_Mt = myPV->B0(Mw,s,Mt,Mt); 
        C0_s_Mt_Mw_Mt = myPV->C0(s,Mt,Mw,Mt);
        C0_s_0_Mw_0 = myPV->C0(s,0.0,Mw,0.0);
        delete myPV;      
    }
    
    complex FWa(0.0,0.0,false);    
    FWa = 2.0*(Rw + 1.0)*(Rw + 1.0)*s*(C0_s_Mt_Mw_Mt - C0_s_0_Mw_0)
          + (2.0*Rw + 3.0)*(- B0_Mw_s_Mt_Mt + log_Rw + M_PI*complex::i() + 2.0)
          - wt*( (3.0*Rw + 2.0 - wt - wt*wt*Rw)*Mw*Mw*C0_s_Mt_Mw_Mt
                 + (Rw + 1.0/2.0 + wt*Rw)*(1.0 - B0_Mw_s_Mt_Mt)
                 - (2.0*Rw + 1.0/2.0 - 2.0/(wt - 1.0) 
                    + 3.0/2.0/(wt - 1.0)/(wt - 1.0) + wt*Rw)*log_wt
                 + 3.0/2.0/(wt - 1.0) + 3.0/4.0 );
    return FWa;
}


complex EWSMOneLoopEW::FbarWa_t(const double s, const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;
    double Mt = cache.Mt();
    double wt = Mt*Mt/Mw/Mw; 

    /* Logarithm and two- and three-point one-loop functions */
    double log_wt = - 2.0*cache.logMZtoMTOP() - cache.log_cW2(Mw);
    complex B0_Mw_s_Mt_Mt;
    complex C0_s_Mt_Mw_Mt;    
    if (s==Mz*Mz) {
        B0_Mw_s_Mt_Mt = cache.B0_Mw_Mz2_Mt_Mt(Mw);
        C0_s_Mt_Mw_Mt = cache.C0_Mz2_Mt_Mw_Mt(Mw);
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        B0_Mw_s_Mt_Mt = myPV->B0(Mw,s,Mt,Mt);     
        C0_s_Mt_Mw_Mt = myPV->C0(s,Mt,Mw,Mt);    
        delete myPV;   
    }
    
    complex FbarWa(0.0,0.0,false);        
    FbarWa = - wt*( (Rw + 2.0 - wt*(2.0 - wt)*Rw)*Mw*Mw*C0_s_Mt_Mw_Mt
                    - (1.0/2.0 - Rw + wt*Rw)*(- B0_Mw_s_Mt_Mt + 1.0) 
                    + wt*Rw*log_wt );
    return FbarWa;
}


complex EWSMOneLoopEW::FWn_t(const double s, const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;
    double Mt = cache.Mt();
    double wt = Mt*Mt/Mw/Mw; 

    /* Logarithm and two- and three-point one-loop functions */
    double log_wt = - 2.0*cache.logMZtoMTOP() - cache.log_cW2(Mw);  
    complex B0_Mw_s_Mw_Mw;
    complex C0_s_Mw_Mt_Mw, C0_s_Mw_0_Mw;
    if (s==Mz*Mz) {
        B0_Mw_s_Mw_Mw = cache.B0_Mw_Mz2_Mw_Mw(Mw);
        C0_s_Mw_Mt_Mw = cache.C0_Mz2_Mw_Mt_Mw(Mw);
        C0_s_Mw_0_Mw = cache.C0_Mz2_Mw_0_Mw(Mw);  
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        B0_Mw_s_Mw_Mw = myPV->B0(Mw,s,Mw,Mw); 
        C0_s_Mw_Mt_Mw = myPV->C0(s,Mw,Mt,Mw);
        C0_s_Mw_0_Mw = myPV->C0(s,Mw,0.0,Mw);    
        delete myPV;     
    }
    
    complex FWn(0.0,0.0,false);        
    FWn = - 2.0*(Rw + 2.0)*Mw*Mw*(C0_s_Mw_Mt_Mw - C0_s_Mw_0_Mw)
          + wt*( (3.0*Rw + 5.0/2.0 - 2.0/Rw - wt*(2.0 - 1.0/2.0/Rw) 
                  + wt*wt*(1.0/2.0 - Rw))*Mw*Mw*C0_s_Mw_Mt_Mw
                  + (Rw + 1.0 - 1.0/4.0/Rw - wt*(1.0/2.0 - Rw))*(B0_Mw_s_Mw_Mw - 1.0)
                  + (2.0*Rw + 1.0/2.0 - 2.0/(wt - 1.0) 
                     + 3.0/2.0/(wt - 1.0)/(wt - 1.0) - wt*(1.0/2.0 - Rw))*log_wt
                  - 3.0/2.0/(wt- 1.0) + 1.0/4.0 
                  - 1.0/2.0/Rw );
    return FWn;
}


complex EWSMOneLoopEW::FZ(const double s, const double Mw_i) const {
    return ( FZa_0(s, Mw_i) );
}


template<typename T> 
complex EWSMOneLoopEW::FW(const double s, const T f, const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double cW2 = cache.cW2(Mw);

    if (typeid(f)==typeid(StandardModel::quark)) {
        StandardModel::quark qprime;
        switch(f) {
            case StandardModel::UP:
                qprime = StandardModel::DOWN;
                break;
            case StandardModel::DOWN:
                qprime = StandardModel::UP;
                break;
            case StandardModel::CHARM:
                qprime = StandardModel::STRANGE;
                break;
            case StandardModel::STRANGE:
                qprime = StandardModel::CHARM;
                break;
            case StandardModel::TOP:
                throw "TOP is not allowed in EWSMOneLoopEW::FW(s,q)";
            case StandardModel::BOTTOM:
                qprime = StandardModel::TOP;
                break;
            default:
                throw "Error in EWSMOneLoopEW::FW()";  
        }         
        complex FW(0.0,0.0,false);
        FW = cW2*FWn_0(s,Mw) - cache.sigmaf(qprime, Mw)/2.0*FWa_0(s, Mw) 
             - FbarWa_0(s)/2.0;
        if (f==StandardModel::BOTTOM) {
            FW += cW2*FWn_t(s,Mw) - cache.sigmaf(qprime, Mw)/2.0*FWa_t(s, Mw) 
                  - FbarWa_t(s, Mw)/2.0;
        }
        return FW;
    } else if (typeid(f)==typeid(StandardModel::lepton)) {
        StandardModel::lepton lprime;
        switch(f) {
            case StandardModel::NEUTRINO_1:
                lprime = StandardModel::ELECTRON;
                break;
            case StandardModel::ELECTRON:
                lprime = StandardModel::NEUTRINO_1;
                break;
            case StandardModel::NEUTRINO_2:
                lprime = StandardModel::MU;
                break;
            case StandardModel::MU:
                lprime = StandardModel::NEUTRINO_2;
                break;
            case StandardModel::NEUTRINO_3:
                lprime = StandardModel::TAU;
                break;
            case StandardModel::TAU:
                lprime = StandardModel::NEUTRINO_3;
                break;
            default:
                throw "Error in EWSMOneLoopEW::FW()";  
        }
        return ( cW2*FWn_0(s, Mw) - cache.sigmaf(lprime, Mw)/2.0*FWa_0(s, Mw) 
                - FbarWa_0(s)/2.0 );
    } else {
        throw "Error in EWSMOneLoopEW::FW()";
    }
}
   

//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW::TEST_FWn(const double s, const double mf, 
                                const double Mw_i) const {
    double Mw = cache.Mw(Mw_i);
    double Rw = Mw*Mw/s;
    double wf = mf*mf/Mw/Mw; 

    /* Logarithm and two- and three-point one-loop functions */
    double log_wf = log(wf);  
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double A0_Mw = myPV->A0(Mw, Mw);
    double A0_mf = myPV->A0(Mw, mf);
    complex B0_Mw_s_Mw_Mw = myPV->B0(Mw,s,Mw,Mw); 
    complex B0_Mw_0_mf_Mw = myPV->B0(Mw,0.0,mf,Mw); 
    complex C0_s_Mw_mf_Mw = myPV->C0(s,Mw,mf,Mw);
    complex C0_s_Mw_0_Mw = myPV->C0(s,Mw,0.0,Mw);    
    delete myPV;     
    
    complex FWn(0.0,0.0,false);        
    /* Eq.(5.586) in Bardin and Passarino's book */
    FWn = ((2.0 + wf)*(1.0 - wf)*(1.0 - wf)*Rw + 4.0 - 5.0/2.0*wf + 2.0*wf*wf
            - wf*wf*wf/2.0 + wf*(2.0 - wf/2.0)/Rw)*Mw*Mw*C0_s_Mw_mf_Mw
           - (- (2.0 + wf)*(1.0 - wf)*Rw - 3.0 + 3.0/2.0*wf - wf*wf/2.0)
             *(B0_Mw_s_Mw_Mw - B0_Mw_0_mf_Mw)
           - (2.0/3.0 - wf/2.0 + (3.0/2.0 - wf/4.0)/Rw + 1.0/12.0/Rw/Rw)*B0_Mw_s_Mw_Mw
           - A0_mf/Mw/Mw - (2.0/3.0 + 1.0/6.0/Rw)*A0_Mw/Mw/Mw
           - 2.0/3.0 - wf/2.0 + (4.0/9.0 + wf/4.0)/Rw - 1.0/18.0/Rw/Rw;
    
    /* Eq.(5.500) in Bardin and Passarino's book */    
    double wW;
    if (mf==0.0) {
        wW = 3.0/2.0;
    } else {
        wW = 5.0/4.0*wf + 3.0*(1.0 - wf/2.0)*wf*wf/(wf - 1.0)/(wf - 1.0)*log_wf
             - 3.0/2.0/(wf - 1.0);        
    }

    return ( - FWn + wW );
}





