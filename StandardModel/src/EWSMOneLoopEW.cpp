/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSMOneLoopEW.h"


EWSMOneLoopEW::EWSMOneLoopEW(const EWSMcache& cache_i) 
: cache(cache_i) 
{
}


////////////////////////////////////////////////////////////////////////

double EWSMOneLoopEW::DeltaAlpha_l(const double s) const 
{  
    double Mz = cache.Mz();

    double oneLoop[3];
    oneLoop[0] = - PiGammaGamma_fer_l(Mz, s, StandardModel::ELECTRON).real() 
                 + PiGammaGamma_fer_l(Mz, 0.0, StandardModel::ELECTRON).real();
    oneLoop[1] = - PiGammaGamma_fer_l(Mz, s, StandardModel::MU).real() 
                 + PiGammaGamma_fer_l(Mz, 0.0, StandardModel::MU).real();
    oneLoop[2] = - PiGammaGamma_fer_l(Mz, s, StandardModel::TAU).real() 
                 + PiGammaGamma_fer_l(Mz, 0.0, StandardModel::TAU).real();
    
    return( cache.ale()/4.0/M_PI
            *(oneLoop[0] + oneLoop[1] + oneLoop[2]) );
}


double EWSMOneLoopEW::DeltaAlpha_5q(const double s) const
{
    double Mz = cache.Mz();

    double oneLoop[5];
    /* Qf and Nc are included in PiGammaGamma_fer_q(). */
    oneLoop[0] = - PiGammaGamma_fer_q(Mz, s, StandardModel::UP).real()
                 + PiGammaGamma_fer_q(Mz, 0.0, StandardModel::UP).real();
    oneLoop[1] = - PiGammaGamma_fer_q(Mz, s, StandardModel::DOWN).real()
                 + PiGammaGamma_fer_q(Mz, 0.0, StandardModel::DOWN).real();
    oneLoop[2] = - PiGammaGamma_fer_q(Mz, s, StandardModel::CHARM).real()
                 + PiGammaGamma_fer_q(Mz, 0.0, StandardModel::CHARM).real();
    oneLoop[3] = - PiGammaGamma_fer_q(Mz, s, StandardModel::STRANGE).real()
                 + PiGammaGamma_fer_q(Mz, 0.0, StandardModel::STRANGE).real();
    oneLoop[4] = - PiGammaGamma_fer_q(Mz, s, StandardModel::BOTTOM).real()
                 + PiGammaGamma_fer_q(Mz, 0.0, StandardModel::BOTTOM).real();

    return( cache.ale()/4.0/M_PI
            *(oneLoop[0] + oneLoop[1] + oneLoop[2] + oneLoop[3] + oneLoop[4]) );
}


double EWSMOneLoopEW::DeltaAlpha_t(const double s) const 
{   
    double xt = s/cache.Mt()/cache.Mt();
    double tmp = 1.0 + xt*0.1071;
    tmp *= -4.0/45.0*cache.ale()/M_PI*xt;
    return tmp;
}


double EWSMOneLoopEW::DeltaRho(const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    return ( - cache.ale()/4.0/M_PI/cache.sW2(Mw)*DeltaRhobar(cache.Mz(),Mw) );
}


double EWSMOneLoopEW::DeltaR_rem(const double Mw_i) const 
{
    double Mz = cache.Mz();
    double Mz2 = Mz*Mz;
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);
    
    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw);
    
    double PiGammaGamma_t_0 = PiGammaGamma_fer_q(Mz,0.0,StandardModel::TOP).real();
    double PiGammaGamma_l5q_Mz2 = PiGammaGamma_fer(Mz,Mz2).real() 
                                  - PiGammaGamma_fer_q(Mz,Mz2,StandardModel::TOP).real();

    double DR_rem = - 2.0/3.0*sW2 + sW2*PiGammaGamma_t_0
                    + sW2*PiGammaGamma_l5q_Mz2 + DeltaRhobarW(Mz,Mw)
                    + (4.0 - 25.0/4.0*cW2 + 3.0/4.0*cW2*cW2 + 9.0*cW2/4.0/sW2)
                      *log_cW2
                    + 11.0/2.0 - 5.0/8.0*cW2*(1.0 + cW2);
    DR_rem *= cache.ale()/4.0/M_PI/sW2;
    return DR_rem;    
}


double EWSMOneLoopEW::DeltaRbar_rem(const double Mw_i) const 
{
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
                                        const double Mw_i) const 
{
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


complex EWSMOneLoopEW::deltaRho_rem_l(const StandardModel::lepton l, const double Mw_i) const 
{
    double Mz = cache.Mz(); 
    double Mw = cache.Mw(Mw_i);
    complex uf = ( 3.0*cache.vl(l,Mw)*cache.vl(l,Mw) + cache.al(l)*cache.al(l) )
                 /4.0/cache.cW2(Mw)*FZ(Mz*Mz,Mw) + FW_l(Mz*Mz,l,Mw);
    return ( deltaRho_rem_tmp(uf,Mw) );   
}


complex EWSMOneLoopEW::deltaRho_rem_q(const StandardModel::quark q, const double Mw_i) const 
{
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    
    double Mz = cache.Mz(); 
    double Mw = cache.Mw(Mw_i);
    complex uf = ( 3.0*cache.vq(q,Mw)*cache.vq(q,Mw) + cache.aq(q)*cache.aq(q) )
                 /4.0/cache.cW2(Mw)*FZ(Mz*Mz,Mw) + FW_q(Mz*Mz,q,Mw);    
    return ( deltaRho_rem_tmp(uf,Mw) );   
}


complex EWSMOneLoopEW::deltaKappa_rem_tmp(const double deltaf, const complex uf,
                                          const double Mw_i) const 
{
    double Mz = cache.Mz(); 
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);

    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw); 

    complex dKappa_rem(0.0,0.0,false);
    dKappa_rem = ( PiZgamma_bos(Mz,Mz*Mz,Mw) + PiZgamma_fer(Mz,Mz*Mz,Mw) )
                 + deltaf*deltaf/4.0/cW2*FZ(Mz*Mz,Mw) - uf
                 + (1.0/12.0/cW2 + 4.0/3.0)*log_cW2;
    dKappa_rem *= cache.ale()/4.0/M_PI/sW2;
    return dKappa_rem;    
}


complex EWSMOneLoopEW::deltaKappa_rem_l(const StandardModel::lepton l, const double Mw_i) const
{
    double Mz = cache.Mz(); 
    double Mw = cache.Mw(Mw_i);
    complex uf = ( 3.0*cache.vl(l,Mw)*cache.vl(l,Mw) + cache.al(l)*cache.al(l) )
                 /4.0/cache.cW2(Mw)*FZ(Mz*Mz,Mw) + FW_l(Mz*Mz,l,Mw);
    return ( deltaKappa_rem_tmp(cache.deltal(l, Mw), uf, Mw) );       
}


complex EWSMOneLoopEW::deltaKappa_rem_q(const StandardModel::quark q, const double Mw_i) const 
{
    if(q==StandardModel::TOP) return ( complex(0.0,0.0,false) );
    
    double Mz = cache.Mz(); 
    double Mw = cache.Mw(Mw_i);
    complex uf = ( 3.0*cache.vq(q,Mw)*cache.vq(q,Mw) + cache.aq(q)*cache.aq(q) )
                 /4.0/cache.cW2(Mw)*FZ(Mz*Mz,Mw) + FW_q(Mz*Mz,q,Mw);
    return ( deltaKappa_rem_tmp(cache.deltaq(q, Mw), uf, Mw) );       
}


double EWSMOneLoopEW::rho_GammaW_tmp(const double Qi, const double Qj, 
                                     const double Mw_i) const 
{ 
    double QiQj = Qi*Qj;
    double Mz = cache.Mz();
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);

    /* Logarithm and one-loop functions */
    double log_cW2 = cache.log_cW2(Mw); 
    complex B0_Mw2_Mw2_Mz2_Mw2 = cache.B0_Mw2_Mw2_Mz2_Mw2(Mw);
    complex C0_Mw2_Mw2_0_Mz2 = cache.C0_Mw2_Mw2_0_Mz2(Mw);
    
    double V1 = FZa_0(Mw*Mw,Mw).real() - 3.0/2.0;
    double V2 = - 2.0*(2.0 + cW2)*Mz*Mz*C0_Mw2_Mw2_0_Mz2.real()
                - (1.0/12.0/cW2/cW2 + 5.0/3.0/cW2 + 1.0)*B0_Mw2_Mw2_Mz2_Mw2.real()
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
                                   const double Mw_i) const 
{
    if ( ((int)li+(int)lj+3)%2 ) 
        throw std::runtime_error("Error in EWSMOneLoopEW::rho_GammaW_l()"); 
    double Mw = cache.Mw(Mw_i);
    return ( rho_GammaW_tmp(cache.Ql(li), cache.Ql(lj), Mw) );
}


double EWSMOneLoopEW::rho_GammaW_q(const StandardModel::quark qi, 
                                   const StandardModel::quark qj, 
                                   const double Mw_i) const 
{
    if ( ((int)qi+(int)qj+3)%2 ) 
        throw std::runtime_error("Error in EWSMOneLoopEW::rho_GammaW_q()"); 
    double Mw = cache.Mw(Mw_i);
    return ( rho_GammaW_tmp(cache.Qq(qi), cache.Qq(qj), Mw) );
}


//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW::SigmaWW_bos(const double mu, const double s,
                                   const double Mw_i) const 
{
    double mu2 = mu*mu;
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;
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
    double A0_Mw2, A0_Mz2, A0_mh2;
    complex B0_s_Mz2_Mw2, B0_s_0_Mw2, B0_s_mh2_Mw2;
    complex B0p_s_Mz2_Mw2, B0p_s_mh2_Mw2; /* for s==0.0 */
    if (mu==Mz && s==0.0) {
        A0_Mw2 = cache.A0_Mz2_Mw2(Mw);
        A0_Mz2 = cache.A0_Mz2_Mz2();
        A0_mh2 = cache.A0_Mz2_mh2();
        B0_s_Mz2_Mw2 = cache.B0_Mz2_0_Mz2_Mw2(Mw);
        B0_s_0_Mw2 = cache.B0_Mz2_0_0_Mw2(Mw);
        B0_s_mh2_Mw2 = cache.B0_Mz2_0_mh2_Mw2(Mw);
        B0p_s_Mz2_Mw2 = cache.B0p_Mz2_0_Mz2_Mw2(Mw);
        B0p_s_mh2_Mw2 = cache.B0p_Mz2_0_mh2_Mw2(Mw);
    } else if (mu==Mz && s==Mw2) {
        A0_Mw2 = cache.A0_Mz2_Mw2(Mw);
        A0_Mz2 = cache.A0_Mz2_Mz2();
        A0_mh2 = cache.A0_Mz2_mh2();
        B0_s_Mz2_Mw2 = cache.B0_Mz2_Mw2_Mz2_Mw2(Mw);
        B0_s_0_Mw2 = cache.B0_Mz2_Mw2_0_Mw2(Mw);
        B0_s_mh2_Mw2 = cache.B0_Mz2_Mw2_mh2_Mw2(Mw);
    } else {
        A0_Mw2 = cache.getPV().A0(mu2, Mw2);
        A0_Mz2 = cache.getPV().A0(mu2, Mz2);
        A0_mh2 = cache.getPV().A0(mu2, mh2);
        B0_s_Mz2_Mw2 = cache.getPV().B0(mu2, s, Mz2, Mw2);
        B0_s_0_Mw2 = cache.getPV().B0(mu2, s, 0.0, Mw2);
        B0_s_mh2_Mw2 = cache.getPV().B0(mu2, s, mh2, Mw2);
        B0p_s_Mz2_Mw2 = cache.getPV().B0p(mu2, s, Mz2, Mw2);
        B0p_s_mh2_Mw2 = cache.getPV().B0p(mu2, s, mh2, Mw2);
    }
    
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        Sigma = Mw2*( 2.0/3.0*(1.0/cW2 - 4.0 - 4.0*cW2 + cW4)*B0_s_Mz2_Mw2
                      + (1.0/12.0/cW4 + 2.0/3.0/cW2 - 3.0/2.0 
                          + 2.0/3.0*cW2 + 1.0/12.0*cW4)*Mw2*B0p_s_Mz2_Mw2
                      - 17.0*sW2/6.0*B0_s_0_Mw2 + 5.0/12.0*sW2
                      - 1.0/12.0*(- 10.0 + 2.0*rw)*B0_s_mh2_Mw2
                      + 1.0/12.0*pow(1.0 - rw, 2.0)*Mw2*B0p_s_mh2_Mw2
                      - 1.0/12.0*(24.0 - 2.0*cW2 + cW4)*A0_Mw2/Mw2
                      - 1.0/6.0*A0_mh2/Mw2
                      - 1.0/12.0*(1.0 + 14.0*cW2 + 9.0*cW4)*A0_Mz2/Mw2
                      - 1.0/6.0*(1.0/cW2 + 22.0 + cW2 + cW4 + rw) );     
    } else {
        Sigma = Mw2*( ( (1.0/12.0/cW4 + 2.0/3.0/cW2 - 3.0/2.0 + 2.0/3.0*cW2 
                         + 1.0/12.0*cW4)*RW 
                       + 2.0/3.0*(1.0/cW2 - 4.0 - 4.0*cW2 + cW4)
                       - (3.0/2.0 + 8.0/3.0*cW2 + 3.0/2.0*cW4)/RW 
                       + 2.0/3.0*cW2*(1.0 + cW2)/RW2 + 1.0/12.0*cW4/RW3 )*B0_s_Mz2_Mw2
                     - sW2/6.0*(- 5.0*RW + 17.0 + 17.0/RW - 5.0/RW2)*B0_s_0_Mw2
                     - 1.0/12.0*(- pow(1.0 - rw, 2.0)*RW - 10.0 + 2.0*rw - 1.0/RW)
                       *B0_s_mh2_Mw2
                     - 1.0/12.0*( (1.0/cW2 - 2.0 + cW2 - cW4 + rw)*RW 
                                   + 24.0 - 2.0*cW2 + cW4
                                   + (- 10.0 + cW2 + cW4)/RW - cW4/RW2 )*A0_Mw2/Mw2
                     - 1.0/12.0*( - (1.0/cW2 + 9.0 - 9.0*cW2 - cW4)*RW 
                                     + 1.0 + 14.0*cW2 + 9.0*cW4
                                     + cW2/RW*(1.0 - 9.0*cW2) - cW4/RW2 )*A0_Mz2/Mw2
                     + 1.0/12.0*((mh2 - Mw2)/s - 2.0)*A0_mh2/Mw2
                     - 1.0/6.0*(1.0/cW2 + 22.0 + cW2 + cW4 + rw) 
                     + 1.0/9.0*( (6.0 + 3.0*cW2 + 7.0/2.0*cW4)/RW
                                  - (1.0 + 3.0/2.0*cW2 + 5.0/2.0*cW4)/RW2 
                                  + cW4/2.0/RW3) );    
    }
    return Sigma;
}


complex EWSMOneLoopEW::SigmaWW_fer(const double mu, const double s,
                                   const double Mw_i) const 
{
    double ml2[6], mq2[6];
    for (int i=0; i<6; i++) { 
        ml2[i] = cache.ml2((StandardModel::lepton) i);
        mq2[i] = cache.mq2((StandardModel::quark) i, mu);
    }
    double mu2 = mu*mu;
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();        

    /* Loop functions */
    complex B1_s_ml2_mlprime2[3], B1_s_mq2_mqprime2[3];
    complex B1_s_mlprime2_ml2[3], B1_s_mqprime2_mq2[3];
    complex Bf_s_mlprime2_ml2[3], Bf_s_mqprime2_mq2[3];
    if (mu==Mz && s==0.0) {
        for (int gen=0; gen<3; gen++) {
            B1_s_ml2_mlprime2[gen] = cache.B1_Mz2_0_ml2_mlprime2(gen);
            B1_s_mq2_mqprime2[gen] = cache.B1_Mz2_0_mq2_mqprime2(gen);
            B1_s_mlprime2_ml2[gen] = cache.B1_Mz2_0_mlprime2_ml2(gen);
            B1_s_mqprime2_mq2[gen] = cache.B1_Mz2_0_mqprime2_mq2(gen);
            Bf_s_mlprime2_ml2[gen] = cache.Bf_Mz2_0_mlprime2_ml2(gen);
            Bf_s_mqprime2_mq2[gen] = cache.Bf_Mz2_0_mqprime2_mq2(gen);
        }
    } else if (mu==Mz && s==Mw2) {
        for (int gen=0; gen<3; gen++) {
            B1_s_ml2_mlprime2[gen] = cache.B1_Mz2_Mw2_ml2_mlprime2(gen,Mw);
            B1_s_mq2_mqprime2[gen] = cache.B1_Mz2_Mw2_mq2_mqprime2(gen,Mw);
            B1_s_mlprime2_ml2[gen] = cache.B1_Mz2_Mw2_mlprime2_ml2(gen,Mw);
            B1_s_mqprime2_mq2[gen] = cache.B1_Mz2_Mw2_mqprime2_mq2(gen,Mw);
            Bf_s_mlprime2_ml2[gen] = cache.Bf_Mz2_Mw2_mlprime2_ml2(gen,Mw);
            Bf_s_mqprime2_mq2[gen] = cache.Bf_Mz2_Mw2_mqprime2_mq2(gen,Mw);
        }
    } else {
        for (int gen=0; gen<3; gen++) {
            B1_s_ml2_mlprime2[gen] = cache.getPV().B1(mu2,s,ml2[2*gen],ml2[2*gen+1]);
            B1_s_mq2_mqprime2[gen] = cache.getPV().B1(mu2,s,mq2[2*gen],mq2[2*gen+1]);
            B1_s_mlprime2_ml2[gen] = cache.getPV().B1(mu2,s,ml2[2*gen+1],ml2[2*gen]);
            B1_s_mqprime2_mq2[gen] = cache.getPV().B1(mu2,s,mq2[2*gen+1],mq2[2*gen]);
            Bf_s_mlprime2_ml2[gen] = cache.getPV().Bf(mu2,s,ml2[2*gen+1],ml2[2*gen]);
            Bf_s_mqprime2_mq2[gen] = cache.getPV().Bf(mu2,s,mq2[2*gen+1],mq2[2*gen]);
        }
    }
    
    complex Sigma(0.0,0.0,false);
    double mf2, mfprime2;
    for (int gen=0; gen<3; gen++) {
        mf2 = ml2[2*gen];
        mfprime2 = ml2[2*gen+1];
        if(s!=0.0) Sigma += - s*Bf_s_mlprime2_ml2[gen];
        Sigma += mfprime2*B1_s_ml2_mlprime2[gen] + mf2*B1_s_mlprime2_ml2[gen];
        //
        mf2 = mq2[2*gen];
        mfprime2 = mq2[2*gen+1];
        if(s!=0.0) Sigma += 3.0*( - s*Bf_s_mqprime2_mq2[gen] );
        Sigma += 3.0*( mfprime2*B1_s_mq2_mqprime2[gen] + mf2*B1_s_mqprime2_mq2[gen] );
    }
    return Sigma;
}


complex EWSMOneLoopEW::SigmaZZ_bos(const double mu, const double s,
                                   const double Mw_i) const 
{
    double mu2 = mu*mu;
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
    double A0_Mw2, A0_Mz2, A0_mh2;
    complex B0_s_Mw2_Mw2, B0_s_mh2_Mz2;
    if (mu==Mz && s==Mz2) {
        A0_Mw2 = cache.A0_Mz2_Mw2(Mw);
        A0_Mz2 = cache.A0_Mz2_Mz2();
        A0_mh2 = cache.A0_Mz2_mh2();
        B0_s_Mw2_Mw2 = cache.B0_Mz2_Mz2_Mw2_Mw2(Mw);
        B0_s_mh2_Mz2 = cache.B0_Mz2_Mz2_mh2_Mz2();
    } else {
        A0_Mw2 = cache.getPV().A0(mu2, Mw2);
        A0_Mz2 = cache.getPV().A0(mu2, Mz2);
        A0_mh2 = cache.getPV().A0(mu2, mh2);
        B0_s_Mw2_Mw2 = cache.getPV().B0(mu2, s, Mw2, Mw2);
        B0_s_mh2_Mz2 = cache.getPV().B0(mu2, s, mh2, Mz2);
    }        

    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW::SigmaZZ_bos(s=0.0)"); 
    } else {
        Sigma = Mw2*( - cW4*(4.0 + 17.0/3.0/RW - 4.0/3.0/RW2 - 1.0/12.0/RW3 )
                        *B0_s_Mw2_Mw2
                      + 1.0/12.0*( (1.0/cW4 - 2.0/cW2*rw + rw*rw)*RW
                                   + 10.0/cW2 - 2.0*rw + 1.0/RW )*B0_s_mh2_Mz2
                      - cW2*(4.0 - 4.0/3.0/RW - 1.0/6.0/RW2)*A0_Mw2/Mz2
                      + 1.0/12.0*((Mz2 - mh2)/s + 1.0)*(A0_Mz2 - A0_mh2)/cW2/Mz2
                      - 1.0/12.0*A0_mh2/cW2/Mz2
                      - ( 1.0/6.0/cW2 + 4.0*cW4 + 1.0/6.0*rw
                          - (1.0/18.0 + 4.0/3.0*cW4)/RW 
                          + 1.0/9.0*cW4*(5.0 - 1.0/2.0/RW)/RW2 ) );
    }
    return Sigma;
}


complex EWSMOneLoopEW::SigmaZZ_fer(const double mu, const double s, 
                                   const double Mw_i) const
{
    double ml2[6], mq2[6];
    for (int i=0; i<6; i++) { 
        ml2[i] = cache.ml2((StandardModel::lepton) i);
        mq2[i] = cache.mq2((StandardModel::quark) i, mu);
    }
    double mu2 = mu*mu;
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_ml2_ml2[6], Bf_s_mq2_mq2[6];
    complex B0_s_ml2_ml2[6], B0_s_mq2_mq2[6];
    if (mu==Mz && s==Mz2) {
        for (int i=0; i<6; i++) {
            Bf_s_ml2_ml2[i] = cache.Bf_Mz2_Mz2_ml2_ml2((StandardModel::lepton) i);
            Bf_s_mq2_mq2[i] = cache.Bf_Mz2_Mz2_mq2_mq2((StandardModel::quark) i);
            B0_s_ml2_ml2[i] = cache.B0_Mz2_Mz2_ml2_ml2((StandardModel::lepton) i);
            B0_s_mq2_mq2[i] = cache.B0_Mz2_Mz2_mq2_mq2((StandardModel::quark) i);
        }
    } else {
        for (int i=0; i<6; i++) {
            Bf_s_ml2_ml2[i] = cache.getPV().Bf(mu2,s,ml2[i],ml2[i]);
            Bf_s_mq2_mq2[i] = cache.getPV().Bf(mu2,s,mq2[i],mq2[i]);
            B0_s_ml2_ml2[i] = cache.getPV().B0(mu2,s,ml2[i],ml2[i]);
            B0_s_mq2_mq2[i] = cache.getPV().B0(mu2,s,mq2[i],mq2[i]);
        }
    }
    
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW::SigmaZZ_fer(s=0.0)"); 
    } else {
        double mf2, vf2, af2;
        for (int i=0; i<6; i++) {
            mf2 = ml2[i];
            vf2 = pow(cache.vl((StandardModel::lepton) i, Mw), 2.0);
            af2 = pow(cache.al((StandardModel::lepton) i), 2.0);
            if(s!=0.0) Sigma += - (vf2 + af2)*s*Bf_s_ml2_ml2[i];
            Sigma += - 2.0*af2*mf2*B0_s_ml2_ml2[i];
            //
            mf2 = mq2[i];
            vf2 = pow(cache.vq((StandardModel::quark) i, Mw), 2.0);
            af2 = pow(cache.aq((StandardModel::quark) i), 2.0);
            if(s!=0.0) Sigma += - 3.0*(vf2 + af2)*s*Bf_s_mq2_mq2[i];
            Sigma += - 3.0*2.0*af2*mf2*B0_s_mq2_mq2[i];
        }
    }   
    return Sigma;
}


complex EWSMOneLoopEW::PiGammaGamma_bos(const double mu, const double s,
                                        const double Mw_i) const 
{
    double mu2 = mu*mu;
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;
    double RW = pow(Mw, 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;    
    
    /* Loop functions */
    double A0_Mw2;
    complex B0_s_Mw2_Mw2;
    if (mu==Mz && s==Mz2) {
        A0_Mw2 = cache.A0_Mz2_Mw2(Mw);
        B0_s_Mw2_Mw2 = cache.B0_Mz2_Mz2_Mw2_Mw2(Mw);
    } else {
        A0_Mw2 = cache.getPV().A0(mu2, Mw2);
        B0_s_Mw2_Mw2 = cache.getPV().B0(mu2, s, Mw2, Mw2);
    }
    
    complex Pi(0.0,0.0,false);
    if (s==0.0) {
        Pi = 7.0*log(Mw2/mu/mu) - 2.0/3.0;
    } else {
        Pi = - RW*( (4.0 + 17.0/3.0/RW - 4.0/3.0/RW2 - 1.0/12.0/RW3)*B0_s_Mw2_Mw2
                    + (4.0 - 4.0/3.0/RW - 1.0/6.0/RW2)*(A0_Mw2/Mw2 + 1.0)
                    - 1.0/18.0/RW2*(1.0/RW - 13.0) );
    }
    return Pi;
}


complex EWSMOneLoopEW::PiGammaGamma_fer_l(const double mu, const double s, 
                                          const StandardModel::lepton l) const 
{
    // Neutrinos do not contribute, since Qf=0.
    if ( (l==StandardModel::NEUTRINO_1) || (l==StandardModel::NEUTRINO_2)
            || (l==StandardModel::NEUTRINO_3) )
        return 0.0;

    double mu2 = mu*mu;
    double mf2 = cache.ml2(l);
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_mf2_mf2;
    if (mu==Mz && s==Mz2) {
        if (mf2==0.0) {
            Bf_s_mf2_mf2 = 0.0;
        } else {
            Bf_s_mf2_mf2 = cache.Bf_Mz2_Mz2_ml2_ml2(l);
        }
    } else if (mu==Mz && s==0.0) {        
        if (mf2==0.0) {
            Bf_s_mf2_mf2 = 0.0;
        } else {
            Bf_s_mf2_mf2 = cache.Bf_Mz2_0_ml2_ml2(l);
        }
    } else {
        if (mf2==0.0) {
            Bf_s_mf2_mf2 = 0.0;
        } else {
            Bf_s_mf2_mf2 = cache.getPV().Bf(mu2,s,mf2,mf2);
        }
    }
    
    double Qf = cache.Ql(l);
    return ( - 4.0*Qf*Qf*Bf_s_mf2_mf2);
}


complex EWSMOneLoopEW::PiGammaGamma_fer_q(const double mu, const double s, 
                                          const StandardModel::quark q) const 
{
    double mu2 = mu*mu;
    double mf2 = cache.mq2(q, mu);
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_mf2_mf2;
    if (mu==Mz && s==Mz2) {
        if (mf2==0.0) {
            Bf_s_mf2_mf2 = 0.0;
        } else {
            Bf_s_mf2_mf2 = cache.Bf_Mz2_Mz2_mq2_mq2(q);
        }
    } else if (mu==Mz && s==0.0) {        
        if (mf2==0.0) {
            Bf_s_mf2_mf2 = 0.0;
        } else {
            Bf_s_mf2_mf2 = cache.Bf_Mz2_0_mq2_mq2(q);
        }
    } else {
        if (mf2==0.0) {
            Bf_s_mf2_mf2 = 0.0;
        } else {
            Bf_s_mf2_mf2 = cache.getPV().Bf(mu2,s,mf2,mf2);
        }
    }
    
    double Qf = cache.Qq(q);
    return ( - 4.0*3.0*Qf*Qf*Bf_s_mf2_mf2);
}


complex EWSMOneLoopEW::PiGammaGamma_fer(const double mu, const double s) const 
{
    complex Pi(0.0,0.0,false);
    for (int i=0; i<6; i++) {
        Pi += PiGammaGamma_fer_l(mu, s, (StandardModel::lepton) i);
        Pi += PiGammaGamma_fer_q(mu, s, (StandardModel::quark) i);        
    }
    return Pi;
}


complex EWSMOneLoopEW::PiZgamma_bos(const double mu, const double s,
                                    const double Mw_i) const
{
    double Mw = cache.Mw(Mw_i);
    double cW2 = cache.cW2(Mw);
    return ( PiGammaGamma_bos(mu,s,Mw)*cW2 );
}


complex EWSMOneLoopEW::PiZgamma_fer(const double mu, const double s,
                                    const double Mw_i) const
{
    double ml2[6], mq2[6];
    for (int i=0; i<6; i++) { 
        ml2[i] = cache.ml2((StandardModel::lepton) i);
        mq2[i] = cache.mq2((StandardModel::quark) i, mu);
    }
    double mu2 = mu*mu;
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    
    /* Loop functions */
    complex Bf_s_ml2_ml2[6], Bf_s_mq2_mq2[6];
    if (mu==Mz && s==Mz2) {
        for (int i=0; i<6; i++) {
            if (i==0 || i==2 || i==4 )
                Bf_s_ml2_ml2[i] = 0.0; // Neutrinos do not contribute, since Ql=0.
            else
                Bf_s_ml2_ml2[i] = cache.Bf_Mz2_Mz2_ml2_ml2((StandardModel::lepton) i);
            Bf_s_mq2_mq2[i] = cache.Bf_Mz2_Mz2_mq2_mq2((StandardModel::quark) i);
        }
    } else {
        for (int i=0; i<6; i++) {
            if (i==0 || i==2 || i==4 )
                Bf_s_ml2_ml2[i] = 0.0; // Neutrinos do not contribute, since Ql=0.
            else    
                Bf_s_ml2_ml2[i] = cache.getPV().Bf(mu2,s,ml2[i],ml2[i]);
            Bf_s_mq2_mq2[i] = cache.getPV().Bf(mu2,s,mq2[i],mq2[i]);
        }
    }

    complex Pi(0.0,0.0,false);
    double Ql, Qq;
    for (int i=0; i<6; i++) {
        Ql = cache.Ql((StandardModel::lepton) i);
        Pi += - (fabs(Ql) - 4.0*sW2*Ql*Ql)*Bf_s_ml2_ml2[i];
        //
        Qq = cache.Qq((StandardModel::quark) i);
        Pi += - 3.0*(fabs(Qq) - 4.0*sW2*Qq*Qq)*Bf_s_mq2_mq2[i];
    }   
    return Pi;
}
 

//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW::SigmaPrime_WW_bos_Mw2(const double mu, 
                                             const double Mw_i) const 
{
    double mu2 = mu*mu;
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;
    double mh = cache.mh();
    double mh2 = mh*mh;
    double sW2 = cache.sW2(Mw_i);
    double cW2 = cache.cW2(Mw_i);
    double cW4 = cW2*cW2;
    double rw = mh2/Mw2;
    //double rz = mh2/Mz2;
    
    /* Loop functions */
    double A0_Mw2, A0_Mz2, A0_mh2;
    complex B0_Mw2_Mz2_Mw2, B0_Mw2_0_Mw2, B0_Mw2_mh2_Mw2;
    complex B0p_Mw2_Mz2_Mw2, B0p_Mw2_0_Mw2, B0p_Mw2_mh2_Mw2;
    if (mu==Mw) {
        A0_Mw2 = cache.A0_Mw2_Mw2(Mw);
        A0_Mz2 = cache.A0_Mw2_Mz2(Mw);
        A0_mh2 = cache.A0_Mw2_mh2(Mw);
        B0_Mw2_Mz2_Mw2 = cache.B0_Mw2_Mw2_Mz2_Mw2(Mw);
        B0_Mw2_0_Mw2 = cache.B0_Mw2_Mw2_0_Mw2(Mw);
        B0_Mw2_mh2_Mw2 = cache.B0_Mw2_Mw2_mh2_Mw2(Mw);
        B0p_Mw2_Mz2_Mw2 = cache.B0p_Mw2_Mw2_Mz2_Mw2(Mw);
        B0p_Mw2_0_Mw2 = cache.B0p_Mw2_Mw2_0_Mw2(Mw);
        B0p_Mw2_mh2_Mw2 = cache.B0p_Mw2_Mw2_mh2_Mw2(Mw);
    } else {
        A0_Mw2 = cache.getPV().A0(mu2, Mw2);
        A0_Mz2 = cache.getPV().A0(mu2, Mz2);
        A0_mh2 = cache.getPV().A0(mu2, mh2);
        B0_Mw2_Mz2_Mw2 = cache.getPV().B0(mu2, Mw2, Mz2, Mw2);
        B0_Mw2_0_Mw2 = cache.getPV().B0(mu2, Mw2, 0.0, Mw2);
        B0_Mw2_mh2_Mw2 = cache.getPV().B0(mu2, Mw2, mh, Mw2);
        B0p_Mw2_Mz2_Mw2 = cache.getPV().B0p(mu2, Mw2, Mz2, Mw2);
        B0p_Mw2_0_Mw2 = cache.getPV().B0p(mu2, Mw2, 0.0, Mw2);
        B0p_Mw2_mh2_Mw2 = cache.getPV().B0p(mu2, Mw2, mh2, Mw2);
    }
    
    complex Sigma(0.0,0.0,false);
    Sigma = - (1.0/12.0/cW4 + 2.0/3.0/cW2 + 2.0*cW2)*B0_Mw2_Mz2_Mw2
            + (1.0/12.0/cW4 + 4.0/3.0/cW2 - 17.0/3.0 - 4.0*cW2)*Mw2*B0p_Mw2_Mz2_Mw2
            - 2.0*sW2*B0_Mw2_0_Mw2 - 4.0*sW2*Mw2*B0p_Mw2_0_Mw2
            + rw/6.0*(1.0 - rw/2.0)*B0_Mw2_mh2_Mw2
            + (1.0 - rw/3.0 + rw*rw/12.0)*Mw2*B0p_Mw2_mh2_Mw2
            + (1.0/cW2 + 8.0 + rw)/12.0*A0_Mw2/Mw2
            - (1.0/cW2 + 9.0 - 8.0*cW2 - 12.0*cW4)/12.0*A0_Mz2/Mw2
            - (rw - 1.0)/12.0*A0_mh2/Mw2 + 4.0/9.0;
    return Sigma;    
}


complex EWSMOneLoopEW::SigmaPrime_WW_fer_Mw2(const double mu,
                                             const double Mw_i) const 
{
    double ml2[6], mq2[6];
    for (int i=0; i<6; i++) { 
        ml2[i] = cache.ml2((StandardModel::lepton) i);
        mq2[i] = cache.mq2((StandardModel::quark) i, mu);
    }
    double mu2 = mu*mu;
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;

    /* Loop functions */
    complex Bf_Mw2_mlprime2_ml2[3], Bf_Mw2_mqprime2_mq2[3];
    complex Bfp_Mw2_mlprime2_ml2[3], Bfp_Mw2_mqprime2_mq2[3];
    complex B1p_Mw2_mlprime2_ml2[3], B1p_Mw2_mqprime2_mq2[3];
    complex B1p_Mw2_ml2_mlprime2[3], B1p_Mw2_mq2_mqprime2[3];
    if (mu==Mw) {
        for (int gen=0; gen<3; gen++) {
            Bf_Mw2_mlprime2_ml2[gen] = cache.Bf_Mw2_Mw2_mlprime2_ml2(gen,Mw);
            Bf_Mw2_mqprime2_mq2[gen] = cache.Bf_Mw2_Mw2_mqprime2_mq2(gen,Mw);
            Bfp_Mw2_mlprime2_ml2[gen] = cache.Bfp_Mw2_Mw2_mlprime2_ml2(gen,Mw);
            Bfp_Mw2_mqprime2_mq2[gen] = cache.Bfp_Mw2_Mw2_mqprime2_mq2(gen,Mw);
            B1p_Mw2_ml2_mlprime2[gen] = cache.B1p_Mw2_Mw2_ml2_mlprime2(gen,Mw);
            B1p_Mw2_mq2_mqprime2[gen] = cache.B1p_Mw2_Mw2_mq2_mqprime2(gen,Mw);
            B1p_Mw2_mlprime2_ml2[gen] = cache.B1p_Mw2_Mw2_mlprime2_ml2(gen,Mw);
            B1p_Mw2_mqprime2_mq2[gen] = cache.B1p_Mw2_Mw2_mqprime2_mq2(gen,Mw);
        }
    } else {
        for (int gen=0; gen<3; gen++) {
            Bf_Mw2_mlprime2_ml2[gen] = cache.getPV().Bf(mu2,Mw2,ml2[2*gen+1],ml2[2*gen]);
            Bf_Mw2_mqprime2_mq2[gen] = cache.getPV().Bf(mu2,Mw2,mq2[2*gen+1],mq2[2*gen]);
            Bfp_Mw2_mlprime2_ml2[gen] = cache.getPV().Bfp(mu2,Mw2,ml2[2*gen+1],ml2[2*gen]);
            Bfp_Mw2_mqprime2_mq2[gen] = cache.getPV().Bfp(mu2,Mw2,mq2[2*gen+1],mq2[2*gen]);
            B1p_Mw2_ml2_mlprime2[gen] = cache.getPV().B1p(mu2,Mw2,ml2[2*gen],ml2[2*gen+1]);
            B1p_Mw2_mq2_mqprime2[gen] = cache.getPV().B1p(mu2,Mw2,mq2[2*gen],mq2[2*gen+1]);
            B1p_Mw2_mlprime2_ml2[gen] = cache.getPV().B1p(mu2,Mw2,ml2[2*gen+1],ml2[2*gen]);
            B1p_Mw2_mqprime2_mq2[gen] = cache.getPV().B1p(mu2,Mw2,mq2[2*gen+1],mq2[2*gen]);
        }        
    }

    complex Sigma(0.0,0.0,false);
    double mf2, mfprime2;
    for (int gen=0; gen<3; gen++) {
        mf2 = ml2[2*gen];
        mfprime2 = ml2[2*gen+1];
        Sigma += - (Bf_Mw2_mlprime2_ml2[gen] + Mw2*Bfp_Mw2_mlprime2_ml2[gen]);
        Sigma += mfprime2*B1p_Mw2_ml2_mlprime2[gen] + mf2*B1p_Mw2_mlprime2_ml2[gen];
        //
        mf2 = mq2[2*gen];
        mfprime2 = mq2[2*gen+1];
        Sigma += - 3.0*(Bf_Mw2_mqprime2_mq2[gen] + Mw2*Bfp_Mw2_mqprime2_mq2[gen]);
        Sigma += 3.0*( mfprime2*B1p_Mw2_mq2_mqprime2[gen] + mf2*B1p_Mw2_mqprime2_mq2[gen] );
    }
    return Sigma;    
}


complex EWSMOneLoopEW::SigmaPrime_ZZ_bos_Mz2(const double mu,
                                             const double Mw_i) const 
{
    double mu2 = mu*mu;
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
    double A0_Mw2, A0_Mz2, A0_mh2;
    complex B0_Mz2_Mw2_Mw2, B0_Mz2_mh2_Mz2;
    complex B0p_Mz2_Mw2_Mw2, B0p_Mz2_mh2_Mz2;
    if (mu==Mz) {
        A0_Mw2 = cache.A0_Mz2_Mw2(Mw);
        A0_Mz2 = cache.A0_Mz2_Mz2();
        A0_mh2 = cache.A0_Mz2_mh2();
        B0_Mz2_Mw2_Mw2 = cache.B0_Mz2_Mz2_Mw2_Mw2(Mw);
        B0_Mz2_mh2_Mz2 = cache.B0_Mz2_Mz2_mh2_Mz2();
        B0p_Mz2_Mw2_Mw2 = cache.B0p_Mz2_Mz2_Mw2_Mw2(Mw);
        B0p_Mz2_mh2_Mz2 = cache.B0p_Mz2_Mz2_mh2_Mz2();
    } else {
        A0_Mw2 = cache.getPV().A0(mu2, Mw2);
        A0_Mz2 = cache.getPV().A0(mu2, Mz2);
        A0_mh2 = cache.getPV().A0(mu2, mh2);
        B0_Mz2_Mw2_Mw2 = cache.getPV().B0(mu2, Mz2, Mw2, Mw2);
        B0_Mz2_mh2_Mz2 = cache.getPV().B0(mu2, Mz2, mh2, Mz2);
        B0p_Mz2_Mw2_Mw2 = cache.getPV().B0p(mu2, Mz2, Mw2, Mw2);
        B0p_Mz2_mh2_Mz2 = cache.getPV().B0p(mu2, Mz2, mh2, Mz2);
    }

    complex Sigma(0.0,0.0,false);
    Sigma = (1.0/4.0/cW2 + 8.0/3.0 - 17.0/3.0*cW2)*B0_Mz2_Mw2_Mw2
            + (1.0/12.0/cW2 + 4.0/3.0 - 17.0/3.0*cW2 - 4.0*cW4)*Mz2*B0p_Mz2_Mw2_Mw2
            + rw/6.0*(1.0 - rz/2.0)*B0_Mz2_mh2_Mz2
            + (1.0 - rz/3.0 + rz*rz/12.0)*Mz2/cW2*B0p_Mz2_mh2_Mz2
            + (1.0 + 4.0*cW2)/3.0/cW2*A0_Mw2/Mz2
            - (1.0 - rz)/12.0/cW2*(A0_Mz2 - A0_mh2)/Mz2
            + 2.0/9.0/cW2 - 10.0/9.0 + 4.0/3.0*cW2;
    Sigma *= cW2;
    return Sigma;        
}


complex EWSMOneLoopEW::SigmaPrime_ZZ_fer_Mz2(const double mu, const double Mw_i) const 
{
    double ml2[6], mq2[6];
    for (int i=0; i<6; i++) { 
        ml2[i] = cache.ml2((StandardModel::lepton) i);
        mq2[i] = cache.mq2((StandardModel::quark) i, mu);
    }
    double mu2 = mu*mu;
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_Mz2_ml2_ml2[6], Bf_Mz2_mq2_mq2[6];
    complex Bfp_Mz2_ml2_ml2[6], Bfp_Mz2_mq2_mq2[6];
    complex B0p_Mz2_ml2_ml2[6], B0p_Mz2_mq2_mq2[6];
    if (mu==Mz) {
         for (int i=0; i<6; i++) {
             Bf_Mz2_ml2_ml2[i] = cache.Bf_Mz2_Mz2_ml2_ml2((StandardModel::lepton) i);
             Bf_Mz2_mq2_mq2[i] = cache.Bf_Mz2_Mz2_mq2_mq2((StandardModel::quark) i);
             Bfp_Mz2_ml2_ml2[i] = cache.Bfp_Mz2_Mz2_ml2_ml2((StandardModel::lepton) i);
             Bfp_Mz2_mq2_mq2[i] = cache.Bfp_Mz2_Mz2_mq2_mq2((StandardModel::quark) i);
             B0p_Mz2_ml2_ml2[i] = cache.B0p_Mz2_Mz2_ml2_ml2((StandardModel::lepton) i);
             B0p_Mz2_mq2_mq2[i] = cache.B0p_Mz2_Mz2_mq2_mq2((StandardModel::quark) i);
         }        
    } else {
        for (int i=0; i<6; i++) {
            Bf_Mz2_ml2_ml2[i] = cache.getPV().Bf(mu2,Mz2,ml2[i],ml2[i]);
            Bf_Mz2_mq2_mq2[i] = cache.getPV().Bf(mu2,Mz2,mq2[i],mq2[i]);
            Bfp_Mz2_ml2_ml2[i] = cache.getPV().Bfp(mu2,Mz2,ml2[i],ml2[i]);
            Bfp_Mz2_mq2_mq2[i] = cache.getPV().Bfp(mu2,Mz2,mq2[i],mq2[i]);
            B0p_Mz2_ml2_ml2[i] = cache.getPV().B0p(mu2,Mz2,ml2[i],ml2[i]);
            B0p_Mz2_mq2_mq2[i] = cache.getPV().B0p(mu2,Mz2,mq2[i],mq2[i]);
        }        
   }
    
    complex Sigma(0.0,0.0,false);
    double mf2, vf2, af2;
    for (int i=0; i<6; i++) {
        mf2 = ml2[i];
        vf2 = pow(cache.vl((StandardModel::lepton) i, Mw), 2.0);
        af2 = pow(cache.al((StandardModel::lepton) i), 2.0);
        Sigma += - (vf2 + af2)*(Bf_Mz2_ml2_ml2[i] + Mz2*Bfp_Mz2_ml2_ml2[i])
                 - 2.0*af2*mf2*B0p_Mz2_ml2_ml2[i];
        //
        mf2 = mq2[i];
        vf2 = pow(cache.vq((StandardModel::quark) i, Mw), 2.0);
        af2 = pow(cache.aq((StandardModel::quark) i), 2.0);
        Sigma += - 3.0*(vf2 + af2)*(Bf_Mz2_mq2_mq2[i] + Mz2*Bfp_Mz2_mq2_mq2[i])
                 - 6.0*af2*mf2*B0p_Mz2_mq2_mq2[i];
    }
    return Sigma;    
}


//////////////////////////////////////////////////////////////////////// 

double EWSMOneLoopEW::DeltaRhobar(const double mu, const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();    
    return ( (SigmaWW_bos(mu,Mw*Mw,Mw).real() + SigmaWW_fer(mu,Mw*Mw,Mw).real() 
              - SigmaZZ_bos(mu,Mz*Mz,Mw).real() - SigmaZZ_fer(mu,Mz*Mz,Mw).real())
             /Mw/Mw );
}


double EWSMOneLoopEW::DeltaRhobarW(const double mu, const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    return ( (SigmaWW_bos(mu,0.0,Mw).real() + SigmaWW_fer(mu,0.0,Mw).real() 
              - SigmaWW_bos(mu,Mw*Mw,Mw).real() - SigmaWW_fer(mu,Mw*Mw,Mw).real())
             /Mw/Mw );
}


//////////////////////////////////////////////////////////////////////// 

double EWSMOneLoopEW::TEST_DeltaRhobar_bos(const double Mw_i) const 
{
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
    complex B0_Mw2_Mz2_Mw2_Mw2 = cache.B0_Mz2_Mz2_Mw2_Mw2(Mw) + log_cW2;
    complex B0_Mw2_Mz2_mh2_Mz2 = cache.B0_Mz2_Mz2_mh2_Mz2() + log_cW2;
    complex B0_Mw2_Mw2_Mz2_Mw2 = cache.B0_Mz2_Mw2_Mz2_Mw2(Mw) + log_cW2;
    complex B0_Mw2_Mw2_mh2_Mw2 = cache.B0_Mz2_Mw2_mh2_Mw2(Mw) + log_cW2;
    
    double DRhobar;    
    DRhobar = - (1.0/12.0/cW2 + 4.0/3.0 - 17.0/3.0*cW2 - 4.0*cW4)
                 *(B0_Mw2_Mz2_Mw2_Mw2.real() - 1.0/cW2*B0_Mw2_Mw2_Mz2_Mw2.real())
              + (1.0 - 1.0/3.0*rw + 1.0/12.0*rw*rw)*B0_Mw2_Mw2_mh2_Mw2.real()
              - (1.0 - 1.0/3.0*rz + 1.0/12.0*rz*rz)/cW2*B0_Mw2_Mz2_mh2_Mz2.real()
              + 1.0/12.0*sW2*rw*rw*(log(rw) - 1.0)
              - (1.0/12.0/cW4 + 1.0/2.0/cW2 - 2.0 + 1.0/12.0*rw)*log_cW2
              - 1.0/12.0/cW4 - 19.0/36.0/cW2 - 133.0/18.0 + 8.0*cW2;
    return DRhobar;
} 


double EWSMOneLoopEW::TEST_DeltaRhobarW_bos(const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double sW2 = cache.sW2(Mw);
    double cW2 = cache.cW2(Mw);
    double cW4 = cW2*cW2;
    double mh = cache.mh();
    double rw = pow(mh/Mw, 2.0);

    /* Logarithm and one-loop functions */
    double log_cW2 = cache.log_cW2(Mw); 
    
    /* B0 functions for mu=Mw */
    complex B0_Mw2_Mz2_Mw2 = cache.B0_Mz2_Mw2_Mz2_Mw2(Mw) + log_cW2;
    complex B0_Mw2_mh2_Mw2 = cache.B0_Mz2_Mw2_mh2_Mw2(Mw) + log_cW2;
    
    double DRhobarW;    
    DRhobarW = - (1.0/12.0/cW4 + 4.0/3.0/cW2 - 17.0/3.0 - 4.0*cW2)*B0_Mw2_Mz2_Mw2.real()
               - (1.0 - 1.0/3.0*rw + 1.0/12.0*rw*rw)*B0_Mw2_mh2_Mw2.real()
               + (3.0/4.0/(1.0-rw) + 1.0/4.0 - 1.0/12.0*rw)*rw*log(rw)
               + (1.0/12.0/cW4 + 17.0/12.0/cW2 - 3.0/sW2 + 1.0/4.0)*log_cW2
               + 1.0/12.0/cW4 + 11.0/8.0/cW2 + 139.0/36.0 - 177.0/24.0*cW2 
               + 5.0/8.0*cW4 - 1.0/12.0*rw*(7.0/2.0 - rw);     
    return DRhobarW;
} 

   
//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW::FZa_0(const double s, const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();   
    double Rz = Mz*Mz/s;

    /* Logarithm and three-point one-loop functions */
    double log_Rz;
    complex C0_s_0_Mz2_0;
    if (s==Mz*Mz) {
        log_Rz = 0.0;
        C0_s_0_Mz2_0 = cache.C0_Mz2_0_Mz2_0();
    } else if (s==Mw*Mw) {
        log_Rz = - cache.log_cW2(Mw);
        C0_s_0_Mz2_0 = cache.C0_Mw2_0_Mz2_0(Mw);
    } else {
        log_Rz = log(Rz);
        C0_s_0_Mz2_0 = cache.getPV().C0(s,0.0,Mz*Mz,0.0);
    }
        
    complex FZa(0.0,0.0,false);
    FZa = 2.0*pow((Rz + 1.0), 2.0)*s*C0_s_0_Mz2_0
          - (2.0*Rz + 3.0)*(log_Rz + M_PI*complex::i()) - 2.0*Rz - 7.0/2.0;
    return FZa;    
}


complex EWSMOneLoopEW::FWa_0(const double s, const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;

    /* Logarithm and three-point one-loop functions */
    double log_Rw;
    complex C0_s_0_Mw2_0;
    if (s==Mz*Mz) {    
        log_Rw = cache.log_cW2(Mw);
        C0_s_0_Mw2_0 = cache.C0_Mz2_0_Mw2_0(Mw);
    } else {
        log_Rw = log(Rw);
        C0_s_0_Mw2_0 = cache.getPV().C0(s,0.0,Mw*Mw,0.0);
    }
    
    complex FWa(0.0,0.0,false);
    FWa = 2.0*pow((Rw + 1.0), 2.0)*s*C0_s_0_Mw2_0
          - (2.0*Rw + 3.0)*(log_Rw + M_PI*complex::i()) - 2.0*Rw - 7.0/2.0;
    return FWa;     
}


complex EWSMOneLoopEW::FbarWa_0(const double s) const 
{
    return ( complex(0.0,0.0,false) );
}


complex EWSMOneLoopEW::FWn_0(const double s, const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;

    /* Two- and three-point one-loop functions */
    complex B0_Mw_s_Mw2_Mw2;
    complex C0_s_Mw2_0_Mw2;
    if (s==Mz*Mz) {
        B0_Mw_s_Mw2_Mw2 = cache.B0_Mw2_Mz2_Mw2_Mw2(Mw);
        C0_s_Mw2_0_Mw2 = cache.C0_Mz2_Mw2_0_Mw2(Mw);
    } else {    
        B0_Mw_s_Mw2_Mw2 = cache.getPV().B0(Mw2,s,Mw2,Mw2);
        C0_s_Mw2_0_Mw2 = cache.getPV().C0(s,Mw2,0.0,Mw2);
    }
 
    complex FWn(0.0,0.0,false);    
    FWn = - 2.0*(Rw + 2.0)*Mw*Mw*C0_s_Mw2_0_Mw2
          - (2.0*Rw + 7.0/3.0 - 3.0/2.0/Rw - 1.0/12.0/Rw/Rw)*B0_Mw_s_Mw2_Mw2
          + 2.0*Rw + 9.0/2.0 - 11.0/18.0/Rw + 1.0/18.0/Rw/Rw;  
    return FWn;
}
    

complex EWSMOneLoopEW::FWa_t(const double s, const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;
    double Mt = cache.Mt();
    double Mt2 = Mt*Mt;
    double wt = Mt2/Mw2;

    /* Logarithm and two- and three-point one-loop functions */
    double log_wt = - 2.0*cache.logMZtoMTOP() - cache.log_cW2(Mw);
    double log_Rw;
    complex B0_Mw2_s_Mt2_Mt2;
    complex C0_s_Mt2_Mw2_Mt2, C0_s_0_Mw2_0;
    if (s==Mz*Mz) {
        log_Rw = cache.log_cW2(Mw);
        B0_Mw2_s_Mt2_Mt2 = cache.B0_Mw2_Mz2_Mt2_Mt2(Mw);
        C0_s_Mt2_Mw2_Mt2 = cache.C0_Mz2_Mt2_Mw2_Mt2(Mw);
        C0_s_0_Mw2_0 = cache.C0_Mz2_0_Mw2_0(Mw);
    } else {
        log_Rw = log(Rw);
        B0_Mw2_s_Mt2_Mt2 = cache.getPV().B0(Mw2,s,Mt2,Mt2);
        C0_s_Mt2_Mw2_Mt2 = cache.getPV().C0(s,Mt2,Mw2,Mt2);
        C0_s_0_Mw2_0 = cache.getPV().C0(s,0.0,Mw2,0.0);
    }
    
    complex FWa(0.0,0.0,false);    
    FWa = 2.0*(Rw + 1.0)*(Rw + 1.0)*s*(C0_s_Mt2_Mw2_Mt2 - C0_s_0_Mw2_0)
          + (2.0*Rw + 3.0)*(- B0_Mw2_s_Mt2_Mt2 + log_Rw + M_PI*complex::i() + 2.0)
          - wt*( (3.0*Rw + 2.0 - wt - wt*wt*Rw)*Mw*Mw*C0_s_Mt2_Mw2_Mt2
                 + (Rw + 1.0/2.0 + wt*Rw)*(1.0 - B0_Mw2_s_Mt2_Mt2)
                 - (2.0*Rw + 1.0/2.0 - 2.0/(wt - 1.0) 
                    + 3.0/2.0/(wt - 1.0)/(wt - 1.0) + wt*Rw)*log_wt
                 + 3.0/2.0/(wt - 1.0) + 3.0/4.0 );
    return FWa;
}


complex EWSMOneLoopEW::FbarWa_t(const double s, const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;
    double Mt = cache.Mt();
    double Mt2 = Mt*Mt;
    double wt = Mt2/Mw2;

    /* Logarithm and two- and three-point one-loop functions */
    double log_wt = - 2.0*cache.logMZtoMTOP() - cache.log_cW2(Mw);
    complex B0_Mw2_s_Mt2_Mt2;
    complex C0_s_Mt2_Mw2_Mt2;
    if (s==Mz*Mz) {
        B0_Mw2_s_Mt2_Mt2 = cache.B0_Mw2_Mz2_Mt2_Mt2(Mw);
        C0_s_Mt2_Mw2_Mt2 = cache.C0_Mz2_Mt2_Mw2_Mt2(Mw);
    } else {
        B0_Mw2_s_Mt2_Mt2 = cache.getPV().B0(Mw2,s,Mt2,Mt2);
        C0_s_Mt2_Mw2_Mt2 = cache.getPV().C0(s,Mt2,Mw2,Mt2);
    }
    
    complex FbarWa(0.0,0.0,false);        
    FbarWa = - wt*( (Rw + 2.0 - wt*(2.0 - wt)*Rw)*Mw*Mw*C0_s_Mt2_Mw2_Mt2
                    - (1.0/2.0 - Rw + wt*Rw)*(- B0_Mw2_s_Mt2_Mt2 + 1.0)
                    + wt*Rw*log_wt );
    return FbarWa;
}


complex EWSMOneLoopEW::FWn_t(const double s, const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Mz = cache.Mz();   
    double Rw = Mw*Mw/s;
    double Mt = cache.Mt();
    double Mt2 = Mt*Mt;
    double wt = Mt2/Mw2;

    /* Logarithm and two- and three-point one-loop functions */
    double log_wt = - 2.0*cache.logMZtoMTOP() - cache.log_cW2(Mw);  
    complex B0_Mw2_s_Mw2_Mw2;
    complex C0_s_Mw2_Mt2_Mw2, C0_s_Mw2_0_Mw2;
    if (s==Mz*Mz) {
        B0_Mw2_s_Mw2_Mw2 = cache.B0_Mw2_Mz2_Mw2_Mw2(Mw);
        C0_s_Mw2_Mt2_Mw2 = cache.C0_Mz2_Mw2_Mt2_Mw2(Mw);
        C0_s_Mw2_0_Mw2 = cache.C0_Mz2_Mw2_0_Mw2(Mw);
    } else {
        B0_Mw2_s_Mw2_Mw2 = cache.getPV().B0(Mw2,s,Mw2,Mw2);
        C0_s_Mw2_Mt2_Mw2 = cache.getPV().C0(s,Mw2,Mt2,Mw2);
        C0_s_Mw2_0_Mw2 = cache.getPV().C0(s,Mw2,0.0,Mw2);
    }
    
    complex FWn(0.0,0.0,false);        
    FWn = - 2.0*(Rw + 2.0)*Mw*Mw*(C0_s_Mw2_Mt2_Mw2 - C0_s_Mw2_0_Mw2)
          + wt*( (3.0*Rw + 5.0/2.0 - 2.0/Rw - wt*(2.0 - 1.0/2.0/Rw) 
                  + wt*wt*(1.0/2.0 - Rw))*Mw*Mw*C0_s_Mw2_Mt2_Mw2
                  + (Rw + 1.0 - 1.0/4.0/Rw - wt*(1.0/2.0 - Rw))*(B0_Mw2_s_Mw2_Mw2 - 1.0)
                  + (2.0*Rw + 1.0/2.0 - 2.0/(wt - 1.0) 
                     + 3.0/2.0/(wt - 1.0)/(wt - 1.0) - wt*(1.0/2.0 - Rw))*log_wt
                  - 3.0/2.0/(wt- 1.0) + 1.0/4.0 
                  - 1.0/2.0/Rw );
    return FWn;
}


complex EWSMOneLoopEW::FZ(const double s, const double Mw_i) const 
{
    return ( FZa_0(s, Mw_i) );
}


complex EWSMOneLoopEW::FW_l(const double s, const StandardModel::lepton l, 
                            const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double cW2 = cache.cW2(Mw);

    StandardModel::lepton lprime;
    switch(l) {
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
            throw std::runtime_error("Error in EWSMOneLoopEW::FW_l()");   
    }
    return ( cW2*FWn_0(s, Mw) - cache.sigmal(lprime, Mw)/2.0*FWa_0(s, Mw) 
            - FbarWa_0(s)/2.0 );
}
   

complex EWSMOneLoopEW::FW_q(const double s, const StandardModel::quark q, 
                            const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double cW2 = cache.cW2(Mw);

    StandardModel::quark qprime;
    switch(q) {
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
            throw std::runtime_error("TOP is not allowed in EWSMOneLoopEW::FW_q(s,q)"); 
        case StandardModel::BOTTOM:
            qprime = StandardModel::TOP;
            break;
        default:
            throw std::runtime_error("Error in EWSMOneLoopEW::FW_q()");   
    }         
    complex FW(0.0,0.0,false);
    FW = cW2*FWn_0(s,Mw) - cache.sigmaq(qprime, Mw)/2.0*FWa_0(s, Mw) 
            - FbarWa_0(s)/2.0;
    if (q==StandardModel::BOTTOM) {
        FW += cW2*FWn_t(s,Mw) - cache.sigmaq(qprime, Mw)/2.0*FWa_t(s, Mw) 
                - FbarWa_t(s, Mw)/2.0;
    }
    return FW;
}


//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW::TEST_FWn(const double s, const double mf, 
                                const double Mw_i) const 
{
    double Mw = cache.Mw(Mw_i);
    double Mw2 = Mw*Mw;
    double Rw = Mw2/s;
    double mf2 = mf*mf;
    double wf = mf2/Mw2;

    /* Logarithm and two- and three-point one-loop functions */
    double log_wf = log(wf);  
    double A0_Mw2 = cache.getPV().A0(Mw2, Mw2);
    double A0_mf2 = cache.getPV().A0(Mw2, mf2);
    complex B0_Mw2_s_Mw2_Mw2 = cache.getPV().B0(Mw2,s,Mw2,Mw2);
    complex B0_Mw2_0_mf2_Mw2 = cache.getPV().B0(Mw2,0.0,mf2,Mw2);
    complex C0_s_Mw2_mf2_Mw2 = cache.getPV().C0(s,Mw2,mf2,Mw2);
    complex C0_s_Mw2_0_Mw2 = cache.getPV().C0(s,Mw2,0.0,Mw2);
    
    complex FWn(0.0,0.0,false);        
    /* Eq.(5.586) in Bardin and Passarino's book */
    FWn = ((2.0 + wf)*(1.0 - wf)*(1.0 - wf)*Rw + 4.0 - 5.0/2.0*wf + 2.0*wf*wf
            - wf*wf*wf/2.0 + wf*(2.0 - wf/2.0)/Rw)*Mw*Mw*C0_s_Mw2_mf2_Mw2
           - (- (2.0 + wf)*(1.0 - wf)*Rw - 3.0 + 3.0/2.0*wf - wf*wf/2.0)
             *(B0_Mw2_s_Mw2_Mw2 - B0_Mw2_0_mf2_Mw2)
           - (2.0/3.0 - wf/2.0 + (3.0/2.0 - wf/4.0)/Rw + 1.0/12.0/Rw/Rw)*B0_Mw2_s_Mw2_Mw2
           - A0_mf2/Mw/Mw - (2.0/3.0 + 1.0/6.0/Rw)*A0_Mw2/Mw/Mw
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





