/* 
 * File:   OneLoopEW.cpp
 * Author: mishima
 */

#include "OneLoopEW.h"


OneLoopEW::OneLoopEW(const EWSMcommon& EWSMC_i) : EWSMC(EWSMC_i) {
}

//OneLoopEW::OneLoopEW(const OneLoopEW& orig) {
//}

OneLoopEW::~OneLoopEW() {
}


////////////////////////////////////////////////////////////////////////

double OneLoopEW::DeltaAlpha_l() const {   
    double xl[3] = { pow(EWSMC.GetSM().getMz()
                     /EWSMC.GetSM().getLeptons(EWSMC.GetSM().ELECTRON).getMass(), 2.0), 
                     pow(EWSMC.GetSM().getMz()
                     /EWSMC.GetSM().getLeptons(EWSMC.GetSM().MU).getMass(), 2.0), 
                     pow(EWSMC.GetSM().getMz()
                     /EWSMC.GetSM().getLeptons(EWSMC.GetSM().TAU).getMass(), 2.0) };
    double log_l[3] = { 2.0*EWSMC.GetLogMZtoME(), 
                        2.0*EWSMC.GetLogMZtoMMU(), 
                        2.0*EWSMC.GetLogMZtoMTAU() };  

    double oneLoop[3];
    for (int i = 0; i < 3; i++) {
        oneLoop[i] = - 5.0/9.0 + log_l[i]/3.0 - 2.0/xl[i];
    }
            
    return( EWSMC.GetSM().getAle()/M_PI
            *(oneLoop[0] + oneLoop[1] + oneLoop[2]) );
}

double OneLoopEW::DeltaAlpha_t() const {   
    double xt = pow(EWSMC.GetSM().getMz()
                    /EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass(), 2.0);
    double tmp = 1.0 + xt*0.1071;
    tmp *= -4.0/45.0*EWSMC.GetSM().getAle()/M_PI*xt;
    return tmp;
}

double OneLoopEW::DeltaRho() const {
    double Mz = EWSMC.GetSM().getMz();
    return ( - EWSMC.GetSM().getAle()/4.0/M_PI/EWSMC.GetSW2()*DeltaRhobar(Mz) );
}

double OneLoopEW::DeltaR_rem() const {
    double Mz = EWSMC.GetSM().getMz();
    double Mz2 = Mz*Mz;
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    
    /* Logarithm */
    double log_cW2 = EWSMC.GetLog_cW2();
    
    double PiGammaGamma_t_0 = PiGammaGamma_fer(Mz,0.0,StandardModel::TOP).real();
    double PiGammaGamma_l5q_Mz2 = PiGammaGamma_fer(Mz,Mz2).real() 
                                  - PiGammaGamma_fer(Mz,Mz2,StandardModel::TOP).real();

    double DR_rem = - 2.0/3.0*sW2 + sW2*PiGammaGamma_t_0
                    + sW2*PiGammaGamma_l5q_Mz2 + DeltaRhobarW(Mz)
                    + (4.0 - 25.0/4.0*cW2 + 3.0/4.0*cW2*cW2 + 9.0*cW2/4.0/sW2)*log_cW2
                    + 11.0/2.0 - 5.0/8.0*cW2*(1.0 + cW2);
    DR_rem *= EWSMC.GetSM().getAle()/4.0/M_PI/sW2;
    return DR_rem;    
}

double OneLoopEW::DeltaRbar_rem() const {
    double Mz = EWSMC.GetSM().getMz();
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();

    /* Logarithm */
    double log_cW2 = EWSMC.GetLog_cW2();    
    
    double DRbar_rem = - 2.0/3.0*sW2 + DeltaRhobarW(Mz)
                       + (4.0 - 25.0/4.0*cW2 + 3.0/4.0*cW2*cW2 + 9.0*cW2/4.0/sW2)*log_cW2
                       + 11.0/2.0 - 5.0/8.0*cW2*(1.0 + cW2);
    DRbar_rem *= EWSMC.GetSM().getAle()/4.0/M_PI/sW2;
    return DRbar_rem;     
}

complex OneLoopEW::deltaRho_rem_tmp(const complex uf) const {
    double Mz = EWSMC.GetSM().getMz();  
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();

    /* Logarithm */
    double log_cW2 = EWSMC.GetLog_cW2(); 

    complex dRho_rem(0.0,0.0,false);
    dRho_rem = - ( SigmaPrime_ZZ_bos_Mz2(Mz).real() 
                   + SigmaPrime_ZZ_fer_Mz2(Mz).real() )/cW2
               - DeltaRhobarW(Mz) + 2.0*uf
               - (1.0/6.0/cW2 - 1.0/3.0 + 3.0/4.0*cW2*(1.0 + cW2) + 9.0*cW2/4.0/sW2)
                 *log_cW2
               - 11.0/2.0 + 5.0/8.0*cW2*(1.0 + cW2);
    dRho_rem *= EWSMC.GetSM().getAle()/4.0/M_PI/sW2;
    return dRho_rem;  
}

complex OneLoopEW::deltaRho_rem_l(const StandardModel::lepton l) const {
    double Mz = EWSMC.GetSM().getMz(); 
    complex uf = ( 3.0*EWSMC.vf(l)*EWSMC.vf(l) + EWSMC.af(l)*EWSMC.af(l) )
                 /4.0/EWSMC.GetCW2()*FZ(Mz*Mz) + FW(Mz*Mz, l);
    return ( deltaRho_rem_tmp(uf) );   
}

complex OneLoopEW::deltaRho_rem_q(const StandardModel::quark q) const {
    double Mz = EWSMC.GetSM().getMz(); 
    complex uf = ( 3.0*EWSMC.vf(q)*EWSMC.vf(q) + EWSMC.af(q)*EWSMC.af(q) )
                 /4.0/EWSMC.GetCW2()*FZ(Mz*Mz) + FW(Mz*Mz, q);
    return ( deltaRho_rem_tmp(uf) );      
}

complex OneLoopEW::deltaKappa_rem_tmp(const double deltaf, const complex uf) const {
    double Mz = EWSMC.GetSM().getMz();  
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();

    /* Logarithm */
    double log_cW2 = EWSMC.GetLog_cW2(); 

    complex dKappa_rem(0.0,0.0,false);
    dKappa_rem = - ( PiZgamma_bos(Mz,Mz*Mz) + PiZgamma_fer(Mz,Mz*Mz) )
                 + deltaf*deltaf/4.0/cW2*FZ(Mz*Mz) - uf
                 + (1.0/12.0/cW2 + 4.0/3.0)*log_cW2;
    dKappa_rem *= EWSMC.GetSM().getAle()/4.0/M_PI/sW2;
    return dKappa_rem;    
}

complex OneLoopEW::deltaKappa_rem_l(const StandardModel::lepton l) const {
    double Mz = EWSMC.GetSM().getMz(); 
    complex uf = ( 3.0*EWSMC.vf(l)*EWSMC.vf(l) + EWSMC.af(l)*EWSMC.af(l) )
                 /4.0/EWSMC.GetCW2()*FZ(Mz*Mz) + FW(Mz*Mz, l);
    return ( deltaKappa_rem_tmp(EWSMC.deltaf(l), uf) );       
}

complex OneLoopEW::deltaKappa_rem_q(const StandardModel::quark q) const {
    double Mz = EWSMC.GetSM().getMz(); 
    complex uf = ( 3.0*EWSMC.vf(q)*EWSMC.vf(q) + EWSMC.af(q)*EWSMC.af(q) )
                 /4.0/EWSMC.GetCW2()*FZ(Mz*Mz) + FW(Mz*Mz, q);
    return ( deltaKappa_rem_tmp(EWSMC.deltaf(q), uf) );     
}


//////////////////////////////////////////////////////////////////////// 

complex OneLoopEW::SigmaWW_bos(const double mu, const double s) const {
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    double Mz = EWSMC.GetSM().getMz();    
    //double Mz2 = Mz*Mz;
    double mh = EWSMC.GetSM().getMHl();
    double mh2 = mh*mh;
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double RW = pow(EWSMC.GetMw(), 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;
    double rw = pow(EWSMC.GetSM().getMHl()/EWSMC.GetMw(), 2.0);
    
    /* Loop functions */
    double A0_Mw, A0_Mz, A0_mh;
    complex B0_s_Mz_Mw, B0_s_0_Mw, B0_s_mh_Mw;
    complex B0p_s_Mz_Mw, B0p_s_mh_Mw; /* for s==0.0 */
    if (mu==Mz && s==0.0) {
        A0_Mw = EWSMC.GetA0_Mz_Mw();
        A0_Mz = EWSMC.GetA0_Mz_Mz();
        A0_mh = EWSMC.GetA0_Mz_mh();
        B0_s_Mz_Mw = EWSMC.GetB0_Mz_0_Mz_Mw();
        B0_s_0_Mw = EWSMC.GetB0_Mz_0_0_Mw();
        B0_s_mh_Mw = EWSMC.GetB0_Mz_0_mh_Mw();
        B0p_s_Mz_Mw = EWSMC.GetB0p_Mz_0_Mz_Mw();
        B0p_s_mh_Mw = EWSMC.GetB0p_Mz_0_mh_Mw();         
    } else if (mu==Mz && s==Mw2) {
        A0_Mw = EWSMC.GetA0_Mz_Mw();
        A0_Mz = EWSMC.GetA0_Mz_Mz();
        A0_mh = EWSMC.GetA0_Mz_mh();
        B0_s_Mz_Mw = EWSMC.GetB0_Mz_Mw2_Mz_Mw();
        B0_s_0_Mw = EWSMC.GetB0_Mz_Mw2_0_Mw();
        B0_s_mh_Mw = EWSMC.GetB0_Mz_Mw2_mh_Mw();
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

complex OneLoopEW::SigmaWW_fer(const double mu, const double s) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = EWSMC.GetSM().getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = EWSMC.GetSM().getQuarks((StandardModel::quark) i).getMass();
    }
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    double Mz = EWSMC.GetSM().getMz();        

    /* Loop functions */
    complex B1_s_ml_mlprime[3], B1_s_mq_mqprime[3]; 
    complex B1_s_mlprime_ml[3], B1_s_mqprime_mq[3];    
    complex Bf_s_mlprime_ml[3], Bf_s_mqprime_mq[3];
    if (mu==Mz && s==0.0) {
        for (int gen=0; gen<3; gen++) {
            B1_s_ml_mlprime[gen] = EWSMC.GetB1_Mz_0_ml_mlprime(gen);
            B1_s_mq_mqprime[gen] = EWSMC.GetB1_Mz_0_mq_mqprime(gen);
            B1_s_mlprime_ml[gen] = EWSMC.GetB1_Mz_0_mlprime_ml(gen);
            B1_s_mqprime_mq[gen] = EWSMC.GetB1_Mz_0_mqprime_mq(gen);
            Bf_s_mlprime_ml[gen] = EWSMC.GetBf_Mz_0_mlprime_ml(gen);
            Bf_s_mqprime_mq[gen] = EWSMC.GetBf_Mz_0_mqprime_mq(gen);            
        }
    } else if (mu==Mz && s==Mw2) {
        for (int gen=0; gen<3; gen++) {
            B1_s_ml_mlprime[gen] = EWSMC.GetB1_Mz_Mw2_ml_mlprime(gen);
            B1_s_mq_mqprime[gen] = EWSMC.GetB1_Mz_Mw2_mq_mqprime(gen);   
            B1_s_mlprime_ml[gen] = EWSMC.GetB1_Mz_Mw2_mlprime_ml(gen);
            B1_s_mqprime_mq[gen] = EWSMC.GetB1_Mz_Mw2_mqprime_mq(gen);           
            Bf_s_mlprime_ml[gen] = EWSMC.GetBf_Mz_Mw2_mlprime_ml(gen);
            Bf_s_mqprime_mq[gen] = EWSMC.GetBf_Mz_Mw2_mqprime_mq(gen);           
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

complex OneLoopEW::SigmaZZ_bos(const double mu, const double s) const {
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;
    double mh = EWSMC.GetSM().getMHl();
    double mh2 = mh*mh;
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double RW = pow(EWSMC.GetMw(), 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;
    double rw = pow(EWSMC.GetSM().getMHl()/EWSMC.GetMw(), 2.0);
    
    /* Loop functions */
    double A0_Mw, A0_Mz, A0_mh;
    complex B0_s_Mw_Mw, B0_s_mh_Mz;
    if (mu==Mz && s==Mz2) {
        A0_Mw = EWSMC.GetA0_Mz_Mw();
        A0_Mz = EWSMC.GetA0_Mz_Mz();
        A0_mh = EWSMC.GetA0_Mz_mh();
        B0_s_Mw_Mw = EWSMC.GetB0_Mz_Mz2_Mw_Mw();
        B0_s_mh_Mz = EWSMC.GetB0_Mz_Mz2_mh_Mz();
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        A0_Mw = myPV->A0(mu, Mw);
        A0_Mz = myPV->A0(mu, Mz);
        A0_mh = myPV->A0(mu, mh);
        B0_s_Mw_Mw = myPV->B0(mu, s, Mw, Mw);
        B0_s_mh_Mz = myPV->B0(mu, s, EWSMC.GetSM().getMHl(), Mz);
        delete myPV;
    }        

    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw "Missing codes for OneLoopEW::SigmaZZ_bos(s=0.0)";
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

complex OneLoopEW::SigmaZZ_fer(const double mu, const double s) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = EWSMC.GetSM().getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = EWSMC.GetSM().getQuarks((StandardModel::quark) i).getMass();
    }
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_ml_ml[6], Bf_s_mq_mq[6];
    complex B0_s_ml_ml[6], B0_s_mq_mq[6];    
    if (mu==Mz && s==Mz2) {
        for (int i=0; i<6; i++) {
            Bf_s_ml_ml[i] = EWSMC.GetBf_Mz_Mz2_ml_ml((StandardModel::lepton) i);
            Bf_s_mq_mq[i] = EWSMC.GetBf_Mz_Mz2_mq_mq((StandardModel::quark) i);           
            B0_s_ml_ml[i] = EWSMC.GetB0_Mz_Mz2_ml_ml((StandardModel::lepton) i);
            B0_s_mq_mq[i] = EWSMC.GetB0_Mz_Mz2_mq_mq((StandardModel::quark) i);           
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
        throw "Missing codes for OneLoopEW::SigmaZZ_fer(s=0.0)";
    } else {
        double ml2, vl2, al2, mq2, vq2, aq2;
        for (int i=0; i<6; i++) {
            ml2 = ml[i]*ml[i];
            vl2 = pow(EWSMC.vf((StandardModel::lepton) i), 2.0);
            al2 = pow(EWSMC.af((StandardModel::lepton) i), 2.0);            
            if(s!=0.0) Sigma += - (vl2 + al2)*s*Bf_s_ml_ml[i];
            Sigma += - 2.0*al2*ml2*B0_s_ml_ml[i];
            //
            mq2 = mq[i]*mq[i];
            vq2 = pow(EWSMC.vf((StandardModel::quark) i), 2.0);
            aq2 = pow(EWSMC.af((StandardModel::quark) i), 2.0);
            if(s!=0.0) Sigma += - 3.0*(vq2 + aq2)*s*Bf_s_mq_mq[i];
            Sigma += - 3.0*2.0*aq2*mq2*B0_s_mq_mq[i];
        }
    }   
 
    /* added O(alpha^2) contribution from the Z-gamma mixing */
    //double Mw2 = pow(EWSMC.GetMw(), 2.0);
    //Sigma += - Mw2*pow(PiZgamma_fer(mu,Mz2), 2.0)*EWSMC.GetSM().getAle()/4.0/M_PI;
    
    return Sigma;
}

complex OneLoopEW::PiGammaGamma_bos(const double mu, const double s) const {
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    //double Mz = EWSMC.GetSM().getMz();    
    //double Mz2 = Mz*Mz;
    //double cW2 = EWSMC.GetCW2();
    double RW = pow(EWSMC.GetMw(), 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;    
    
    /* Loop functions */
    double A0_Mw;
    complex B0_s_Mw_Mw;
    PVfunctions* myPV;
    myPV = new PVfunctions();
    A0_Mw = myPV->A0(mu, Mw);
    B0_s_Mw_Mw = myPV->B0(mu, s, Mw, Mw);
    delete myPV; 
    
    complex Pi(0.0,0.0,false);
    if (s==0.0) {
        throw "Missing codes for OneLoopEW::PiGammaGamma_bos(s=0.0)";
    } else {
        Pi = - RW*( (4.0 + 17.0/3.0/RW - 4.0/3.0/RW2 - 1.0/12.0/RW3)*B0_s_Mw_Mw
                    + (4.0 - 4.0/3.0/RW - 1.0/6.0/RW2)*(A0_Mw/Mw2 + 1.0)
                    - 1.0/18.0/RW2*(1.0/RW - 13.0) );
    }
    return Pi;
}

complex OneLoopEW::PiGammaGamma_fer(const double mu, const double s, 
                                    const StandardModel::lepton l) const {
    double ml = EWSMC.GetSM().getLeptons(l).getMass();
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_ml_ml;
    if (mu==Mz && s==Mz2) {
        if (ml==0.0) { 
            Bf_s_ml_ml = 0.0; 
        } else {
            Bf_s_ml_ml = EWSMC.GetBf_Mz_Mz2_ml_ml(l);
        }
    } else {
        if (mu==Mz && s==0.0) {        
            if (ml==0.0) {
                Bf_s_ml_ml = 0.0; 
            } else {
                Bf_s_ml_ml = EWSMC.GetBf_Mz_0_ml_ml(l);
            }
        } else {
            PVfunctions* myPV;
            myPV = new PVfunctions();
            if (ml==0.0) {
                Bf_s_ml_ml = 0.0; 
            } else {
                Bf_s_ml_ml = myPV->Bf(mu,s,ml,ml);
            }
            delete myPV;
        }
    }
    
    double Ql = EWSMC.Qf(l);
    return ( - 4.0*Ql*Ql*Bf_s_ml_ml);
}

complex OneLoopEW::PiGammaGamma_fer(const double mu, const double s, 
                                    const StandardModel::quark q) const {
    double mq = EWSMC.GetSM().getQuarks(q).getMass();
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_mq_mq;
    if (mu==Mz && s==Mz2) {
        if (mq==0.0) { 
            Bf_s_mq_mq = 0.0; 
        } else {
            Bf_s_mq_mq = EWSMC.GetBf_Mz_Mz2_mq_mq(q);
        }
    } else {
        if (mu==Mz && s==0.0) {        
            if (mq==0.0) {
                Bf_s_mq_mq = 0.0; 
            } else {
                Bf_s_mq_mq = EWSMC.GetBf_Mz_0_mq_mq(q);
            }
        } else {
            PVfunctions* myPV;
            myPV = new PVfunctions();
            if (mq==0.0) {
                Bf_s_mq_mq = 0.0; 
            } else {
                Bf_s_mq_mq = myPV->Bf(mu,s,mq,mq);
            }
            delete myPV;
        }
    }
    
    double Qq = EWSMC.Qf(q);
    return ( - 4.0*3.0*Qq*Qq*Bf_s_mq_mq);
}

complex OneLoopEW::PiGammaGamma_fer(const double mu, const double s) const {
    complex Pi(0.0,0.0,false);
    for (int i=0; i<6; i++) {
        Pi += PiGammaGamma_fer(mu, s, (StandardModel::lepton) i);
        Pi += PiGammaGamma_fer(mu, s, (StandardModel::quark) i);        
    }
    return Pi;
}

complex OneLoopEW::PiZgamma_bos(const double mu, const double s) const {
    double cW2 = EWSMC.GetCW2();
    return ( - PiGammaGamma_bos(mu,s)*cW2);
}

complex OneLoopEW::PiZgamma_fer(const double mu, const double s) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = EWSMC.GetSM().getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = EWSMC.GetSM().getQuarks((StandardModel::quark) i).getMass();
    }
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_ml_ml[6], Bf_s_mq_mq[6];
    if (mu==Mz && s==Mz2) {
        for (int i=0; i<6; i++) {
            Bf_s_ml_ml[i] = EWSMC.GetBf_Mz_Mz2_ml_ml((StandardModel::lepton) i);
            Bf_s_mq_mq[i] = EWSMC.GetBf_Mz_Mz2_mq_mq((StandardModel::quark) i);           
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
        Ql = EWSMC.Qf((StandardModel::lepton) i);
        Pi += (fabs(Ql) - 4.0*EWSMC.GetSW2()*Ql*Ql)*Bf_s_ml_ml[i];
        //
        Qq = EWSMC.Qf((StandardModel::quark) i);
        Pi += 3.0*(fabs(Qq) - 4.0*EWSMC.GetSW2()*Qq*Qq)*Bf_s_mq_mq[i];
    }   
    return Pi;
}
 

//////////////////////////////////////////////////////////////////////// 

complex OneLoopEW::SigmaPrime_WW_bos_Mw2(const double mu) const {
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    double Mz = EWSMC.GetSM().getMz();    
    //double Mz2 = Mz*Mz;
    double mh = EWSMC.GetSM().getMHl();
    double mh2 = mh*mh;
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double rw = mh2/Mw2;
    //double rz = mh2/Mz2;
    
    /* Loop functions */
    double A0_Mw, A0_Mz, A0_mh;
    complex B0_Mw2_Mz_Mw, B0_Mw2_0_Mw, B0_Mw2_mh_Mw;
    complex B0p_Mw2_Mz_Mw, B0p_Mw2_0_Mw, B0p_Mw2_mh_Mw;
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
    
    complex Sigma(0.0,0.0,false);
    Sigma = - (1.0/12.0/cW4 + 2.0/3.0/cW2 + 2.0*cW2)*B0_Mw2_Mz_Mw
            + (1.0/12.0/cW4 + 4.0/3.0/cW2 - 17.0/3.0 - 4.0*cW2)*Mw2*B0p_Mw2_Mz_Mw
            - 29.0/6.0*sW2*B0_Mw2_0_Mw - 4.0*sW2*Mw2*B0p_Mw2_0_Mw
            + rw/6.0*(1.0 - rw/2.0)*B0_Mw2_mh_Mw
            + (1.0 - rw/3.0 + rw*rw/12.0)*Mw2*B0p_Mw2_0_Mw
            + (1.0/cW2 + 8.0 + rw)/12.0*A0_Mw/Mw2
            - (1.0/cW2 + 9.0 - 8.0*cW2 - 12.0*cW4)/12.0*A0_Mz/Mw2
            - (rw - 1.0)/12.0*A0_mh/Mw2 + 4.0/9.0;
    return Sigma;    
}

complex OneLoopEW::SigmaPrime_WW_fer_Mw2(const double mu) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = EWSMC.GetSM().getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = EWSMC.GetSM().getQuarks((StandardModel::quark) i).getMass();
    }
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    //double Mz = EWSMC.GetSM().getMz();        

    /* Loop functions */
    complex Bf_Mw2_mlprime_ml[3], Bf_Mw2_mqprime_mq[3];
    complex Bfp_Mw2_mlprime_ml[3], Bfp_Mw2_mqprime_mq[3];    
    complex B1p_Mw2_mlprime_ml[3], B1p_Mw2_mqprime_mq[3];  
    complex B1p_Mw2_ml_mlprime[3], B1p_Mw2_mq_mqprime[3];    
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

complex OneLoopEW::SigmaPrime_ZZ_bos_Mz2(const double mu) const {
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;
    double mh = EWSMC.GetSM().getMHl();
    double mh2 = mh*mh;
    //double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double rw = mh2/Mw2;
    double rz = mh2/Mz2;
    
    /* Loop functions */
    double A0_Mw, A0_Mz, A0_mh;
    complex B0_Mz2_Mw_Mw, B0_Mz2_mh_Mz;
    complex B0p_Mz2_Mw_Mw, B0p_Mz2_mh_Mz;
    if (mu==Mz) {
        A0_Mw = EWSMC.GetA0_Mz_Mw();
        A0_Mz = EWSMC.GetA0_Mz_Mz();
        A0_mh = EWSMC.GetA0_Mz_mh();
        B0_Mz2_Mw_Mw = EWSMC.GetB0_Mz_Mz2_Mw_Mw();
        B0_Mz2_mh_Mz = EWSMC.GetB0_Mz_Mz2_mh_Mz();
        B0p_Mz2_Mw_Mw = EWSMC.GetB0p_Mz_Mz2_Mw_Mw();
        B0p_Mz2_mh_Mz = EWSMC.GetB0p_Mz_Mz2_mh_Mz();
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

complex OneLoopEW::SigmaPrime_ZZ_fer_Mz2(const double mu) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = EWSMC.GetSM().getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = EWSMC.GetSM().getQuarks((StandardModel::quark) i).getMass();
    }
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_Mz2_ml_ml[6], Bf_Mz2_mq_mq[6];
    complex Bfp_Mz2_ml_ml[6], Bfp_Mz2_mq_mq[6];
    complex B0p_Mz2_ml_ml[6], B0p_Mz2_mq_mq[6]; 
    if (mu==Mz) {
         for (int i=0; i<6; i++) {
             Bf_Mz2_ml_ml[i] = EWSMC.GetBf_Mz_Mz2_ml_ml((StandardModel::lepton) i);
             Bf_Mz2_mq_mq[i] = EWSMC.GetBf_Mz_Mz2_mq_mq((StandardModel::quark) i);           
             Bfp_Mz2_ml_ml[i] = EWSMC.GetBfp_Mz_Mz2_ml_ml((StandardModel::lepton) i);
             Bfp_Mz2_mq_mq[i] = EWSMC.GetBfp_Mz_Mz2_mq_mq((StandardModel::quark) i);           
             B0p_Mz2_ml_ml[i] = EWSMC.GetB0p_Mz_Mz2_ml_ml((StandardModel::lepton) i);
             B0p_Mz2_mq_mq[i] = EWSMC.GetB0p_Mz_Mz2_mq_mq((StandardModel::quark) i);           
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
        vl2 = pow(EWSMC.vf((StandardModel::lepton) i), 2.0);
        al2 = pow(EWSMC.af((StandardModel::lepton) i), 2.0);            
        Sigma += - (vl2 + al2)*(Bf_Mz2_ml_ml[i] + Mz2*Bfp_Mz2_ml_ml[i])
                 - 2.0*al2*ml2*B0p_Mz2_ml_ml[i];
        //
        mq2 = mq[i]*mq[i];
        vq2 = pow(EWSMC.vf((StandardModel::quark) i), 2.0);
        aq2 = pow(EWSMC.af((StandardModel::quark) i), 2.0);
        Sigma += - 3.0*(vq2 + aq2)*(Bf_Mz2_mq_mq[i] + Mz2*Bfp_Mz2_mq_mq[i])
                 - 6.0*aq2*mq2*B0p_Mz2_mq_mq[i];
    }
    return Sigma;    
}


//////////////////////////////////////////////////////////////////////// 

double OneLoopEW::DeltaRhobar(const double mu) const {
    double Mw = EWSMC.GetMw();
    double Mz = EWSMC.GetSM().getMz();    
    return ( (SigmaWW_bos(mu,Mw*Mw).real() + SigmaWW_fer(mu,Mw*Mw).real() 
              - SigmaZZ_bos(mu,Mz*Mz).real() - SigmaZZ_fer(mu,Mz*Mz).real())
             /Mw/Mw );
}

double OneLoopEW::DeltaRhobarW(const double mu) const {
    double Mw = EWSMC.GetMw();
    return ( (SigmaWW_bos(mu,0.0).real() + SigmaWW_fer(mu,0.0).real() 
              - SigmaWW_bos(mu,Mw*Mw).real() - SigmaWW_fer(mu,Mw*Mw).real())
             /Mw/Mw );
}


//////////////////////////////////////////////////////////////////////// 

double OneLoopEW::TEST_DeltaRhobar_bos() const {
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double rw = pow(EWSMC.GetSM().getMHl()/EWSMC.GetMw(), 2.0);
    double rz = pow(EWSMC.GetSM().getMHl()/EWSMC.GetSM().getMz(), 2.0);
    
    /* B0 functions for mu=Mw */
    complex B0_Mw_Mz2_Mw_Mw = EWSMC.GetB0_Mz_Mz2_Mw_Mw() + log(cW2);
    complex B0_Mw_Mz2_mh_Mz = EWSMC.GetB0_Mz_Mz2_mh_Mz() + log(cW2);
    complex B0_Mw_Mw2_Mz_Mw = EWSMC.GetB0_Mz_Mw2_Mz_Mw() + log(cW2);
    complex B0_Mw_Mw2_mh_Mw = EWSMC.GetB0_Mz_Mw2_mh_Mw() + log(cW2);
    
    double DRhobar;    
    DRhobar = - (1.0/12.0/cW2 + 4.0/3.0 - 17.0/3.0*cW2 - 4.0*cW4)
                 *(B0_Mw_Mz2_Mw_Mw.real() - 1.0/cW2*B0_Mw_Mw2_Mz_Mw.real())
              + (1.0 - 1.0/3.0*rw + 1.0/12.0*rw*rw)*B0_Mw_Mw2_mh_Mw.real()
              - (1.0 - 1.0/3.0*rz + 1.0/12.0*rz*rz)/cW2*B0_Mw_Mz2_mh_Mz.real()
              + 1.0/12.0*sW2*rw*rw*(log(rw) - 1.0)
              - (1.0/12.0/cW4 + 1.0/2.0/cW2 - 2.0 + 1.0/12.0*rw)*log(cW2)
              - 1.0/12.0/cW4 - 19.0/36.0/cW2 - 133.0/18.0 + 8.0*cW2;
    return DRhobar;
} 

double OneLoopEW::TEST_DeltaRhobarW_bos() const {
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double rw = pow(EWSMC.GetSM().getMHl()/EWSMC.GetMw(), 2.0);
    
    /* B0 functions for mu=Mw */
    complex B0_Mw2_Mz_Mw = EWSMC.GetB0_Mz_Mw2_Mz_Mw() + log(cW2);
    complex B0_Mw2_mh_Mw = EWSMC.GetB0_Mz_Mw2_mh_Mw() + log(cW2);
    
    double DRhobarW;    
    DRhobarW = - (1.0/12.0/cW4 + 4.0/3.0/cW2 - 17.0/3.0 - 4.0*cW2)*B0_Mw2_Mz_Mw.real()
               - (1.0 - 1.0/3.0*rw + 1.0/12.0*rw*rw)*B0_Mw2_mh_Mw.real()
               + (3.0/4.0/(1.0-rw) + 1.0/4.0 - 1.0/12.0*rw)*rw*log(rw)
               + (1.0/12.0/cW4 + 17.0/12.0/cW2 - 3.0/sW2 + 1.0/4.0)*log(cW2)
               + 1.0/12.0/cW4 + 11.0/8.0/cW2 + 139.0/36.0 - 177.0/24.0*cW2 
               + 5.0/8.0*cW4 - 1.0/12.0*rw*(7.0/2.0 - rw);     
    return DRhobarW;
} 

   
//////////////////////////////////////////////////////////////////////// 

complex OneLoopEW::FZa_0(const double s) const {
    double Mz = EWSMC.GetSM().getMz();   
    double Rz = Mz*Mz/s;

    /* Three-point one-loop functions */
    complex C0_s_0_Mz_0;
    if (s==Mz*Mz) {
        C0_s_0_Mz_0 = EWSMC.GetC0_Mz2_0_Mz_0();
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        C0_s_0_Mz_0 = myPV->C0(s,0.0,Mz,0.0);
        delete myPV;
    }
    
    /* Logarithm */
    double log_Rz;
    if (s==Mz*Mz) { 
        log_Rz = 0.0;
    } else {
        log_Rz = log(Rz);
    }
    
    complex FZa(0.0,0.0,false);
    FZa = 2.0*pow((Rz + 1.0), 2.0)*s*C0_s_0_Mz_0
          - (2.0*Rz + 3.0)*(log_Rz + M_PI*complex::i()) - 2.0*Rz - 7.0/2.0;
    return FZa;    
}

complex OneLoopEW::FWa_0(const double s) const {
    double Mw = EWSMC.GetMw();
    double Mz = EWSMC.GetSM().getMz();   
    double Rw = Mw*Mw/s;

    /* Three-point one-loop functions */
    complex C0_s_0_Mw_0;
    if (s==Mz*Mz) {    
        C0_s_0_Mw_0 = EWSMC.GetC0_Mz2_0_Mw_0(); 
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        C0_s_0_Mw_0 = myPV->C0(s,0.0,Mw,0.0);
        delete myPV;    
    }

    /* Logarithm */
    double log_Rw;
    if (s==Mz*Mz) {  
        log_Rw = EWSMC.GetLog_cW2();
    } else {
        log_Rw = log(Rw);
    }
    
    complex FWa(0.0,0.0,false);
    FWa = 2.0*pow((Rw + 1.0), 2.0)*s*C0_s_0_Mw_0
          - (2.0*Rw + 3.0)*(log_Rw + M_PI*complex::i()) - 2.0*Rw - 7.0/2.0;
    return FWa;     
}

complex OneLoopEW::FbarWa_0(const double s) const {
    complex FbarWa(0.0,0.0,false);
    return FbarWa;    
}

complex OneLoopEW::FWn_0(const double s) const {
    double Mw = EWSMC.GetMw();
    double Mz = EWSMC.GetSM().getMz();   
    double Rw = Mw*Mw/s;

    /* Two- and three-point one-loop functions */
    complex B0_Mw_s_Mw_Mw;
    complex C0_s_Mw_0_Mw;
    if (s==Mz*Mz) {
        B0_Mw_s_Mw_Mw = EWSMC.GetB0_Mw_Mz2_Mw_Mw();
        C0_s_Mw_0_Mw = EWSMC.GetC0_Mz2_Mw_0_Mw();
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
    
complex OneLoopEW::FWa_t(const double s) const {
    double Mw = EWSMC.GetMw();
    double Mz = EWSMC.GetSM().getMz();   
    double Rw = Mw*Mw/s;
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();
    double wt = Mt*Mt/Mw/Mw; 

    /* Two- and three-point one-loop functions */
    complex B0_Mw_s_Mt_Mt;
    complex C0_s_Mt_Mw_Mt, C0_s_0_Mw_0;
    if (s==Mz*Mz) {
        B0_Mw_s_Mt_Mt = EWSMC.GetB0_Mw_Mz2_Mt_Mt();
        C0_s_Mt_Mw_Mt = EWSMC.GetC0_Mz2_Mt_Mw_Mt();
        C0_s_0_Mw_0 = EWSMC.GetC0_Mz2_0_Mw_0();
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        B0_Mw_s_Mt_Mt = myPV->B0(Mw,s,Mt,Mt); 
        C0_s_Mt_Mw_Mt = myPV->C0(s,Mt,Mw,Mt);
        C0_s_0_Mw_0 = myPV->C0(s,0.0,Mw,0.0);
        delete myPV;      
    }
    
    /* Logarithm */
    double log_wt = - 2.0*EWSMC.GetLogMZtoMTOP() - EWSMC.GetLog_cW2();
    double log_Rw;
    if (s==Mz*Mz) {
        log_Rw = EWSMC.GetLog_cW2();
    } else {
        log_Rw = log(Rw);
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

complex OneLoopEW::FbarWa_t(const double s) const {
    double Mw = EWSMC.GetMw();
    double Mz = EWSMC.GetSM().getMz();   
    double Rw = Mw*Mw/s;
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();
    double wt = Mt*Mt/Mw/Mw; 

    /* Two- and three-point one-loop functions */
    complex B0_Mw_s_Mt_Mt;
    complex C0_s_Mt_Mw_Mt;    
    if (s==Mz*Mz) {
        B0_Mw_s_Mt_Mt = EWSMC.GetB0_Mw_Mz2_Mt_Mt();
        C0_s_Mt_Mw_Mt = EWSMC.GetC0_Mz2_Mt_Mw_Mt();
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        B0_Mw_s_Mt_Mt = myPV->B0(Mw,s,Mt,Mt);     
        C0_s_Mt_Mw_Mt = myPV->C0(s,Mt,Mw,Mt);    
        delete myPV;   
    }
    
    /* Logarithm */
    double log_wt = - 2.0*EWSMC.GetLogMZtoMTOP() - EWSMC.GetLog_cW2();
    
    complex FbarWa(0.0,0.0,false);        
    FbarWa = - wt*( (Rw + 2.0 - wt*(2.0 - wt)*Rw)*Mw*Mw*C0_s_Mt_Mw_Mt
                    - (1.0/2.0 - Rw + wt*Rw)*(- B0_Mw_s_Mt_Mt + 1.0) 
                    + wt*Rw*log_wt );
    return FbarWa;
}

complex OneLoopEW::FWn_t(const double s) const {
    double Mw = EWSMC.GetMw();
    double Mz = EWSMC.GetSM().getMz();   
    double Rw = Mw*Mw/s;
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();
    double wt = Mt*Mt/Mw/Mw; 

    /* Two- and three-point one-loop functions */
        complex B0_Mw_s_Mw_Mw;
        complex C0_s_Mw_Mt_Mw, C0_s_Mw_0_Mw;
    if (s==Mz*Mz) {
        B0_Mw_s_Mw_Mw = EWSMC.GetB0_Mw_Mz2_Mw_Mw();
        C0_s_Mw_Mt_Mw = EWSMC.GetC0_Mz2_Mw_Mt_Mw();
        C0_s_Mw_0_Mw = EWSMC.GetC0_Mz2_Mw_0_Mw();  
    } else {
        PVfunctions* myPV;
        myPV = new PVfunctions();
        B0_Mw_s_Mw_Mw = myPV->B0(Mw,s,Mw,Mw); 
        C0_s_Mw_Mt_Mw = myPV->C0(s,Mw,Mt,Mw);
        C0_s_Mw_0_Mw = myPV->C0(s,Mw,0.0,Mw);    
        delete myPV;     
    }

    /* Logarithm */
    double log_wt = - 2.0*EWSMC.GetLogMZtoMTOP() - EWSMC.GetLog_cW2();  
    
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

complex OneLoopEW::FZ(const double s) const {
    return ( FZa_0(s) );
}

complex OneLoopEW::FW(const double s, const StandardModel::lepton l) const {
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
            throw "Error in OneLoopEW::FW(s,l)";  
    }
    return ( EWSMC.GetCW2()*FWn_0(s) - EWSMC.sigmaf(lprime)/2.0*FWa_0(s) 
             - FbarWa_0(s)/2.0 );
}

complex OneLoopEW::FW(const double s, const StandardModel::quark q) const {
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
            throw "TOP is not allowed in OneLoopEW::FW(s,q)";
        case StandardModel::BOTTOM:
            qprime = StandardModel::TOP;
            break;
        default:
            throw "Error in OneLoopEW::FW(s,q)";  
    }
    
    complex FW(0.0,0.0,false);
    FW = EWSMC.GetCW2()*FWn_0(s) - EWSMC.sigmaf(qprime)/2.0*FWa_0(s) 
         - FbarWa_0(s)/2.0;
    if (q==StandardModel::BOTTOM) {
        FW += EWSMC.GetCW2()*FWn_t(s) - EWSMC.sigmaf(qprime)/2.0*FWa_t(s) 
              - FbarWa_t(s)/2.0;
    }
    return FW;
}
   

//////////////////////////////////////////////////////////////////////// 

complex OneLoopEW::TEST_FWn(const double s, const double mf) const {
    double Mw = EWSMC.GetMw();
    //double Mz = EWSMC.GetSM().getMz();   
    double Rw = Mw*Mw/s;
    double wf = mf*mf/Mw/Mw; 

    /* Two- and three-point one-loop functions */
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double A0_Mw = myPV->A0(Mw, Mw);
    double A0_mf = myPV->A0(Mw, mf);
    complex B0_Mw_s_Mw_Mw = myPV->B0(Mw,s,Mw,Mw); 
    complex B0_Mw_0_mf_Mw = myPV->B0(Mw,0.0,mf,Mw); 
    complex C0_s_Mw_mf_Mw = myPV->C0(s,Mw,mf,Mw);
    complex C0_s_Mw_0_Mw = myPV->C0(s,Mw,0.0,Mw);    
    delete myPV;     

    /* Logarithm */
    double log_wf = log(wf);  
    
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





