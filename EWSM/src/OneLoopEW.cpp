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
    double Mw2 = pow(EWSMC.GetMw(), 2.0);
    double Mz2 = pow(EWSMC.GetSM().getMz(), 2.0);
    double DeltaRho = ( SigmaWW_bos(Mw2).real() + SigmaWW_fer(Mw2).real() 
                        - SigmaZZ_bos(Mz2).real() - SigmaZZ_fer(Mz2).real() )/Mw2;
    DeltaRho *= - EWSMC.GetSM().getAle()/4.0/M_PI/EWSMC.GetSW2();

    
    std::cout << "SigmaWW_bos(Mw2)/Mw2 = " << SigmaWW_bos(Mw2)/Mw2 << std::endl;
    std::cout << "SigmaZZ_bos(Mz2)/Mw2 = " << SigmaZZ_bos(Mz2)/Mw2 << std::endl;
    std::cout << "PiZgamma_bos = " << PiZgamma_bos(Mz2) << std::endl;
    std::cout << std::endl;

    std::cout << "SigmaWW_fer(Mw2)/Mw2 = " << SigmaWW_fer(Mw2)/Mw2 << std::endl;    
    std::cout << "SigmaZZ_fer(Mz2)/Mw2 = " << SigmaZZ_fer(Mz2)/Mw2 << std::endl;    
    std::cout << "PiZgamma_fer = " << PiZgamma_fer(Mz2) << std::endl;
    std::cout << std::endl;    

    std::cout << "DeltaRhobar_bos = " << SigmaWW_bos(Mw2).real()/Mw2 - SigmaZZ_bos(Mz2).real()/Mw2 << std::endl;
    std::cout << "DeltaRhobarW_bos() = " << DeltaRhobarW_bos() << std::endl;
    std::cout << std::endl;        
    
    std::cout << "DeltaRhobar_fer = " << SigmaWW_fer(Mw2).real()/Mw2 - SigmaZZ_fer(Mz2).real()/Mw2 << std::endl;
    //std::cout << "DeltaRhobarW_fer = " << DeltaRhobarW_bos() << std::endl;
    std::cout << std::endl;  
    
    return DeltaRho;
}

double OneLoopEW::DeltaR_rem() const {
    /* !! Write codes !!*/
    return (0.0);    
}

double OneLoopEW::DeltaRbar_rem() const {
    /* !! Write codes !!*/
    return (0.0);    
}

complex OneLoopEW::deltaRho_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;    
}

complex OneLoopEW::deltaRho_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex OneLoopEW::deltaKappa_rem_l(const StandardModel::lepton l) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;      
}

complex OneLoopEW::deltaKappa_rem_q(const StandardModel::quark q) const {
    /* !! Write codes !!*/
    complex a(0.0,0.0,false);
    return a;    
}


//////////////////////////////////////////////////////////////////////// 

complex OneLoopEW::SigmaWW_bos(const double s) const {
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    double Mz = EWSMC.GetSM().getMz();    
    //double Mz2 = Mz*Mz;
    double mh2 = pow(EWSMC.GetSM().getMHl(), 2.0);
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double RW = pow(EWSMC.GetMw(), 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;
    double rw = pow(EWSMC.GetSM().getMHl()/EWSMC.GetMw(), 2.0);
    
    /* Loop functions */
    double A0_Mw = EWSMC.GetA0_Mw();
    double A0_Mz = EWSMC.GetA0_Mz();
    double A0_mh = EWSMC.GetA0_mh();
    complex B0_s_Mz_Mw, B0_s_0_Mw, B0_s_mh_Mw;

    complex Sigma;
    if (s==0.0) {
       /* This formula is wrong! */
       /* Sigma = Mw2*( 2.0/3.0*(1.0/cW2 - 4.0 - 4.0*cW2 + cW4)
                        *(-A0_Mz + A0_Mw)/Mz2/(1.0 - cW2)
                      + 17.0*sW2/6.0*A0_Mw/Mw2
                      - 1.0/12.0*(- 10.0 + 2.0*rw)*(-A0_mh + A0_Mw)/Mw2/(rw - 1.0)
                      - 1.0/12.0*(24.0 - 2.0*cW2 + cW4)*A0_Mw/Mw2
                      - 1.0/6.0*A0_mh/Mw2
                      - 1.0/12.0*(1.0 + 14.0*cW2 + 9.0*cW4)*A0_Mz/Mw2
                      - 1.0/6.0*(1.0/cW2 + 22.0 + cW2 + cW4 + rw) ); */
        throw "Codes for OneLoopEW::SigmaWW_bos(s=0.0) is missing";
    } else {
        if (s==Mw2) {
            B0_s_Mz_Mw = EWSMC.GetB0_Mw2_Mz_Mw();
            B0_s_0_Mw = EWSMC.GetB0_Mw2_0_Mw();
            B0_s_mh_Mw = EWSMC.GetB0_Mw2_mh_Mw();
        } else {
            PVfunctions* myPV;
            myPV = new PVfunctions();
            B0_s_Mz_Mw = myPV->B0(Mz, s, Mz, Mw);
            B0_s_0_Mw = myPV->B0(Mz, s, 0.0, Mw);
            B0_s_mh_Mw = myPV->B0(Mz, s, EWSMC.GetSM().getMHl(), Mw);
            delete myPV;
        }
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

complex OneLoopEW::SigmaZZ_bos(const double s) const {
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;
    double mh2 = pow(EWSMC.GetSM().getMHl(), 2.0);
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double RW = pow(EWSMC.GetMw(), 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;
    double rw = pow(EWSMC.GetSM().getMHl()/EWSMC.GetMw(), 2.0);
    
    /* Loop functions */
    double A0_Mw = EWSMC.GetA0_Mw();
    double A0_Mz = EWSMC.GetA0_Mz();
    double A0_mh = EWSMC.GetA0_mh();
    complex B0_s_Mw_Mw, B0_s_mh_Mz;

    complex Sigma;
    if (s==0.0) {
        throw "Codes for OneLoopEW::SigmaZZ_bos(s=0.0) is missing";
    } else {
        if (s==Mz2) {
            B0_s_Mw_Mw = EWSMC.GetB0_Mz2_Mw_Mw();
            B0_s_mh_Mz = EWSMC.GetB0_Mz2_mh_Mz();
        } else {
            PVfunctions* myPV;
            myPV = new PVfunctions();
            B0_s_Mw_Mw = myPV->B0(Mz, s, Mw, Mw);
            B0_s_mh_Mz = myPV->B0(Mz, s, EWSMC.GetSM().getMHl(), Mz);
            delete myPV;
        }        
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




complex OneLoopEW::PiGammaGamma_bos(const double s) const {
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;
    //double cW2 = EWSMC.GetCW2();
    double RW = pow(EWSMC.GetMw(), 2.0)/s;
    double RW2 = RW*RW;
    double RW3 = RW2*RW;    
    
    /* Loop functions */
    double A0_Mw = EWSMC.GetA0_Mw();
    complex B0_s_Mw_Mw;
    
    complex Pi;
    if (s==0.0) {
        throw "Codes for OneLoopEW::PiGammaGamma_bos(s=0.0) is missing";
    } else {
        if (s==Mz2) {
            B0_s_Mw_Mw = EWSMC.GetB0_Mz2_Mw_Mw();    
        } else {
            PVfunctions* myPV;
            myPV = new PVfunctions();
            B0_s_Mw_Mw = myPV->B0(Mz, s, Mw, Mw);
            delete myPV;    
        }    
        Pi = RW*( (4.0 + 17.0/3.0/RW - 4.0/3.0/RW2 - 1.0/12.0/RW3)*B0_s_Mw_Mw
                  + (4.0 - 4.0/3.0/RW - 1.0/6.0/RW2)*(A0_Mw/Mw2 + 1.0)
                  - 1.0/18.0/RW2*(1.0/RW - 13.0) );
    }
    return Pi;
}


complex OneLoopEW::PiZgamma_bos(const double s) const {
    double cW2 = EWSMC.GetCW2();
    return (PiGammaGamma_bos(s)*cW2);
}



complex OneLoopEW::SigmaWW_fer(const double s) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = EWSMC.GetSM().getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = EWSMC.GetSM().getQuarks((StandardModel::quark) i).getMass();
    }
    double Mw = EWSMC.GetMw();
    double Mw2 = Mw*Mw;
    double Mz = EWSMC.GetSM().getMz();        

    /* Loop functions */
    complex Bf_s_ml_mlprime[3], B1_s_ml_mlprime[3];
    complex Bf_s_mq_mqprime[3], B1_s_mq_mqprime[3];    
    complex Bf_s_mlprime_ml[3], B1_s_mlprime_ml[3];
    complex Bf_s_mqprime_mq[3], B1_s_mqprime_mq[3];    
    
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw "Codes for OneLoopEW::SigmaWW_fer(s=0.0) is missing";
    } else {
       if (s==Mw2) {
            for (int gen=0; gen<3; gen++) {
                Bf_s_ml_mlprime[gen] = EWSMC.GetBf_Mw2_ml_mlprime(gen);
                Bf_s_mq_mqprime[gen] = EWSMC.GetBf_Mw2_mq_mqprime(gen);           
                B1_s_ml_mlprime[gen] = EWSMC.GetB1_Mw2_ml_mlprime(gen);
                B1_s_mq_mqprime[gen] = EWSMC.GetB1_Mw2_mq_mqprime(gen);   
                Bf_s_mlprime_ml[gen] = EWSMC.GetBf_Mw2_mlprime_ml(gen);
                Bf_s_mqprime_mq[gen] = EWSMC.GetBf_Mw2_mqprime_mq(gen);           
                B1_s_mlprime_ml[gen] = EWSMC.GetB1_Mw2_mlprime_ml(gen);
                B1_s_mqprime_mq[gen] = EWSMC.GetB1_Mw2_mqprime_mq(gen);           
            }
       } else {
            PVfunctions* myPV;
            myPV = new PVfunctions();
            for (int gen=0; gen<3; gen++) {
                Bf_s_ml_mlprime[gen] = myPV->Bf(Mz,s,ml[2*gen],ml[2*gen+1]);
                Bf_s_mq_mqprime[gen] = myPV->Bf(Mz,s,mq[2*gen],mq[2*gen+1]);            
                B1_s_ml_mlprime[gen] = myPV->B1(Mz,s,ml[2*gen],ml[2*gen+1]);
                B1_s_mq_mqprime[gen] = myPV->B1(Mz,s,mq[2*gen],mq[2*gen+1]);
                Bf_s_mlprime_ml[gen] = myPV->Bf(Mz,s,ml[2*gen+1],ml[2*gen]);
                Bf_s_mqprime_mq[gen] = myPV->Bf(Mz,s,mq[2*gen+1],mq[2*gen]);            
                B1_s_mlprime_ml[gen] = myPV->B1(Mz,s,ml[2*gen+1],ml[2*gen]);
                B1_s_mqprime_mq[gen] = myPV->B1(Mz,s,mq[2*gen+1],mq[2*gen]);
            }
            delete myPV;
        }
        
        double ml2, mlprime2, mq2, mqprime2;
        for (int gen=0; gen<3; gen++) {
            ml2 = ml[2*gen]*ml[2*gen];
            mlprime2 = ml[2*gen+1]*ml[2*gen+1];
            Sigma += - s*Bf_s_ml_mlprime[gen];
            Sigma += mlprime2*B1_s_ml_mlprime[gen] + ml2*B1_s_mlprime_ml[gen];
            //
            mq2 = mq[2*gen]*mq[2*gen];
            mqprime2 = mq[2*gen+1]*mq[2*gen+1];
            Sigma += 3.0*( - s*Bf_s_mq_mqprime[gen] );
            Sigma += 3.0*( mqprime2*B1_s_mq_mqprime[gen] + mq2*B1_s_mqprime_mq[gen] );
        }
    }   
    return Sigma;
}

complex OneLoopEW::SigmaZZ_fer(const double s) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = EWSMC.GetSM().getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = EWSMC.GetSM().getQuarks((StandardModel::quark) i).getMass();
    }
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_ml_ml[6], B0_s_ml_ml[6];
    complex Bf_s_mq_mq[6], B0_s_mq_mq[6];    
    
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw "Codes for OneLoopEW::SigmaZZ_fer(s=0.0) is missing";
    } else {
       if (s==Mz2) {
            for (int i=0; i<6; i++) {
                Bf_s_ml_ml[i] = EWSMC.GetBf_Mz2_ml_ml((StandardModel::lepton) i);
                Bf_s_mq_mq[i] = EWSMC.GetBf_Mz2_mq_mq((StandardModel::quark) i);           
                B0_s_ml_ml[i] = EWSMC.GetB0_Mz2_ml_ml((StandardModel::lepton) i);
                B0_s_mq_mq[i] = EWSMC.GetB0_Mz2_mq_mq((StandardModel::quark) i);           
            }
       } else {
            PVfunctions* myPV;
            myPV = new PVfunctions();
            for (int i=0; i<6; i++) {
                Bf_s_ml_ml[i] = myPV->Bf(Mz,s,ml[i],ml[i]);
                Bf_s_mq_mq[i] = myPV->Bf(Mz,s,mq[i],mq[i]);            
                B0_s_ml_ml[i] = myPV->B0(Mz,s,ml[i],ml[i]);
                B0_s_mq_mq[i] = myPV->B0(Mz,s,mq[i],mq[i]);            
            }
            delete myPV;
        }
        
        double ml2, vl2, al2, mq2, vq2, aq2;
        for (int i=0; i<6; i++) {
            ml2 = ml[i]*ml[i];
            vl2 = pow(EWSMC.vf((StandardModel::lepton) i), 2.0);
            al2 = pow(EWSMC.af((StandardModel::lepton) i), 2.0);            
            Sigma += - (vl2 + al2)*s*Bf_s_ml_ml[i] - 2.0*al2*ml2*B0_s_ml_ml[i];
            //
            mq2 = mq[i]*mq[i];
            vq2 = pow(EWSMC.vf((StandardModel::quark) i), 2.0);
            aq2 = pow(EWSMC.af((StandardModel::quark) i), 2.0);
            Sigma += 3.0*(- (vq2 + aq2)*s*Bf_s_mq_mq[i] - 2.0*aq2*mq2*B0_s_mq_mq[i]);
        
            //double Mw = EWSMC.GetMw();
            //std::cout << "l" << i << ": " 
            //          << - (vl2 + al2)*s*Bf_s_ml_ml[i]/Mw2 << " + "
            //          << - 2.0*al2*ml2*B0_s_ml_ml[i]/Mw2 << " = " 
            //          << (- (vl2 + al2)*s*Bf_s_ml_ml[i] - 2.0*al2*ml2*B0_s_ml_ml[i])/Mw2
            //          << std::endl;
            //std::cout << "q" << i << ": " 
            //          << 3.0*(- (vq2 + aq2)*s*Bf_s_mq_mq[i])/Mw2 << " + "
            //          << 3.0*(- 2.0*aq2*mq2*B0_s_mq_mq[i])/Mw2 << " = "
            //          << 3.0*(- (vq2 + aq2)*s*Bf_s_mq_mq[i] - 2.0*aq2*mq2*B0_s_mq_mq[i])/Mw2
            //          << std::endl;
        }
    }   
 
    // added O(alpha^2) contribution from the Z-gamma mixing 
    double Mw2 = pow(EWSMC.GetMw(), 2.0);
    Sigma += - Mw2*pow(PiZgamma_fer(Mz2), 2.0)*EWSMC.GetSM().getAle()/4.0/M_PI;
    
    return Sigma;
}




complex OneLoopEW::PiGammaGamma_fer(const double s) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = EWSMC.GetSM().getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = EWSMC.GetSM().getQuarks((StandardModel::quark) i).getMass();
    }
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_ml_ml[6], Bf_s_mq_mq[6];
    
    complex Pi(0.0,0.0,false);
    if (s==0.0) {
        throw "Codes for OneLoopEW::PiZgamma_fer(s=0.0) is missing";
    } else {
       if (s==Mz2) {
            for (int i=0; i<6; i++) {
                Bf_s_ml_ml[i] = EWSMC.GetBf_Mz2_ml_ml((StandardModel::lepton) i);
                Bf_s_mq_mq[i] = EWSMC.GetBf_Mz2_mq_mq((StandardModel::quark) i);           
            }
       } else {
            PVfunctions* myPV;
            myPV = new PVfunctions();
            for (int i=0; i<6; i++) {
                Bf_s_ml_ml[i] = myPV->Bf(Mz,s,ml[i],ml[i]);
                Bf_s_mq_mq[i] = myPV->Bf(Mz,s,mq[i],mq[i]);            
            }
            delete myPV;
       }
        
       double Ql, Qq;
       for (int i=0; i<6; i++) {
           Ql = EWSMC.Qf((StandardModel::lepton) i);
           Pi += 4.0*Ql*Ql*Bf_s_ml_ml[i];
           //
           Qq = EWSMC.Qf((StandardModel::quark) i);
           Pi += 4.0*3.0*Qq*Qq*Bf_s_mq_mq[i];
        }
    }   
    return Pi;
}


complex OneLoopEW::PiZgamma_fer(const double s) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = EWSMC.GetSM().getLeptons((StandardModel::lepton) i).getMass();
        mq[i] = EWSMC.GetSM().getQuarks((StandardModel::quark) i).getMass();
    }
    double Mz = EWSMC.GetSM().getMz();    
    double Mz2 = Mz*Mz;

    /* Loop functions */
    complex Bf_s_ml_ml[6], Bf_s_mq_mq[6];
    
    complex Pi(0.0,0.0,false);
    if (s==0.0) {
        throw "Codes for OneLoopEW::PiZgamma_fer(s=0.0) is missing";
    } else {
       if (s==Mz2) {
            for (int i=0; i<6; i++) {
                Bf_s_ml_ml[i] = EWSMC.GetBf_Mz2_ml_ml((StandardModel::lepton) i);
                Bf_s_mq_mq[i] = EWSMC.GetBf_Mz2_mq_mq((StandardModel::quark) i);           
            }
       } else {
            PVfunctions* myPV;
            myPV = new PVfunctions();
            for (int i=0; i<6; i++) {
                Bf_s_ml_ml[i] = myPV->Bf(Mz,s,ml[i],ml[i]);
                Bf_s_mq_mq[i] = myPV->Bf(Mz,s,mq[i],mq[i]);            
            }
            delete myPV;
       }
        
       double Ql, Qq;
       for (int i=0; i<6; i++) {
           Ql = EWSMC.Qf((StandardModel::lepton) i);
           Pi += (fabs(Ql) - 4.0*EWSMC.GetSW2()*Ql*Ql)*Bf_s_ml_ml[i];
           //
           Qq = EWSMC.Qf((StandardModel::quark) i);
           Pi += 3.0*(fabs(Qq) - 4.0*EWSMC.GetSW2()*Qq*Qq)*Bf_s_mq_mq[i];
        }
    }   
    return Pi;
}
 

double OneLoopEW::TEST_DeltaRhobar_bos() const {
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double rw = pow(EWSMC.GetSM().getMHl()/EWSMC.GetMw(), 2.0);
    double rz = pow(EWSMC.GetSM().getMHl()/EWSMC.GetSM().getMz(), 2.0);
    
    /* B0 functions for mu=Mw */
    complex B0_Mw_Mz2_Mw_Mw = EWSMC.GetB0_Mz2_Mw_Mw() + log(cW2);
    complex B0_Mw_Mz2_mh_Mz = EWSMC.GetB0_Mz2_mh_Mz() + log(cW2);
    complex B0_Mw_Mw2_Mz_Mw = EWSMC.GetB0_Mw2_Mz_Mw() + log(cW2);
    complex B0_Mw_Mw2_mh_Mw = EWSMC.GetB0_Mw2_mh_Mw() + log(cW2);
    
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


double OneLoopEW::DeltaRhobarW_bos() const {
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double cW4 = cW2*cW2;
    double rw = pow(EWSMC.GetSM().getMHl()/EWSMC.GetMw(), 2.0);
    
    complex B0_Mw2_Mz_Mw = EWSMC.GetB0_Mw2_Mz_Mw();
    complex B0_Mw2_mh_Mw = EWSMC.GetB0_Mw2_mh_Mw();
    
    double DRhobarW;    
    DRhobarW = - (1.0/12.0/cW4 + 4.0/3.0/cW2 - 17.0/3.0 - 4.0*cW2)*B0_Mw2_Mz_Mw.real()
               - (1.0 - 1.0/3.0*rw + 1.0/12.0*rw*rw)*B0_Mw2_mh_Mw.real()
               + (3.0/4.0/(1.0-rw) + 1.0/4.0 - 1.0/12.0*rw)*rw*log(rw)
               + (1.0/12.0/cW4 + 17.0/12.0/cW2 - 3.0/sW2 + 1.0/4.0)*log(cW2)
               + 1.0/12.0/cW4 + 11.0/8.0/cW2 + 139.0/36.0 - 177.0/24.0*cW2 
               + 5.0/8.0*cW4 - 1.0/12.0*rw*(7.0/2.0 - rw);     
    return DRhobarW;
} 

    
    
    