/* 
 * File:   EWSMOneLoopEW_HV.cpp
 * Author: mishima
 */

#include <cmath>
#include <stdexcept>
#include "EWSMOneLoopEW_HV.h"
#include "EWSMOneLoopEW.h"


EWSMOneLoopEW_HV::EWSMOneLoopEW_HV(const StandardModel& SM_i) : SM(SM_i) {
}


//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW_HV::SigmaWW_bos(const double mu, const double s,
                                      const double Mw) const {
    double Mz = SM.getMz(), mh = SM.getMHl();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz, mh2 = mh*mh;
    double cW2 = Mw2/Mz2, cW4 = cW2*cW2; 
    double sW2 = 1.0 - cW2, sW4 = sW2*sW2;
    double w = - s/Mw2, h = - s/mh2;
    double wh = mh2/Mw2, zh = mh2/Mz2;
        
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW_HV::SigmaWW_bos(s=0)");         
    } else {
        /* Loop functions */
        double A0_Mw = PV.A0(mu, Mw);
        double A0_mh = PV.A0(mu, mh);
        double A0_Mz = PV.A0(mu, Mz);
        complex B0_s_Mw_Mz = PV.B0(mu, s, Mw, Mz);
        complex B0_s_mh_Mw = PV.B0(mu, s, mh, Mw);
        complex B0_s_0_Mw = PV.B0(mu, s, 0.0, Mw);
        
        Sigma = Mw2/12.0
                *( - (sW4/cW4*(1.0 + 8.0*cW2)/w - 10.0/cW2 + 54.0 
                      + 16.0*cW2 + (1.0 - 40.0*cW2)*w)*B0_s_Mw_Mz
                   - ((1.0-wh)*(1.0-wh)/w + 2.0*wh - 10.0 + w)*B0_s_mh_Mw
                   - 8.0*sW2*(1.0/w + 2.0 - 5.0*w)*B0_s_0_Mw
                   + ((wh + 1.0/cW2 - 2.0)/w + 36.0/wh - 14.0)*A0_Mw/Mw2
                   - (1.0/h - 1.0/w - 7.0)*A0_mh/Mw2
                   - (sW2/cW2*(1.0 + 8.0*cW2)/w - 18.0/zh - 1.0 + 16.0*cW2)
                     *A0_Mz/Mw2
                   + 12.0/wh*(1.0/cW4 + 2.0) 
                   - 2.0*(1.0/cW2 + 18.0 + wh - 2.0*w/3.0) );    
    }
    return Sigma;
}


complex EWSMOneLoopEW_HV::SigmaWW_fer(const double mu, const double muForMq, 
                                      const double s) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = this->ml((StandardModel::lepton) i);
        mq[i] = this->mq((StandardModel::quark) i, muForMq); 
    }
    
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW_HV::SigmaWW_fer(s=0)");         
    } else {
        /* Loop functions */
        complex B1_s_ml_mlprime[3], B1_s_mq_mqprime[3]; 
        complex B1_s_mlprime_ml[3], B1_s_mqprime_mq[3];    
        complex Bf_s_mlprime_ml[3], Bf_s_mqprime_mq[3];
        for (int gen=0; gen<3; gen++) {
            B1_s_ml_mlprime[gen] = PV.B1(mu,s,ml[2*gen],ml[2*gen+1]);
            B1_s_mq_mqprime[gen] = PV.B1(mu,s,mq[2*gen],mq[2*gen+1]);
            B1_s_mlprime_ml[gen] = PV.B1(mu,s,ml[2*gen+1],ml[2*gen]);
            B1_s_mqprime_mq[gen] = PV.B1(mu,s,mq[2*gen+1],mq[2*gen]);
            Bf_s_mlprime_ml[gen] = PV.Bf(mu,s,ml[2*gen+1],ml[2*gen]);
            Bf_s_mqprime_mq[gen] = PV.Bf(mu,s,mq[2*gen+1],mq[2*gen]);            
        }
    
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
    }
    return Sigma;
}


complex EWSMOneLoopEW_HV::SigmaZZ_bos(const double mu, const double s,
                                      const double Mw) const {
    double Mz = SM.getMz(), mh = SM.getMHl();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz, mh2 = mh*mh;
    double cW2 = Mw2/Mz2, cW4 = cW2*cW2, cW6 = cW4*cW2;
    double z = - s/Mz2, h = - s/mh2;
    double wh = mh2/Mw2, zh = mh2/Mz2;
        
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW_HV::SigmaZZ_bos(s=0)");         
    } else {
        /* Loop functions */
        double A0_Mw = PV.A0(mu, Mw);
        double A0_mh = PV.A0(mu, mh);
        double A0_Mz = PV.A0(mu, Mz);
        complex B0_s_Mw_Mw = PV.B0(mu, s, Mw, Mw);
        complex B0_s_mh_Mz = PV.B0(mu, s, mh, Mz);
        
        Sigma = Mz2/12.0
                *( (4.0*cW2*(5.0 - 8.0*cW2 - 12.0*cW4) 
                     - (1.0 - 4.0*cW2 - 36.0*cW4)*z)*B0_s_Mw_Mw
                   - ((1.0-zh)*(1.0-zh)/z + 2.0*zh - 10.0 + z)*B0_s_mh_Mz
                   - (1.0/h - 1.0/z - 7.0)*A0_mh/Mz2
                   + 2.0*(18.0/wh + 1.0 + 8.0*cW2 - 24.0*cW4)*A0_Mw/Mz2
                   + (1.0/h - 1.0/z + 18.0/zh + 1.0)*A0_Mz/Mz2
                   + 2.0*( 6.0*(1.0 + 2.0*cW4)/zh - zh - 1.0 - 2.0*cW2 + 8.0*cW4
                           - 24.0*cW6 - 2.0/3.0*(1.0 - 2.0*cW2)*z ) );
    }
    return Sigma;
}


complex EWSMOneLoopEW_HV::SigmaZZ_fer(const double mu, const double muForMq, 
                                      const double s, const double Mw) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = this->ml((StandardModel::lepton) i);
        mq[i] = this->mq((StandardModel::quark) i, muForMq); 
    }
    
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        throw std::runtime_error("Missing codes for EWSMOneLoopEW_HV::SigmaZZ_fer(s=0)"); 
    } else {
        /* Loop functions */
        complex Bf_s_ml_ml[6], Bf_s_mq_mq[6];
        complex B0_s_ml_ml[6], B0_s_mq_mq[6];    
        for (int i=0; i<6; i++) {
            Bf_s_ml_ml[i] = PV.Bf(mu,s,ml[i],ml[i]);
            Bf_s_mq_mq[i] = PV.Bf(mu,s,mq[i],mq[i]);            
            B0_s_ml_ml[i] = PV.B0(mu,s,ml[i],ml[i]);
            B0_s_mq_mq[i] = PV.B0(mu,s,mq[i],mq[i]);            
        }

        double ml2, vl2, al2, mq2, vq2, aq2;
        for (int i=0; i<6; i++) {
            ml2 = ml[i]*ml[i];
            vl2 = pow(vl((StandardModel::lepton) i, Mw), 2.0);
            al2 = pow(al((StandardModel::lepton) i), 2.0);            
            Sigma += - (vl2 + al2)*s*Bf_s_ml_ml[i];
            Sigma += - 2.0*al2*ml2*B0_s_ml_ml[i];
            //
            mq2 = mq[i]*mq[i];
            vq2 = pow(vq((StandardModel::quark) i, Mw), 2.0);
            aq2 = pow(aq((StandardModel::quark) i), 2.0);
            Sigma += - 3.0*(vq2 + aq2)*s*Bf_s_mq_mq[i];
            Sigma += - 3.0*2.0*aq2*mq2*B0_s_mq_mq[i];
        }
    }   
    return Sigma;
}


complex EWSMOneLoopEW_HV::SigmaGammaGamma_bos(const double mu, const double s,
                                              const double Mw) const {
    return ( s*PiGammaGamma_bos(mu, s, Mw) );
}


complex EWSMOneLoopEW_HV::SigmaGammaGamma_fer(const double mu, const double muForMq,
                                              const double s) const {
    return ( s*PiGammaGamma_fer(mu, muForMq, s) );
}


complex EWSMOneLoopEW_HV::PiGammaGamma_bos(const double mu, const double s,
                                           const double Mw) const {
    double Mw2 = Mw*Mw;
    double w = - s/Mw2;
    
    complex Pi(0.0,0.0,false);
    if (s==0.0) {
        Pi = 3.0*log(Mw2/mu/mu) - 2.0/3.0;
    } else {
        /* Loop functions */
        double A0_Mw = PV.A0(mu, Mw);
        complex B0_s_Mw_Mw = PV.B0(mu, s, Mw, Mw);

        Pi = - Mw2/s*( (4.0 - 3.0*w)*B0_s_Mw_Mw + 4.0*(A0_Mw/Mw2 + 1.0) );
    }
    return Pi;
}


complex EWSMOneLoopEW_HV::PiGammaGamma_fer_l(const double mu, const double s, 
                                             const StandardModel::lepton l) const {
    // Neutrinos do not contribute, since Qf=0.
    if ( (l==StandardModel::NEUTRINO_1) || (l==StandardModel::NEUTRINO_2)
            || (l==StandardModel::NEUTRINO_3) )
        return 0.0;

    double mf = this->ml((StandardModel::lepton) l);
    double Qf = SM.getLeptons(l).getCharge();
    
    /* Loop functions */
    complex Bf_s_mf_mf;
    if (mf==0.0)
        Bf_s_mf_mf = 0.0; 
    else 
        Bf_s_mf_mf = PV.Bf(mu,s,mf,mf);
    
    return ( - 4.0*Qf*Qf*Bf_s_mf_mf );
}


complex EWSMOneLoopEW_HV::PiGammaGamma_fer_q(const double mu, const double muForMq,
                                             const double s, const StandardModel::quark q) const {
    double mf = this->mq((StandardModel::quark) q, muForMq);
    double Qf = SM.getQuarks(q).getCharge();
 
    /* Loop functions */
    complex Bf_s_mf_mf;
    if (mf==0.0)
        Bf_s_mf_mf = 0.0; 
    else 
        Bf_s_mf_mf = PV.Bf(mu,s,mf,mf);
    
    return ( - 4.0*3.0*Qf*Qf*Bf_s_mf_mf );
}


complex EWSMOneLoopEW_HV::PiGammaGamma_fer(const double mu, const double muForMq, 
                                           const double s) const {
    complex Pi(0.0,0.0,false);
    for (int i=0; i<6; i++) {
        Pi += PiGammaGamma_fer_l(mu, s, (StandardModel::lepton) i);
        Pi += PiGammaGamma_fer_q(mu, muForMq, s, (StandardModel::quark) i);        
    }
    return Pi;
}


complex EWSMOneLoopEW_HV::SigmaZgamma_bos(const double mu, const double s,
                                          const double Mw) const {
    double Mz = SM.getMz();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz;
    double cW2 = Mw2/Mz2;
    double w = - s/Mw2;
        
    complex Sigma(0.0,0.0,false);
    if (s==0.0) {
        Sigma = 2.0*Mw2*log(Mw2/mu/mu);
    } else {
        /* Loop functions */
        double A0_Mw = PV.A0(mu, Mw);
        complex B0_s_Mw_Mw = PV.B0(mu, s, Mw, Mw);
        
        Sigma = (4.0*(1.0/3.0 + cW2)/w - 1.0/6.0 - 3.0*cW2)*s*B0_s_Mw_Mw
                + (2.0/3.0 - 4.0*cW2)*Mw2*(A0_Mw/Mw2 + 1.0) - s/9.0;
    }
    return Sigma;
}


complex EWSMOneLoopEW_HV::SigmaZgamma_fer(const double mu, const double muForMq,
                                          const double s, const double Mw) const {
    double ml[6], mq[6];
    for (int i=0; i<6; i++) { 
        ml[i] = this->ml((StandardModel::lepton) i);
        mq[i] = this->mq((StandardModel::quark) i, muForMq); 
    }
    double Mz = SM.getMz();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz;
    double sW2 = 1.0 - Mw2/Mz2;
    
    /* Loop functions */
    complex Bf_s_ml_ml[6], Bf_s_mq_mq[6];
    for (int i=0; i<6; i++) {
        if (i==0 || i==2 || i==4 )
            Bf_s_ml_ml[i] = 0.0; // Neutrinos do not contribute, since Ql=0.
        else
            Bf_s_ml_ml[i] = PV.Bf(mu,s,ml[i],ml[i]);
        Bf_s_mq_mq[i] = PV.Bf(mu,s,mq[i],mq[i]);            
    }

    complex Pi(0.0,0.0,false);
    double Ql, Qq;
    for (int i=0; i<6; i++) {
        Ql = SM.getLeptons((StandardModel::lepton) i).getCharge();
        Pi += - (fabs(Ql) - 4.0*sW2*Ql*Ql)*Bf_s_ml_ml[i];
        //
        Qq = SM.getQuarks((StandardModel::quark) i).getCharge();
        Pi += - 3.0*(fabs(Qq) - 4.0*sW2*Qq*Qq)*Bf_s_mq_mq[i];
    }   
    return ( s*Pi );
}
 

//////////////////////////////////////////////////////////////////////// 

complex EWSMOneLoopEW_HV::F_Hollik(const double s, const double m1, 
                                   const double m2) const {
    double m12 = m1*m1, m22 = m2*m2;
    double mu = SM.getMz(); // The result is independent of mu. 
    
    if (m1!=0.0 && m2!=0.0) {
        if (m1==m2)
            return ( PV.B0(mu,s,m1,m1) + log(m1*m1/mu/mu) );
        else 
            return ( PV.B0(mu,s,m1,m2) + log(m1*m2/mu/mu) - 1.0 
                     + (m12 + m22)/(m12 - m22)*log(m1/m2) );
    } else if (m1==0.0 && m2!=0.0)  
        return ( PV.B0(mu,s,0.0,m2) + log(m2*m2/mu/mu) - 1.0 );    
    else if (m1!=0.0 && m2==0.0)  
        return ( PV.B0(mu,s,m1,0.0) + log(m1*m1/mu/mu) - 1.0 );

    else
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::F_Hollik()");
}


complex EWSMOneLoopEW_HV::Fprime_Hollik(const double muIR, const double s, 
                                        const double m1, const double m2) const {
    return ( PV.B0p(muIR,s,m1,m2) );
}


complex EWSMOneLoopEW_HV::SigmaWW_bos_Hollik(const double mu, const double s,
                                             const double Mw) const {
    double Mz = SM.getMz(), mh = SM.getMHl();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz, mh2 = mh*mh;
    double cW2 = Mw2/Mz2;
    double sW2 = 1.0 - cW2;
    double w = Mw2, z = Mz2, h = mh2;
        
    if (mu<=0.0 || s==0.0 || w==0.0 || z==0.0 || h==0.0 || h==w || z==w)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::SigmaWW_bos_Hollik()");

    double DeltaW = -log(Mw2/mu/mu);
    
    complex Sigma = - (19.0/2.0*s + 3.0*w*(1.0 - sW2/cW2))*DeltaW/3.0 // correct, Hollik (90)
                    //- (19.0/2.0*s + 3.0*w*(1.0 - sW2/cW2))*DeltaW // incorrect, Consoli, Hollik, Jegerlehner (89)
                    + (sW2*sW2*z - cW2/3.0*(7.0*(z + w) + 10.0*s - 2.0*(z - w)*(z - w)/s) 
                       - (w + z - s/2.0 - (z - w)*(z - w)/2.0/s)/6.0)
                      *F_Hollik(s, Mz, Mw)
                    + sW2/3.0*(-4.0*w - 10.0*s + 2.0*w*w/s)*F_Hollik(s, 0.0, Mw)
                    + (5.0*w - h + s/2.0 + (h - w)*(h - w)/2.0/s)/6.0
                      *F_Hollik(s, mh, Mw)
                    + (cW2/3.0*(7.0*z + 7.0*w + 10.0*s - 4.0*(z-w)) 
                    //   - sW2*sW2*z +(2.0*w - s/2.0)/6.0)*3.0*z/(z - w)*log(z/w) // Hollik (90)
                       - sW2*sW2*z +(2.0*w - s/2.0)/6.0)*z/(z - w)*log(z/w) // Consoli, Hollik, Jegerlehner (89)
                    - (2.0/3.0*w + s/12.0)*h/(h - w)*log(h/w) 
                    - cW2/3.0*(7.0*z + 7.0*w + 32.0/3.0*s) + sW2*sW2*z
                    + (5.0/3.0*s + 4.0*w - z - h)/6.0
                    - sW2/3.0*(4.0*w + 32.0/3.0*s);
    return Sigma;
}
    

complex EWSMOneLoopEW_HV::SigmaZZ_bos_Hollik(const double mu, const double s,
                                             const double Mw) const {
    double Mz = SM.getMz(), mh = SM.getMHl();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz, mh2 = mh*mh;
    double cW2 = Mw2/Mz2;
    double sW2 = 1.0 - cW2;
    double w = Mw2, z = Mz2, h = mh2;

    if (mu<=0.0 || s==0.0 || w==0.0 || z==0.0 || h==0.0 || h==z)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::SigmaZZ_bos_Hollik()");
    
    double DeltaW = -log(Mw2/mu/mu);
    
    complex Sigma = ( (3.0 - 19.0/6.0/sW2 + 1.0/6.0/cW2)*s 
                       + (4.0 + 1.0/cW2 - 1.0/sW2)*Mz2 )*DeltaW*sW2*cW2
                    + ( ( - cW2*cW2*(40.0*s + 80.0*w) 
                          + (cW2 - sW2)*(cW2 - sW2)*(8.0*w + s) + 12.0*w )
                          *F_Hollik(s, Mw, Mw)
                        + (10.0*z - 2.0*h + s + (h - z)*(h - z)/s)
                          *F_Hollik(s, mh, Mz)
                        - 2.0*h*log(h/w) - 2.0*z*log(z/w)
                        + (10.0*z - 2.0*h + s)
                          *(1.0 - (h + z)/(h - z)*log(mh/Mz) - log(mh*Mz/w))
                        + 2.0/3.0*s*(1.0 + (cW2 - sW2)*(cW2 - sW2) - 4.0*cW2)
                      )/12.0;
    return Sigma;
}


complex EWSMOneLoopEW_HV::SigmaGammaGamma_bos_Hollik(const double mu, const double s, 
                                                     const double Mw) const {
    double Mw2 = Mw*Mw;
    double w = Mw2;
        
    if (mu<=0.0)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::SigmaGammaGamma_bos_Hollik()");    

    double DeltaW = -log(Mw2/mu/mu);

    complex Sigma = - 3.0*s*DeltaW - (3.0*s + 4.0*w)*F_Hollik(s, Mw, Mw);
    return Sigma;
}


complex EWSMOneLoopEW_HV::PiGammaGamma_bos_Hollik(const double mu, const double s, 
                                                  const double Mw) const {
    double Mw2 = Mw*Mw;
    double w = Mw2;
    
    if (mu<=0.0)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::PiGammaGamma_bos_Hollik()");    

    double DeltaW = -log(Mw2/mu/mu);
    double muIR = Mw; // relevant only for Fprime_Hollik(muIR, Mw*Mw, 0, Mw)
    
    complex Pi;
    if (s==0.0) {
        //-- Pi(0) = dSigma(s)/ds|_{s=0} --
        Pi = - 3.0*DeltaW - 3.0*F_Hollik(s, Mw, Mw) 
             - (3.0*s + 4.0*w)*Fprime_Hollik(muIR, s, Mw, Mw);
    } else {
        //-- Pi(s) = Sigma(s)/s --
        Pi = SigmaGammaGamma_bos_Hollik(mu, s, Mw)/s;
    }
    return Pi;
}


complex EWSMOneLoopEW_HV::SigmaZgamma_bos_Hollik(const double mu, const double s,
                                                 const double Mw) const {
    double Mz = SM.getMz();
    double Mw2 = Mw*Mw, Mz2 = Mz*Mz;
    double cW2 = Mw2/Mz2;
    double w = Mw2;
        
    if (mu<=0.0)
        throw std::runtime_error("Missing cases in EWSMOneLoopEW_HV::SigmaZgamma_bos_Hollik()");    
    
    double DeltaW = -log(Mw2/mu/mu);

    complex Sigma = ( (3.0*cW2 + 1.0/6.0)*s + 2.0*w)*DeltaW
                    + ((3.0*cW2 + 1.0/6.0)*s + (4.0*cW2 + 4.0/3.0)*w)
                       *F_Hollik(s, Mw, Mw)
                    + s/9.0;
    return ( -Sigma ); // The minus sign is attributed to the different definition of s_W.
}




