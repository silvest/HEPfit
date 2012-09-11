/* 
 * File:   EWSMOneloopLEP2.cpp
 * Author: giovannigrilli
 */


#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf.h>
#include "EWSMOneLoopLEP2.h"
#include "EWSM.h"
#define EPSILON 0.00001

EWSMOneLoopLEP2::EWSMOneLoopLEP2(const EWSMcache& cache_i, const StandardModel& SM_i): 
                                 cache(cache_i), EWOL(cache_i), SM(SM_i){
}

///////////////////////////////////////////////////////////////////////////////

//should be ok but check!!
complex EWSMOneLoopLEP2::Chi_Z(const double mu, const double s, const double Mw_i,
                               const double W, const double X, const double Y) const {
    complex D_Z;
    double Mz = SM.getMz();
    double Mz2 = Mz*Mz;
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double epsilonZZ = cW2*W - 2.*sW*cW*X + sW2*Y;
    double gammaZ = 2.49465;
    complex i = complex::i();

//    complex SigmaZZ = SM.getAle()/4./sW2/cW2 * (EWOL.SigmaZZ_bos(mu,s,Mw) 
//            + EWOL.SigmaZZ_fer(mu,s,Mw));
//    complex SigmaGammaGamma = SM.getAle() / 4. * s * (EWOL.PiGammaGamma_bos(mu,s,Mw)
//            + EWOL.PiGammaGamma_fer(mu,s)  );
//    complex SigmaGammaZ = s * SM.getAle() / 4. /sW2/cW2 * ( EWOL.PiZgamma_bos(mu,s,Mw) 
//            + EWOL.PiZgamma_fer(mu,s,Mw));     
    
            
//    D_Z = 1./(s - Mz2 + SigmaZZ - SigmaGammaZ * SigmaGammaZ / (s + 
//           SigmaGammaGamma)-epsilonZZ/Mw/Mw);
    
    D_Z = 1./(s-Mz2+i*gammaZ*Mz);
      
    return(s*D_Z);
    
}

//should be ok but check!!
complex EWSMOneLoopLEP2::Chi_gamma(const double mu, const double s, const double Mw_i,
                               const double W, const double X, const double Y) const {
    complex D_gamma;
    
    double Mz = SM.getMz();
    double Mz2 = Mz*Mz;
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    
    double epsilongammagamma = sW2*W+2.*sW*cW*X+cW2*Y;
    
//    complex SigmaZZ = SM.getAle()/4./sW2/cW2 
//            * (EWOL.SigmaZZ_bos(mu,s,Mw) +
//             EWOL.SigmaZZ_fer(mu,s,Mw))
//            ;
//    complex SigmaGammaGamma = SM.getAle()/4.*s
//            * (EWOL.PiGammaGamma_bos(mu,s,Mw) +
//            EWOL.PiGammaGamma_fer(mu,s)  )
//            ;
//    complex SigmaGammaZ = s * SM.getAle()/4./sW2/cW2 
//            * ( EWOL.PiZgamma_bos(mu,s,Mw) +
//             EWOL.PiZgamma_fer(mu,s,Mw))
//            ;
    
            
//    D_gamma = 1./(s +  SigmaGammaGamma - SigmaGammaZ * SigmaGammaZ / (s - Mz2 + 
//           SigmaZZ))-epsilongammagamma/Mw/Mw;
    
    D_gamma = 1./s;
      
        return(s*D_gamma);
    
}

//should be ok but check!!
complex EWSMOneLoopLEP2::Chi_gammaZ(const double mu, const double s, const double Mw_i,
                               const double W, const double X, const double Y) const {
    complex D_Z;
    complex D_gammaZ;
    double Mz = SM.getMz();
    double Mz2 = Mz*Mz;
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double epsilongammaZ = (cW2-sW2)*X+sW*cW*(W-Y);

    complex SigmaZZ = SM.getAle() / 4./sW2/cW2 * (EWOL.SigmaZZ_bos(mu,s,Mw) +
            EWOL.SigmaZZ_fer(mu,s,Mw));
    complex SigmaGammaGamma = SM.getAle() / 4. * s * (EWOL.PiGammaGamma_bos(mu,s,Mw) + 
            EWOL.PiGammaGamma_fer(mu,s)  );
    complex SigmaGammaZ = s * SM.getAle() / 4./sW2/cW2 * ( EWOL.PiZgamma_bos(mu,s,Mw) +
            EWOL.PiZgamma_fer(mu,s,Mw));
    
            
    D_Z = 1./(s - Mz2 + SigmaZZ - SigmaGammaZ * SigmaGammaZ / (s + 
           SigmaGammaGamma));
    
    D_gammaZ = - SigmaGammaZ / (s + SigmaGammaGamma) * D_Z-epsilongammaZ/Mw/Mw;
      
        return(-s*D_gammaZ);
    
}

//finish and check!!!
double EWSMOneLoopLEP2::g_rhofq(const QCD::quark q, const double rho,const double Mw_i) const{
    double Mw = SM.Mw();
    
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3 = SM.getQuarks(q).getIsospin();
    double Q = SM.getQuarks(q).getMass();
   
    if(rho == 0.5) {
        return (-Q*sW/cW);
    } else if (rho == -0.5) {
        return ((I3 - sW2*Q) / (sW*cW));   
    } else {
        throw "Error in EWSMOneLoopLEP2::g_rho(): rho must be 1./2 or -1./2";
    }
}


double EWSMOneLoopLEP2::g_rhofl(const StandardModel::lepton l, const double rho,const double Mw_i) const {
    double Mw = SM.Mw();
    
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3 = SM.getLeptons(l).getIsospin();
    double Q = SM.getLeptons(l).getCharge();
    double g;
   
    if(rho == 0.5) {
       g = (-Q*sW/cW);
    } else if (rho == -0.5) {
       g = ((I3 - sW2*Q) / (sW*cW));   
    } else {
        throw "Error in EWSMOneLoopLEP2::g_rho(): rho must be 1./2 or -1./2";
    }
    
    return g;
    
}

//check!!
double EWSMOneLoopLEP2::g_rhoe(const double rho, const double Mw_i) const {
    double Mw = SM.Mw();
    
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double g;
   
    if(rho == 0.5) {
        g = (sW/cW);
    } else if (rho == -0.5) {
        g = ((-1./2. + sW2) / (sW*cW));   
    } else {
        throw "Error in EWSMOneLoopLEP2::g_rho(): rho must be 1./2 or -1./2";
    }
    
    return g;
}

double EWSMOneLoopLEP2::t(const double mf, const double s, const double cos_theta) const{    
    
    //return (-s/2.*(1.-sqrt(1.-4.*mf*mf/s)*cos_theta-2.*mf*mf/s));
    return -s/2.*(1.-cos_theta);
    
}

double EWSMOneLoopLEP2::u(const double mf, const double s, const double cos_theta) const{
    
    //return (-s/2.*(1.+sqrt(1.-4.*mf*mf/s)*cos_theta-2.*mf*mf/s));
     return -s/2.*(1.+cos_theta);
}


//finish and check!!
complex EWSMOneLoopLEP2::Al(const double mu,const StandardModel::lepton l, 
                            const double rho, const double k, const double s, const double Mw_i,
                            const double W, const double X, const double Y) const{
    double Qf = SM.getLeptons(l).getCharge();
    
    return (-Qf*Chi_gamma(mu,s,Mw_i,W,X,Y) + g_rhoe(k,Mw_i)*g_rhofl(l,rho,Mw_i)*Chi_Z(mu,s,Mw_i,W,X,Y));
    
}

///Mw_i
complex EWSMOneLoopLEP2::Aq(const double mu,const QCD::quark q,
                            const double rho, const double k, const double s, const double Mw_i,
                            const double W, const double X, const double Y) const{
    double Qf = SM.getQuarks(q).getCharge();
    
    return (-Qf*Chi_gamma(mu,s,Mw_i,W,X,Y) + g_rhoe(k,Mw_i)*g_rhofq(q,rho,Mw_i)*Chi_Z(mu,s,Mw_i,W,X,Y));
    
}

//finish and check!!
complex EWSMOneLoopLEP2::Bl(const double mu,StandardModel::lepton l, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const{
    double Qf = SM.getLeptons(l).getCharge();///////////////////////////////////////////
    
    return ((-g_rhofl(l,rho,Mw_i) + Qf*g_rhoe(rho,Mw_i))*Chi_gammaZ(mu,s,Mw_i,W,X,Y));
    
}

//finish and check!!
complex EWSMOneLoopLEP2::Bq(const double mu,QCD::quark q, const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const{
    double Qf = SM.getQuarks(q).getCharge();///////////////////////////////////////////
    
    return ((-g_rhofq(q,rho,Mw_i) + Qf*g_rhoe(rho,Mw_i))*Chi_gammaZ(mu,s,Mw_i,W,X,Y));
    
}


//should be ok
complex EWSMOneLoopLEP2::Lambda2(const double m, const double s)const{
    
    complex i = complex::i();
    complex w = m*m/(s+i*EPSILON);
    //double x = (sqrt(1. - 4. * w) - 1.)/ (sqrt(1. - 4. * w) + 1.);
    
    complex y = 1 + 1./w;
    complex Li2;
    gsl_sf_result re, im;
    gsl_sf_complex_dilog_xy_e(y.real(), y.imag(), &re, &im);
    Li2.real() = re.val;
    Li2.imag() = im.val;
    
    return (-3.5-2.*w-(2.*w+3.)*(log(-w)+2.*(1+w)*(1+w)*(Li2 - M_PI * M_PI / 6.)));
    
}

//should be ok
complex EWSMOneLoopLEP2::Lambda3(const double m, const double s) const{
    
    complex i = complex::i();
    complex w = m*m/(s+i*EPSILON);
    complex x = (sqrt(1. - 4. * w) - 1.)/ (sqrt(1. - 4. * w) + 1.);
    
    return (5./6.-2.*w/3.-(2.*w + 1) / 3. * sqrt(1-4.*w) * log(x) + 
            2. / 3. * w*(w+2.) * log(x)*log(x));
}



complex EWSMOneLoopLEP2::FLgammal(const double s, const double Mw_i) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    
    return (3./4./sW2*Lambda3(Mw,s));
    
}


complex EWSMOneLoopLEP2::FLZl(const double s, const double Mw_i) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    
    return (1./8./sW2/sW/cW*Lambda2(Mw,s) - 3.*cW/4./sW2/sW*Lambda3(Mw,s));
    
}



complex EWSMOneLoopLEP2::Cl(const double mu,const StandardModel::lepton l, 
        const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const{
    double Mz = SM.getMz();
    double Qf = SM.getLeptons(l).getCharge();
    double alpha = SM.getAle();
    
    if(k == -1./2.){
        return (-alpha/4./M_PI * (g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)*Lambda2(Mz,s)+2.*FLgammal(s, Mw_i))*Qf*Chi_Z(mu,s,Mw_i,W,X,Y)/s);
    } else if(k == 1./2.){
        return (-alpha/4./M_PI * (g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)*Lambda2(Mz,s))*Qf*Chi_Z(mu,s,Mw_i,W,X,Y)/s);
    } else {
        throw "Error in EWSMOneLoopLEP2::Cl(): k must be 1./2 or -1./2";
    }
    
    
    
}

complex EWSMOneLoopLEP2::Cq(const double mu,const QCD::quark q, 
        const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const{
    double Mz = SM.getMz();
    double Qf = SM.getQuarks(q).getCharge();
    double alpha = SM.getAle();
    
    if(k == -1./2.){
        return (-alpha/4./M_PI * (g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)*Lambda2(Mz,s)+2.*FLgammal(s, Mw_i))*Qf*Chi_Z(mu,s,Mw_i,W,X,Y)/s);
    } else if(k == 1./2.){
        return (-alpha/4./M_PI * (g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)*Lambda2(Mz,s))*Qf*Chi_Z(mu,s,Mw_i,W,X,Y)/s);
    } else {
        throw "Error in EWSMOneLoopLEP2::Cq(): k must be 1./2 or -1./2";
    }
    
}


complex EWSMOneLoopLEP2::Dl_rho(const double mu,const StandardModel::lepton l, 
        const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const{
    double Mz = SM.getMz();
    double alpha = SM.getAle();
    complex x;
    
    if(k == -1./2.){
        x = alpha/4./M_PI * (g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)
                *Lambda2(Mz,s)+2.*FLZl(s, Mw_i))*Chi_Z(mu,s,Mw_i,W,X,Y)/s*g_rhofl(l,rho,Mw_i);
    } else if(k == 1./2.){
        x = alpha/4./M_PI * (g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)
                *Lambda2(Mz,s))*Chi_Z(mu,s,Mw_i,W,X,Y)/s*g_rhofl(l,rho,Mw_i);
    } else {
        throw "Error in EWSMOneLoopLEP2::Dl_rho(): k must be 1./2 or -1./2";
    }
    
    return x;
    
}


complex EWSMOneLoopLEP2::Dq_rho(const double mu,const QCD::quark q, 
        const double rho, const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const{
    double Mz = SM.getMz();
    
    double alpha = SM.getAle();
    
    if(k == -1./2.){
        return (alpha/4./M_PI * (g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)
                *Lambda2(Mz,s)+2.*FLZl(s, Mw_i))*Chi_Z(mu,s,Mw_i,W,X,Y)/s*g_rhofq(q,rho,Mw_i));
    } else if(k == 1./2.){
        return (alpha/4./M_PI * (g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)*g_rhoe(rho,Mw_i)
                *Lambda2(Mz,s))*Chi_Z(mu,s,Mw_i,W,X,Y)/s*g_rhofq(q,rho,Mw_i));
    } else {
        throw "Error in EWSMOneLoopLEP2::Dq_rho(): k must be 1./2 or -1./2";
    }
    
    
    
}


double EWSMOneLoopLEP2::ve(const double Mw_i) const{
    
    return (0.5*(g_rhoe(0.5,Mw_i)+g_rhoe(-0.5,Mw_i)));
    
}

double EWSMOneLoopLEP2::ae(const double Mw_i) const{
    
    return (0.5*(g_rhoe(-0.5,Mw_i)-g_rhoe(0.5,Mw_i)));
    
}

double EWSMOneLoopLEP2::vq(const QCD::quark q, const double Mw_i) const{
    
    return (0.5*(g_rhofq(q,0.5,Mw_i)+g_rhofq(q,-0.5,Mw_i)));
    
}

double EWSMOneLoopLEP2::aq(const QCD::quark q, const double Mw_i) const{
    
    return (0.5*(g_rhofq(q,-0.5,Mw_i)-g_rhofq(q,0.5,Mw_i)));
    
}

double EWSMOneLoopLEP2::vl(const StandardModel::lepton l, const double Mw_i) const{
    
    return (0.5*(g_rhofl(l,0.5,Mw_i)+g_rhofl(l,-0.5,Mw_i)));
    
}

double EWSMOneLoopLEP2::al(const StandardModel::lepton l, const double Mw_i) const{
    
    return (0.5*(g_rhofl(l,-0.5,Mw_i)-g_rhofl(l,0.5,Mw_i)));
    
}


complex EWSMOneLoopLEP2::C1plus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const{
   
    //complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getLeptons(l).getMass();
    double mf2 = mf*mf;
    
    return (1./(4.*mf2-s)*(PV.B0(mu,s,m1,m2)-0.5*(PV.B0(mu,mf2,m3,m1)
            +PV.B0(mu,mf2,m3,m2))+0.5*(2.*m3*m3+2.*mf2-m1*m1-m2*m2)));
  
}

complex EWSMOneLoopLEP2::C1plus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const{
   
    //complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getQuarks(q).getMass();
    double mf2 = mf*mf;
    
    return (1./(4.*mf2-s)*(PV.B0(mu,s,m1,m2)-0.5*(PV.B0(mu,mf2,m3,m1)
            +PV.B0(mu,mf2,m3,m2))+0.5*(2.*m3*m3+2.*mf2-m1*m1-m2*m2)));
  
}
complex EWSMOneLoopLEP2::C1minus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const{
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getLeptons(l).getMass();
    double mf2 = mf*mf;
    
    return (1./2./s*(PV.B0(mu,mf2,m3,m1)-PV.B0(mu,mf2,m3,m2) + (m2*m2-m1*m1)*C0));
  
}

complex EWSMOneLoopLEP2::C1minus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const {
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getQuarks(q).getMass();
    double mf2 = mf*mf;
    
    return (1./2./s*(PV.B0(mu,mf2,m3,m1)-PV.B0(mu,mf2,m3,m2) + (m2*m2-m1*m1)*C0));
  
}

complex EWSMOneLoopLEP2::C20_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const{
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getLeptons(l).getMass();
    double mf2 = mf*mf;
    
    return (1./2.*m3*m3*C0+0.25+0.25*((m1*m1+m2*m2-2.*m3*m3-2*mf2)*C1plus_l(mu,l,m1,m2,m3,s)
             +(m1*m1-m2*m2)*C1minus_l(mu,l,m1,m2,m3,s))+PV.B0(mu,s,m1,m2));
  
}

complex EWSMOneLoopLEP2::C20_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const {
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getQuarks(q).getMass();
    double mf2 = mf*mf;
    
    return (1./2.*m3*m3*C0+0.25+0.25*((m1*m1+m2*m2-2.*m3*m3-2*mf2)*C1plus_q(mu,q,m1,m2,m3,s)
             +(m1*m1-m2*m2)*C1minus_q(mu,q,m1,m2,m3,s))+PV.B0(mu,s,m1,m2));
  
}

complex EWSMOneLoopLEP2::C2plus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const {
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getLeptons(l).getMass();
    double mf2 = mf*mf;
    
    return (1./(4.*mf2-s)*(-C20_l(mu,l,m1,m2,m3,s)+0.5*(2.*m3*m3+2.*mf2-m1*m1-m2*m2)*C1plus_l(mu,l,m1,m2,m3,s)
            +0.5*PV.B0(mu,s,m1,m2)+0.25*(PV.B1(mu,mf2,m3,m1)+PV.B1(mu,mf2,m3,m2))));
  
}

complex EWSMOneLoopLEP2::C2plus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const{
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getQuarks(q).getMass();
    double mf2 = mf*mf;
    
    return (1./(4.*mf2-s)*(-C20_q(mu,q,m1,m2,m3,s)+0.5*(2.*m3*m3+2.*mf2-m1*m1-m2*m2)*C1plus_q(mu,q,m1,m2,m3,s)
            +0.5*PV.B0(mu,s,m1,m2)+0.25*(PV.B1(mu,mf2,m3,m1)+PV.B1(mu,mf2,m3,m2))));
  
}


complex EWSMOneLoopLEP2::C2minus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const{
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getLeptons(l).getMass();
    double mf2 = mf*mf;
    
    return (1./s*(-C20_l(mu,l,m1,m2,m3,s)+0.5*(m2*m2-m1*m1)*C1minus_l(mu,l,m1,m2,m3,s)
            -0.25*(PV.B1(mu,mf2,m3,m1)+PV.B1(mu,mf2,m3,m2))));
  
}

complex EWSMOneLoopLEP2::C2minus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const {
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getQuarks(q).getMass();
    double mf2 = mf*mf;
    
    return (1./s*(-C20_q(mu,q,m1,m2,m3,s)+0.5*(m2*m2-m1*m1)*C1minus_q(mu,q,m1,m2,m3,s)
            -0.25*(PV.B1(mu,mf2,m3,m1)+PV.B1(mu,mf2,m3,m2))));
  
}

complex EWSMOneLoopLEP2::C2plusminus_l(const double mu, const StandardModel::lepton l,
                    const double m1, const double m2, const double m3,
                    const double s) const{
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getLeptons(l).getMass();
    double mf2 = mf*mf;
    
    return (1./4./s*(2.*(m2*m2-m1*m1)*C1plus_l(mu,l,m1,m2,m3,s)
            +PV.B1(mu,mf2,m3,m2)-PV.B1(mu,mf2,m3,m1)));
  
}

complex EWSMOneLoopLEP2::C2plusminus_q(const double mu, const QCD::quark q,
                    const double m1, const double m2, const double m3,
                    const double s) const{
   
    complex C0 = -PV.C0(s,m1,m3,m2);
    double mf = SM.getQuarks(q).getMass();
    double mf2 = mf*mf;
    
    return (1./4./s*(2.*(m2*m2-m1*m1)*C1plus_q(mu,q,m1,m2,m3,s)
            +PV.B1(mu,mf2,m3,m2)-PV.B1(mu,mf2,m3,m1)));
  
}



double EWSMOneLoopLEP2::muf(const double Mw_i, const double mf) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    //double mf = SM.getQuarks(q).getMass();
    
    return (mf / (2.*cW*sW*Mz));
    
}

//double EWSMOneLoopLEP2::muf_l(const double Mw_i, const StandardModel::lepton l) const{
//    double Mw = SM.Mw();
//    double Mz = SM.getMz();
//    double sW2 = SM.sW2();
//    double cW2 = SM.cW2();
//    double sW = sqrt(sW2);
//    double cW = sqrt(cW2);
//    double mf = SM.getLeptons(l).getMass();
//    
//    return (mf / (2.*cW*sW*Mz));
    
//}


double EWSMOneLoopLEP2::alphaf_q(const double Mw_i, const QCD::quark q) const{
    double Mw = SM.Mw();
    double sW2 =  SM.sW2();
    double cW2 = SM.cW2();
    
    return (aq(q,Mw_i)*(sW2-cW2));
    
}

double EWSMOneLoopLEP2::alphaf_qprime(const double Mw_i, const QCD::quark q) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    QCD::quark qprime;
    
    if(q==SM.UP){
      qprime = SM.DOWN;      /////
    } else if(q == SM.DOWN){
      qprime = SM.UP;          
    } else if(q == SM.CHARM){
      qprime = SM.STRANGE;   
    } else if(q == SM.STRANGE){
      qprime = SM.CHARM;    
    } else if(q == SM.TOP){
      qprime = SM.BOTTOM;    
    } else if(q == SM.BOTTOM){
      qprime = SM.TOP;    
    }
    
    return (aq(qprime,Mw_i)*(sW2-cW2));
    
}


QCD::quark EWSMOneLoopLEP2::qprime(const QCD::quark q) const{
   
    QCD::quark qprime_t;
    
    if(q==SM.UP){
      qprime_t = SM.DOWN;      /////
    } else if(q == SM.DOWN){
      qprime_t = SM.UP;          
    } else if(q == SM.CHARM){
      qprime_t = SM.STRANGE;   
    } else if(q == SM.STRANGE){
      qprime_t = SM.CHARM;    
    } else if(q == SM.TOP){
      qprime_t = SM.BOTTOM;    
    } else if(q == SM.BOTTOM){
      qprime_t = SM.TOP;    
    }
    
    return (qprime_t);
    
}

double EWSMOneLoopLEP2::alphaf_l(const double Mw_i, const StandardModel::lepton l) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    
    return (al(l,Mw_i)*(sW2-cW2));
    
}


double EWSMOneLoopLEP2::C1f_l(const StandardModel::lepton l, const double Mw_i,
                             const double mu) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mf=SM.getLeptons(l).getMass();
    double mf2 = mf*mf;
    //Because we use only e^+e^-\to mu^+mu^- the mass mfprime = mnu = 0
    double mfprime = 0.;     
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    
    double vl2 = vl(l,Mw_i)*vl(l,Mw_i);
    double al2 = al(l,Mw_i)*al(l,Mw_i);
    double muf2 = muf(Mw_i,mf)*muf(Mw_i,mf);
    double mufprime2 = muf(Mw_i,mfprime)*muf(Mw_i,mfprime);

    complex C1f_tmp = SM.getAle() /4./M_PI * ((vl2+al2)*(2.*PV.B1(mu,mf2,mf,Mz)+1.)+
            1./4./sW2*(2.*PV.B1(mu,mf2,mfprime,Mw)+1) + muf2*
            (PV.B1(mu,mf2,mf,MH) +PV.B1(mu,mf2,mf,Mz)+PV.B1(mu,mf2,mfprime,Mw))+
            mufprime2*PV.B1(mu,mf2,mfprime,Mw)+2.*mf2*(vl2*(2.*PV.B1p(mu,mf2,mf,Mz)+ 
            4.*PV.B0p(mu,mf2,mf,Mz))+al2*(2.*PV.B1p(mu,mf2,mf,Mz)- 
            4.*PV.B0p(mu,mf2,mf,Mz)) + 1./2./sW2*PV.B1p(mu,mf2,mfprime,Mw)+
            muf2*(PV.B1p(mu,mf2,mf,MH)+PV.B1p(mu,mf2,mf,Mz)
            +PV.B1p(mu,mf2,mfprime,Mw)+PV.B0p(mu,mf2,mf,Mz)-PV.B0p(mu,mf2,mf,MH))
            +mufprime2*(PV.B1p(mu,mf2,mfprime,Mw)+ 2.*PV.B0p(mu,mf2,mfprime,Mw))));
    
    return (C1f_tmp.real());
    
}

//double EWSMOneLoopLEP2::mqprime(const QCD::quark q) const{
//    
//    double mfprime;
//    
//    if(q==0){
//      mfprime = SM.getQuarks(1).getMass();     /////
//    } else if(q == 1){
//      mfprime = SM.getQuarks(0).getMass();          
//    } else if(q == 2){
//      mfprime = SM.getQuarks(3).getMass();   
//    } else if(q == 3){
//      mfprime = SM.getQuarks(2).getMass();   
//    } else if(q == 4){
//      mfprime = SM.getQuarks(5).getMass();   
//    } else if(q == 5){
//      mfprime = SM.getQuarks(4).getMass();   
//    }
//    
//    return mfprime;
//    
//}

double EWSMOneLoopLEP2::Qqprime(const QCD::quark q) const{
    
    double Qfprime;
    
    if(q==0){
      Qfprime = SM.getQuarks(1).getCharge();     /////
    } else if(q == 1){
      Qfprime = SM.getQuarks(0).getCharge();          
    } else if(q == 2){
      Qfprime = SM.getQuarks(3).getCharge();   
    } else if(q == 3){
      Qfprime = SM.getQuarks(2).getCharge();   
    } else if(q == 4){
      Qfprime = SM.getQuarks(5).getCharge();   
    } else if(q == 5){
      Qfprime = SM.getQuarks(4).getCharge();   
    }
    
    return Qfprime;
    
}

double EWSMOneLoopLEP2::I3qprime(const QCD::quark q) const{
    
    double I3fprime;
    
    if(q==0){
      I3fprime = SM.getQuarks(1).getIsospin();     /////
    } else if(q == 1){
      I3fprime = SM.getQuarks(0).getIsospin();          
    } else if(q == 2){
      I3fprime = SM.getQuarks(3).getIsospin();   
    } else if(q == 3){
      I3fprime = SM.getQuarks(2).getIsospin();   
    } else if(q == 4){
      I3fprime = SM.getQuarks(5).getIsospin();   
    } else if(q == 5){
      I3fprime = SM.getQuarks(4).getIsospin();   
    }
    
    return I3fprime;
    
}

double EWSMOneLoopLEP2::C1f_q(const QCD::quark q, const double Mw_i,
                             const double mu) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    
    double mf=SM.getQuarks(q).getMass();
    double mf2 = mf*mf;
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    
    double vq2 = vq(q,Mw_i)*vq(q,Mw_i);
    double aq2 = aq(q,Mw_i)*aq(q,Mw_i);
    double muf2 = muf(Mw_i,mf)*muf(Mw_i,mf);
    double muqprime2 = muf(Mw_i,mfprime)*muf(Mw_i,mfprime);

    complex C1f_tmp = SM.getAle() /4./M_PI * ((vq2+aq2)*(2.*PV.B1(mu,mf2,mf,Mz)+1.)+
            1./4./sW2*(2.*PV.B1(mu,mf2,mfprime,Mw)+1) + muf2*
            (PV.B1(mu,mf2,mf,MH) +PV.B1(mu,mf2,mf,Mz)+PV.B1(mu,mf2,mfprime,Mw))+
            muqprime2*PV.B1(mu,mf2,mfprime,Mw)+2.*mf2*(vq2*(2.*PV.B1p(mu,mf2,mf,Mz)+ 
            4.*PV.B0p(mu,mf2,mf,Mz))+aq2*(2.*PV.B1p(mu,mf2,mf,Mz)- 
            4.*PV.B0p(mu,mf2,mf,Mz)) + 1./2./sW2*PV.B1p(mu,mf2,mfprime,Mw)+
            muf2*(PV.B1p(mu,mf2,mf,MH)+PV.B1p(mu,mf2,mf,Mz)
            +PV.B1p(mu,mf2,mfprime,Mw)+PV.B0p(mu,mf2,mf,Mz)-PV.B0p(mu,mf2,mf,MH))
            +muqprime2*(PV.B1p(mu,mf2,mfprime,Mw)+ 2.*PV.B0p(mu,mf2,mfprime,Mw))));
    
    return (C1f_tmp.real());
    
}


double EWSMOneLoopLEP2::C2f_l(const StandardModel::lepton l, const double Mw_i,
                             const double mu) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mf=SM.getLeptons(l).getMass();
    double mf2 = mf*mf;
    //Because we use only e^+e^-\to mu^+mu^- the mass mfprime = mnu = 0
    double mfprime = 0.;     
    double Mz = SM.getMz();
    
    double muf2 = muf(Mw_i,mf)*muf(Mw_i,mf);
    double mufprime2 = muf(Mw_i,mfprime)*muf(Mw_i,mfprime);

    complex C2f_tmp = -SM.getAle() /4./M_PI * (-2.*vl(l,Mw_i)*al(l,Mw_i)
            *(2.*PV.B1(mu,mf2,mf,Mz)+1.)-1./4./sW2*(2.*PV.B1(mu,mf2,mfprime,Mw)+1.)+
            (muf2-mufprime2)*PV.B1(mu,mf2,mfprime,Mw));
    
    return (C2f_tmp.real());
    
}

double EWSMOneLoopLEP2::C2f_q(const QCD::quark q, const double Mw_i,
                             const double mu) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    
    double mf=SM.getQuarks(q).getMass();
    double mf2 = mf*mf;
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    
    double Mz = SM.getMz();
    
    double muf2 = muf(Mw_i,mf)*muf(Mw_i,mf);
    double muqprime2 = muf(Mw_i,mfprime)*muf(Mw_i,mfprime);

    complex C2f_tmp = -SM.getAle() /4./M_PI * (-2.*vq(q,Mw_i)*aq(q,Mw_i)
            *(2.*PV.B1(mu,mf2,mf,Mz)+1.)-1./4./sW2*(2.*PV.B1(mu,mf2,mfprime,Mw)+1.)+
            (muf2-muqprime2)*PV.B1(mu,mf2,mfprime,Mw));
    
    return (C2f_tmp.real());
    
}

double EWSMOneLoopLEP2::CVgammal(const StandardModel::lepton l,const double Mw_i, const double mu) const{
    double Qf = SM.getLeptons(l).getCharge();
    
    return (-Qf*C1f_l(l,Mw_i,mu));
    
}

double EWSMOneLoopLEP2::CVgammaq(const QCD::quark q,const double Mw_i, const double mu) const{
    double Qf = SM.getQuarks(q).getCharge();
    
    return (-Qf*C1f_q(q,Mw_i,mu));
    
}

double EWSMOneLoopLEP2::CAgammal(const StandardModel::lepton l,const double Mw_i, const double mu) const{
    double Qf = SM.getLeptons(l).getCharge();
    
    return (-Qf*C2f_l(l,Mw_i,mu));
    
}

double EWSMOneLoopLEP2::CAgammaq(const QCD::quark q,const double Mw_i, const double mu) const{
    double Qf = SM.getQuarks(q).getCharge();
    
    return (-Qf*C2f_q(q,Mw_i,mu));
    
}


double EWSMOneLoopLEP2::CVZl(const StandardModel::lepton l,const double Mw_i, const double mu) const{
    
    return (vl(l,Mw_i)*C1f_l(l,Mw_i,mu)+al(l,Mw_i)*C2f_l(l,Mw_i,mu));
    
}

double EWSMOneLoopLEP2::CVZq(const QCD::quark q,const double Mw_i, const double mu) const{
    
    return (vq(q,Mw_i)*C1f_q(q,Mw_i,mu)+aq(q,Mw_i)*C2f_q(q,Mw_i,mu));
    
}

double EWSMOneLoopLEP2::CAZl(const StandardModel::lepton l,const double Mw_i, const double mu) const{
    
    return (al(l,Mw_i)*C1f_l(l,Mw_i,mu)+vl(l,Mw_i)*C2f_l(l,Mw_i,mu));
    
}

double EWSMOneLoopLEP2::CAZq(const QCD::quark q,const double Mw_i, const double mu) const{
    
    return (aq(q,Mw_i)*C1f_q(q,Mw_i,mu)+vq(q,Mw_i)*C2f_q(q,Mw_i,mu));
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////   VECTOR FORM FACTOR      /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////


complex EWSMOneLoopLEP2::FI_Val(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double Qf = SM.getLeptons(l).getCharge();
    double vl2 = vl(l,Mw_i)*vl(l,Mw_i);
    double al2 = al(l,Mw_i)*al(l,Mw_i);
    double lVp = -Qf*(vl2+al2);
    double mf=SM.getLeptons(l).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(4.*C20_l(mu,l,mf,mf,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_l(mu,l,mf,mf,Mz,s)+2.*s*C2minus_l(mu,l,mf,mf,Mz,s)
            -4.*(4.*mf*mf-s)*C1plus_l(mu,l,mf,mf,Mz,s)+(6.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lVp*C0) );
    
}//ok

complex EWSMOneLoopLEP2::FI_Vaq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double Qf = SM.getQuarks(q).getCharge();
    double vq2 = vq(q,Mw_i)*vq(q,Mw_i);
    double aq2 = aq(q,Mw_i)*aq(q,Mw_i);
    double lVp = -Qf*(vq2+aq2);
    double mf=SM.getQuarks(q).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(4.*C20_q(mu,q,mf,mf,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_q(mu,q,mf,mf,Mz,s)+2.*s*C2minus_q(mu,q,mf,mf,Mz,s)
            -4.*(4.*mf*mf-s)*C1plus_q(mu,q,mf,mf,Mz,s)+(6.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lVp*C0) );
    
}//ok


complex EWSMOneLoopLEP2::FI_Vbl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double lVp = 0.;
    double mf=SM.getLeptons(l).getMass();
    double mfprime =0.;
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    
    return (SM.getAle()/4./M_PI *(lVp*(4.*C20_l(mu,l,mfprime,mfprime,Mw,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_l(mu,l,mfprime,mfprime,Mw,s)+2.*s*C2minus_l(mu,l,mfprime,mfprime,Mw,s)
            -4.*(4.*mf*mf-s)*C1plus_l(mu,l,mfprime,mfprime,Mw,s)+(6.*mf*mf-2.*s)* C0)));
    
}


complex EWSMOneLoopLEP2::FI_Vbq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double lVp = -Qqprime(q)/4./sW2;
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    return (SM.getAle()/4./M_PI *(lVp*(4.*C20_q(mu,q,mfprime,mfprime,Mw,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_q(mu,q,mfprime,mfprime,Mw,s)+2.*s*C2minus_q(mu,q,mfprime,mfprime,Mw,s)
            -4.*(4.*mf*mf-s)*C1plus_q(mu,q,mfprime,mfprime,Mw,s)+(6.*mf*mf-2.*s)* C0)-
            2.*mfprime*mfprime*lVp*C0) );
    
}

complex EWSMOneLoopLEP2::FII_Vcl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = -I3/2./sW2;
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    complex C0 = -PV.C0(s,Mw,0.,Mw);
   
    
    return (SM.getAle()/4./M_PI *(lV* (12.*C20_l(mu,l,Mw,Mw,mfprime,s)-2.
            +2.*(4.*mf*mf-s)*C2plus_l(mu,l,Mw,Mw,mfprime,s)
            +2.*s*C2minus_l(mu,l,Mw,Mw,mfprime,s)-4.*(4.*mf*mf-s)*C1plus_l(mu,l,Mw,Mw,mfprime,s) )      
            ) );
    
}//ok

complex EWSMOneLoopLEP2::FII_Vcq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = -I3/2./sW2;
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    
    return (SM.getAle()/4./M_PI *(lV* (12.*C20_q(mu, q, Mw,Mw,mfprime,s)-2.
            +2.*(4.*mf*mf-s)*C2plus_q(mu,q,Mw,Mw,mfprime,s)
            +2.*s*C2minus_q(mu,q,Mw,Mw,mfprime,s)-4.*(4.*mf*mf-s)*C1plus_q(mu,q,Mw,Mw,mfprime,s) )      
            ) );
    
}


complex EWSMOneLoopLEP2::FIII_Vdl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double Qf = SM.getLeptons(l).getCharge();
    double lVp = -Qf*muf(Mw_i,mf)*muf(Mw_i,mf);
    double MH = SM.getMHl();

    complex C0 = -PV.C0(s,mf,MH,mf);
   
    
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_l(mu,l,mf,mf,MH,s)-  0.5+
            (4.*mf*mf-s)*C2plus_l(mu,l,mf,mf,MH,s) + s*C2minus_l(mu,l,mf,mf,MH,s)
            -mf*mf*C0)-mf*mf*lVp*C0-2.*mf*mf*lVp*C1plus_l(mu,l,mf,mf,MH,s)) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Vdq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Qf = SM.getQuarks(q).getCharge();
    double mf=SM.getQuarks(q).getMass();
    double MH = SM.getMHl();
    double lVp =  -Qf*muf(Mw_i,mf)*muf(Mw_i,mf);
    
    complex C0 = -PV.C0(s,mf,MH,mf);
   
    
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_q(mu,q,mf,mf,MH,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mf,mf,MH,s) + s*C2minus_q(mu,q,mf,mf,MH,s)
            -mf*mf*C0)-mf*mf*lVp*C0-2.*mf*mf*lVp*C1plus_q(mu,q,mf,mf,MH,s)) ) ;
    
}


complex EWSMOneLoopLEP2::FIII_Vel(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double Qf = SM.getLeptons(l).getCharge();
    double lVp = -Qf*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lVprime = -lVp;

    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
   
    
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_l(mu,l,mf,mf,Mz,s)-  0.5+
            (4.*mf*mf-s)*C2plus_l(mu,l,mf,mf,Mz,s) + s*C2minus_l(mu,l,mf,mf,Mz,s)
            -mf*mf*C0)-mf*mf*lVp*C0-2.*mf*mf*lVprime*C1plus_l(mu,l,mf,mf,Mz,s)) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Veq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Qf = SM.getQuarks(q).getCharge();
    double mf=SM.getQuarks(q).getMass();
    double lVp =  -Qf*muf(Mw_i,mf)*muf(Mw_i,mf);
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
   
    
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_q(mu,q,mf,mf,Mz,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mf,mf,Mz,s) + s*C2minus_q(mu,q,mf,mf,Mz,s)
            -mf*mf*C0)-mf*mf*lVp*C0-2.*mf*mf*lVp*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Vfl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    
    return (0.) ;
    
}

complex EWSMOneLoopLEP2::FIII_Vfq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double lVp =  -Qqprime(q)*(muf(Mw_i,mf)*muf(Mw_i,mf)+muf(Mw_i,mfprime)*muf(Mw_i,mfprime));
    double lVprime = 2.*Qqprime(q)*muf(Mw_i,mf)*muf(Mw_i,mfprime);
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
  
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_q(mu,q,mfprime,mfprime,Mw,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mfprime,mfprime,Mw,s) + s*C2minus_q(mu,q,mfprime,mfprime,Mw,s)
            -mf*mf*C0)-mfprime*mfprime*lVp*C0-2.*mf*mfprime*lVprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIV_Vgl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = -2.*I3*(muf(Mw_i,mf)*muf(Mw_i,mf));

    return (SM.getAle()/4./M_PI * (2.*lV*C20_l(mu,l,Mw,Mw,mfprime,s))) ;
    
}

complex EWSMOneLoopLEP2::FIV_Vgq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = -2.*I3*(muf(Mw_i,mf)*muf(Mw_i,mf)+muf(Mw_i,mfprime)*muf(Mw_i,mfprime));
  
    return (SM.getAle()/4./M_PI * (2.*lV*C20_q(mu,q,Mw,Mw,mfprime,s))) ;
    
}

complex EWSMOneLoopLEP2::FV_Vhl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = -I3*mf/sW2/2.;

    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    return (SM.getAle()/4./M_PI * ( -2.*mf*lV*C1minus_l(mu,l,Mw,Mw,mfprime,s)  )) ;
    
}

complex EWSMOneLoopLEP2::FV_Vhq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = -I3*mf/sW2/2.;
    double lVprime= I3*mfprime/sW2/2.;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( -2.*mf*lV*C1minus_q(mu,q,Mw,Mw,mfprime,s)
            +mfprime*lVprime*C0  )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Vil(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = -I3*mf/sW2/2.;

    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    return (SM.getAle()/4./M_PI * ( 2.*mf*lV*C1minus_l(mu,l,Mw,Mw,mfprime,s)  )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Viq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = -I3*mf/sW2/2.;
    double lVprime= I3*mfprime/sW2/2.;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( 2.*mf*lV*C1minus_q(mu,q,Mw,Mw,mfprime,s)
            +mfprime*lVprime*C0  )) ;
    
}






complex EWSMOneLoopLEP2::FI_Vjl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double vl2 = vl(l,Mw_i)*vl(l,Mw_i);
    double al2 = al(l,Mw_i)*al(l,Mw_i);
    double lVp = vl(l,Mw_i)*(vl2+3.*al2);
    double lVm = vl(l,Mw_i)*(vl2-al2);
    
    double mf=SM.getLeptons(l).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(4.*C20_l(mu,l,mf,mf,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_l(mu,l,mf,mf,Mz,s)+2.*s*C2minus_l(mu,l,mf,mf,Mz,s)
            -4.*(4.*mf*mf-s)*C1plus_l(mu,l,mf,mf,Mz,s)+(6.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lVm*C0) );
    
}

complex EWSMOneLoopLEP2::FI_Vjq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double vq2 = vq(q,Mw_i)*vq(q,Mw_i);
    double aq2 = aq(q,Mw_i)*aq(q,Mw_i);
    double lVp = vq(q,Mw_i)*(vq2+3.*aq2);
    double lVm = vq(q,Mw_i)*(vq2-aq2);
    double mf=SM.getQuarks(q).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(4.*C20_q(mu,q,mf,mf,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_q(mu,q,mf,mf,Mz,s)+2.*s*C2minus_q(mu,q,mf,mf,Mz,s)
            -4.*(4.*mf*mf-s)*C1plus_q(mu,q,mf,mf,Mz,s)+(6.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lVm*C0) );
    
}


complex EWSMOneLoopLEP2::FI_Vkl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3fprime = 1./2.;
    double lVp = I3fprime/sW/cW/4./sW2;
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    
    return (SM.getAle()/4./M_PI *(lVp*(4.*C20_l(mu,l,mfprime,mfprime,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_l(mu,l,mfprime,mfprime,Mz,s)+2.*s*C2minus_l(mu,l,mfprime,mfprime,Mz,s)
            -4.*(4.*mf*mf-s)*C1plus_l(mu,l,mfprime,mfprime,Mz,s)+(6.*mf*mf-2.*s)* C0)) );
    
}


complex EWSMOneLoopLEP2::FI_Vkq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double lVp = (I3qprime(q)-Qqprime(q)*sW2)/sW/cW/4./sW2;
    double lVm = -Qqprime(q)*sW/cW;
    double mf = SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    return (SM.getAle()/4./M_PI *(lVp*(4.*C20_q(mu,q,mfprime,mfprime,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_q(mu,q,mfprime,mfprime,Mz,s)+2.*s*C2minus_q(mu,q,mfprime,mfprime,Mz,s)
            -4.*(4.*mf*mf-s)*C1plus_q(mu,q,mfprime,mfprime,Mz,s)+(6.*mf*mf-2.*s)* C0)-
            2.*mfprime*mfprime*lVm*C0) );
    
}

complex EWSMOneLoopLEP2::FII_Vll(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = cW*I3/2./sW2/sW;
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    complex C0 = -PV.C0(s,Mw,0.,Mw);
    
    return (SM.getAle()/4./M_PI *(lV* (12.*C20_l(mu, l, Mw,Mw,mfprime,s)-2.
            +2.*(4.*mf*mf-s)*C2plus_l(mu,l,Mw,Mw,mfprime,s)
            +2.*s*C2minus_l(mu,l,Mw,Mw,mfprime,s)-4.*(4.*mf*mf-s)*C1plus_l(mu,l,Mw,Mw,mfprime,s) )      
            ) );
    
}

complex EWSMOneLoopLEP2::FII_Vlq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = cW*I3/2./sW2/sW;
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    
    return (SM.getAle()/4./M_PI *(lV* (12.*C20_q(mu, q, Mw,Mw,mfprime,s)-2.
            +2.*(4.*mf*mf-s)*C2plus_q(mu,q,Mw,Mw,mfprime,s)
            +2.*s*C2minus_q(mu,q,Mw,Mw,mfprime,s)-4.*(4.*mf*mf-s)*C1plus_q(mu,q,Mw,Mw,mfprime,s) )      
            ) );
    
}


complex EWSMOneLoopLEP2::FIII_Vml(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double lV = vl(l,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double MH = SM.getMHl();

    complex C0 = -PV.C0(s,mf,MH,mf);
    
    return (SM.getAle()/4./M_PI *(lV*(2.*C20_l(mu,l,mf,mf,MH,s)-  0.5+
            (4.*mf*mf-s)*C2plus_l(mu,l,mf,mf,MH,s) + s*C2minus_l(mu,l,mf,mf,MH,s)
            -mf*mf*C0)-mf*mf*lV*C0-2.*mf*mf*lV*C1plus_l(mu,l,mf,mf,MH,s)) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Vmq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getQuarks(q).getMass();
    double MH = SM.getMHl();
    double lVp =  vq(q,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    
    complex C0 = -PV.C0(s,mf,MH,mf);
   
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_q(mu,q,mf,mf,MH,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mf,mf,MH,s) + s*C2minus_q(mu,q,mf,mf,MH,s)
            -mf*mf*C0)-mf*mf*lVp*C0-2.*mf*mf*lVp*C1plus_q(mu,q,mf,mf,MH,s)) ) ;
    
}


complex EWSMOneLoopLEP2::FIII_Vnl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double lVp = vl(l,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lVprime = -lVp;

    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_l(mu,l,mf,mf,Mz,s)-  0.5+
            (4.*mf*mf-s)*C2plus_l(mu,l,mf,mf,Mz,s) + s*C2minus_l(mu,l,mf,mf,Mz,s)
            -mf*mf*C0)-mf*mf*lVp*C0-2.*mf*mf*lVprime*C1plus_l(mu,l,mf,mf,Mz,s)) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Vnq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getQuarks(q).getMass();
    double lVp =  vq(q,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lVprime = -lVp;
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
   
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_q(mu,q,mf,mf,Mz,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mf,mf,Mz,s) + s*C2minus_q(mu,q,mf,mf,Mz,s)
            -mf*mf*C0)-mf*mf*lVp*C0-2.*mf*mf*lVprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Vol(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double I3fprime = 1./2.;
    
    double lVp = muf(Mw_i,mf)*(I3fprime/sW/cW/4./sW2);
    double lVm = 0.;
    double lVprime = 0.;

    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_l(mu,l,mfprime,mfprime,Mw,s)-  0.5+
            (4.*mf*mf-s)*C2plus_l(mu,l,mfprime,mfprime,Mw,s) + s*C2minus_l(mu,l,mfprime,mfprime,Mw,s)
            -mf*mf*C0)-mfprime*mfprime*lVm*C0-2.*mf*mfprime*lVprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Voq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double lVp = muf(Mw_i,mf)*muf(Mw_i,mf)*(I3qprime(q)-sW2*Qqprime(q))/sW/cW 
                 +muf(Mw_i,mfprime)*muf(Mw_i,mfprime)* (-Qqprime(q)*sW/cW);
    double lVm = muf(Mw_i,mf)*muf(Mw_i,mf)*(-Qqprime(q)*sW/cW)
            +muf(Mw_i,mfprime)*muf(Mw_i,mfprime)*(I3qprime(q)-sW2*Qqprime(q))/sW/cW;
    double lVprime = -2.* 0.5*((I3qprime(q)-sW2*Qqprime(q))/sW/cW
                     -Qqprime(q)*sW/cW)*muf(Mw_i,mf)*muf(Mw_i,mfprime);
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
  
    return (SM.getAle()/4./M_PI *(lVp*(2.*C20_q(mu,q,mfprime,mfprime,Mw,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mfprime,mfprime,Mw,s) + s*C2minus_q(mu,q,mfprime,mfprime,Mw,s)
            -mf*mf*C0)-mfprime*mfprime*lVm*C0-2.*mf*mfprime*lVprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIV_Vpl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lV = -2.*alphaf_l(Mw_i,l)*(muf(Mw_i,mf)*muf(Mw_i,mf));

    return (SM.getAle()/4./M_PI * (2.*lV*C20_l(mu,l,Mw,Mw,mfprime,s))) ;
    
}

complex EWSMOneLoopLEP2::FIV_Vpq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lV = -2.*alphaf_q(Mw_i,q)*(muf(Mw_i,mf)*muf(Mw_i,mf)
                +muf(Mw_i,mfprime)*muf(Mw_i,mfprime));
    
    return (SM.getAle()/4./M_PI * (2.*lV*C20_q(mu,q,Mw,Mw,mfprime,s))) ;
    
}

complex EWSMOneLoopLEP2::FV_Vsl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lV = -al(l,Mw_i)*mf;

    return (SM.getAle()/4./M_PI * ( -2.*mf*lV*C1minus_l(mu,l,Mw,Mw,mfprime,s)  )) ;
    
}

complex EWSMOneLoopLEP2::FV_Vsq(const QCD::quark q, const double Mw_i, 
                                 const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lV = -aq(q,Mw_i)*mf;
    double lVprime= -aq(qprime(q),Mw_i)*mfprime;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( -2.*mf*lV*C1minus_q(mu,q,Mw,Mw,mfprime,s)
            +mfprime*lVprime*C0  )) ;
    
}

complex EWSMOneLoopLEP2::FV_Vtl(const StandardModel::lepton l, const double Mw_i, 
                                 const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double MH = SM.getMHl();
    double mf=SM.getLeptons(l).getMass();
    //mf = 0.;///MAYBE NOW IT IS OK
    double lV = -vl(l,Mw_i)*mf/2./sW2/cW2;

    complex C0 = -PV.C0(s,Mz,mf,MH);
   
    return (SM.getAle()/4./M_PI * ( -2.*mf*lV*C1minus_l(mu,l,Mz,MH,mf,s)  )) ;
    
}//NOT OK 

complex EWSMOneLoopLEP2::FV_Vtq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    //mf = 0.;//MAYBE NOW IT IS OK
    double lV = -vq(q,Mw_i)*mf/2./sW2/cW2;
    complex C0 = -PV.C0(s,Mz,mf,MH);
  
    return (SM.getAle()/4./M_PI * ( -2.*mf*lV*C1minus_q(mu,q,Mz,MH,mf,s)
            +mfprime*lV*C0  )) ;
    
}


complex EWSMOneLoopLEP2::FVI_Vul(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lV = -al(l,Mw_i)*mf;

    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    return (SM.getAle()/4./M_PI * ( 2.*mf*lV*C1minus_l(mu,l,Mw,Mw,mfprime,s)  )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Vuq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lV = -aq(q,Mw_i)*mf;
    double lVprime= -aq(qprime(q),Mw_i)*mfprime;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( 2.*mf*lV*C1minus_q(mu,q,Mw,Mw,mfprime,s)
            +mfprime*lVprime*C0  )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Vvl(const StandardModel::lepton l, const double Mw_i, 
                                 const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getLeptons(l).getMass();
    double lV = -vl(l,Mw_i)*mf/2./sW2/cW2;
    //mf = 0.;//MAYBE NOW IT IS OK

    complex C0 = -PV.C0(s,MH,mf,Mz);
   
    return (SM.getAle()/4./M_PI * ( 2.*mf*lV*C1minus_l(mu,l,MH,Mz,mf,s)  )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Vvq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double lV = -vq(q,Mw_i)*mf/2./sW2/cW2;
    //mf = 0.;//MAYBE NOW IT IS OK
    
    complex C0 = -PV.C0(s,MH,mf,Mz);
  
    return (SM.getAle()/4./M_PI * ( 2.*mf*lV*C1minus_q(mu,q,MH,Mz,mf,s)
            +mfprime*lV*C0  )) ;
    
}

complex EWSMOneLoopLEP2::FVgammal_weak(const StandardModel::lepton l, 
                     const double Mw_i, const double mu, const double s) const{
    
    return (CVgammal(l,Mw_i,mu)+ FI_Val(l,Mw_i,mu,s)+FI_Vbl(l,Mw_i,mu,s)
            +FII_Vcl(l,Mw_i,mu,s)+FIII_Vdl(l,Mw_i,mu,s)+FIII_Vel(l,Mw_i,mu,s)
            +FIII_Vfl(l,Mw_i,mu,s)+FIV_Vgl(l,Mw_i,mu,s)+FV_Vhl(l,Mw_i,mu,s)
            +FVI_Vil(l,Mw_i,mu,s));
}


complex EWSMOneLoopLEP2::FVZl_weak(const StandardModel::lepton l, 
                     const double Mw_i, const double mu, const double s) const{

return (CVZl(l,Mw_i,mu)+FI_Vjl(l,Mw_i,mu,s)+FI_Vkl(l,Mw_i,mu,s)
            +FII_Vll(l,Mw_i,mu,s)+FIII_Vml(l,Mw_i,mu,s)+FIII_Vnl(l,Mw_i,mu,s)
            +FIII_Vol(l,Mw_i,mu,s)+FIV_Vpl(l,Mw_i,mu,s)+FV_Vsl(l,Mw_i,mu,s)
            +FV_Vtl(l,Mw_i,mu,s)+FVI_Vul(l,Mw_i,mu,s)+FVI_Vvl(l,Mw_i,mu,s));

}

complex EWSMOneLoopLEP2::FVgammaq_weak(const QCD::quark q, const double Mw_i,
                      const double mu, const double s) const{
    
    return (CVgammaq(q,Mw_i,mu)+ FI_Vaq(q,Mw_i,mu,s)+FI_Vbq(q,Mw_i,mu,s)
            +FII_Vcq(q,Mw_i,mu,s)+FIII_Vdq(q,Mw_i,mu,s)+FIII_Veq(q,Mw_i,mu,s)
            +FIII_Vfq(q,Mw_i,mu,s)+FIV_Vgq(q,Mw_i,mu,s)+FV_Vhq(q,Mw_i,mu,s)
            +FVI_Viq(q,Mw_i,mu,s));
}


complex EWSMOneLoopLEP2::FVZq_weak(const QCD::quark q, const double Mw_i,
                      const double mu, const double s) const{

return (CVZq(q,Mw_i,mu)+FI_Vjq(q,Mw_i,mu,s)+FI_Vkq(q,Mw_i,mu,s)
            +FII_Vlq(q,Mw_i,mu,s)+FIII_Vmq(q,Mw_i,mu,s)+FIII_Vnq(q,Mw_i,mu,s)
            +FIII_Voq(q,Mw_i,mu,s)+FIV_Vpq(q,Mw_i,mu,s)+FV_Vsq(q,Mw_i,mu,s)
            +FV_Vtq(q,Mw_i,mu,s)+FVI_Vuq(q,Mw_i,mu,s)+FVI_Vvq(q,Mw_i,mu,s));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////   MAGNETIC FORM FACTOR      /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////



complex EWSMOneLoopLEP2::FI_Mal(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double Qf = SM.getLeptons(l).getCharge();
    double vl2 = vl(l,Mw_i)*vl(l,Mw_i);
    double al2 = al(l,Mw_i)*al(l,Mw_i);
    double lVp = -Qf*(vl2+al2);
    double lVprime = -Qf*(vl2-al2);
    double mf=SM.getLeptons(l).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(-8.*mf*C2plus_l(mu,l,mf,mf,Mz,s)
            +12.*mf*C1plus_l(mu,l,mf,mf,Mz,s)-4.*mf* C0)
            +lVprime*(-8.*mf*C1plus_l(mu,l,mf,mf,Mz,s)+4.*mf*C0) ));
    
}

complex EWSMOneLoopLEP2::FI_Maq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double Qf = SM.getQuarks(q).getCharge();
    double vq2 = vq(q,Mw_i)*vq(q,Mw_i);
    double aq2 = aq(q,Mw_i)*aq(q,Mw_i);
    double lVp = -Qf*(vq2+aq2);
    double lVprime = -Qf*(vq2-aq2);
    double mf=SM.getQuarks(q).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(-8.*mf*C2plus_q(mu,q,mf,mf,Mz,s)
            +12.*mf*C1plus_q(mu,q,mf,mf,Mz,s)-4.*mf* C0)
            +lVprime*(-8.*mf*C1plus_q(mu,q,mf,mf,Mz,s)+4.*mf*C0)));
    
}


complex EWSMOneLoopLEP2::FI_Mbl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{

    return (0.);
    
}


complex EWSMOneLoopLEP2::FI_Mbq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double lVp = -Qqprime(q)/4./sW2;
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    return (SM.getAle()/4./M_PI *(lVp*(-8.*mf*C2plus_q(mu,q,mfprime,mfprime,Mw,s)
            +12.*mf*C1plus_q(mu,q,mfprime,mfprime,Mw,s)-4.*mf* C0)));
    
}

/////

complex EWSMOneLoopLEP2::FII_Mcl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = -I3/2./sW2;
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    complex C0 = -PV.C0(s,Mw,0.,Mw);
   
    
    return (SM.getAle()/4./M_PI *(lV* (8.*mf*C2plus_l(mu,l,Mw,Mw,mfprime,s)+
            2.*mf*C1plus_l(mu,l,Mw,Mw,mfprime,s) ) ) );
    
}

complex EWSMOneLoopLEP2::FII_Mcq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = -I3/2./sW2;
    double mf=SM.getQuarks(q).getMass();
    double mfprime = SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    
    return (SM.getAle()/4./M_PI *(lV* (8.*mf*C2plus_q(mu,q,Mw,Mw,mfprime,s)+
            2.*mf*C1plus_q(mu,q,Mw,Mw,mfprime,s) ) ) );
    
}


complex EWSMOneLoopLEP2::FIII_Mdl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double Qf = SM.getLeptons(l).getCharge();
    double lVp = -Qf*muf(Mw_i,mf)*muf(Mw_i,mf);
    double MH = SM.getMHl();

    complex C0 = -PV.C0(s,mf,MH,mf);
   
    
    return (SM.getAle()/4./M_PI *lVp*(-4.*mf*C2plus_l(mu,l,mf,mf,MH,s) 
            +4.*mf*C1plus_l(mu,l,mf,mf,MH,s)));
    
}

complex EWSMOneLoopLEP2::FIII_Mdq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Qf = SM.getQuarks(q).getCharge();
    double mf=SM.getQuarks(q).getMass();
    double MH = SM.getMHl();
    double lVp =  -Qf*muf(Mw_i,mf)*muf(Mw_i,mf);
    
    complex C0 = -PV.C0(s,mf,MH,mf);
   
    
    return (SM.getAle()/4./M_PI *lVp*(-4.*mf*C2plus_q(mu,q,mf,mf,MH,s) 
            +4.*mf*C1plus_q(mu,q,mf,mf,MH,s)));
    
}


complex EWSMOneLoopLEP2::FIII_Mel(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double Qf = SM.getLeptons(l).getCharge();
    double lVp = -Qf*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lVprime = -lVp;

    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
   
    
    return (SM.getAle()/4./M_PI *(lVp*(-4.*mf*C2plus_l(mu,l,mf,mf,Mz,s) 
            +2.*mf*C1plus_l(mu,l,mf,mf,Mz,s))+lVprime*2.*mf*C1plus_l(mu,l,mf,mf,Mz,s)));
    
}

complex EWSMOneLoopLEP2::FIII_Meq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Qf = SM.getQuarks(q).getCharge();
    double mf=SM.getQuarks(q).getMass();
    double lVp =  -Qf*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lVprime = -lVp;
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
   
    
    return (SM.getAle()/4./M_PI *(lVp*(-4.*mf*C2plus_q(mu,q,mf,mf,Mz,s) 
            +2.*mf*C1plus_q(mu,q,mf,mf,Mz,s))+lVprime*2.*mf*C1plus_q(mu,q,mf,mf,Mz,s) ));
    
}

complex EWSMOneLoopLEP2::FIII_Mfl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    
    return (0.) ;
    
}

complex EWSMOneLoopLEP2::FIII_Mfq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lVp =  -Qqprime(q)*(muf(Mw_i,mf)*muf(Mw_i,mf)
                  +muf(Mw_i,mfprime)*muf(Mw_i,mfprime));
    double lVprime = 2.*Qqprime(q)*muf(Mw_i,mf)*muf(Mw_i,mfprime);
    
  
    return (SM.getAle()/4./M_PI *(lVp*(-4.*mf*C2plus_q(mu,q,mfprime,mfprime,Mw,s) 
            +2.*mf*C1plus_q(mu,q,mfprime,mfprime,Mw,s))
            +lVprime*2.*mfprime*C1plus_q(mu,q,mfprime,mfprime,Mw,s)));
    
}

complex EWSMOneLoopLEP2::FIV_Mgl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = -2.*I3*(muf(Mw_i,mf)*muf(Mw_i,mf));

    return (SM.getAle()/4./M_PI * (lV*(4.*mf*C2plus_l(mu,l,Mw,Mw,mfprime,s)
            -2.*mf*C1plus_l(mu,l,Mw,Mw,mfprime,s))) );
    
}

complex EWSMOneLoopLEP2::FIV_Mgq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
     double mfprime=SM.getQuarks(qprime(q)).getMass();
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = -2.*I3*(muf(Mw_i,mf)*muf(Mw_i,mf)+muf(Mw_i,mfprime)*muf(Mw_i,mfprime));
    double lVprime = 4.*I3*muf(Mw_i,mf)*muf(Mw_i,mfprime);
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * (lV*(4.*mf*C2plus_q(mu,q,Mw,Mw,mfprime,s)
            -2.*mf*C1plus_q(mu,q,Mw,Mw,mfprime,s))
            +mfprime*lVprime*(2.*C1plus_q(mu,q,Mw,Mw,mfprime,s) - C0)  ) );
    
}

complex EWSMOneLoopLEP2::FV_Mhl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = -I3*mf/sW2/2.;

    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    return (SM.getAle()/4./M_PI * ( lV*(C1minus_l(mu,l,Mw,Mw,mfprime,s)
            + C1plus_l(mu,l,Mw,Mw,mfprime,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FV_Mhq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(q).getMass();
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = -I3*mf/sW2/2.;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( lV*(C1minus_q(mu,q,Mw,Mw,mfprime,s)
            + C1plus_q(mu,q,Mw,Mw,mfprime,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Mil(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = -I3*mf/sW2/2.;

    return (SM.getAle()/4./M_PI * ( lV*(C1minus_l(mu,l,Mw,Mw,mfprime,s)
            - C1plus_l(mu,l,Mw,Mw,mfprime,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Miq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = -I3*mf/sW2/2.;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( lV*(C1minus_q(mu,q,Mw,Mw,mfprime,s)
            - C1plus_q(mu,q,Mw,Mw,mfprime,s) ) )) ;
    
}




complex EWSMOneLoopLEP2::FI_Mjl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double vl2 = vl(l,Mw_i)*vl(l,Mw_i);
    double al2 = al(l,Mw_i)*al(l,Mw_i);
    double lVp = vl(l,Mw_i)*(vl2+3.*al2);
    double lVprime = vl(l,Mw_i)*(vl2-al2);
    
    double mf=SM.getLeptons(l).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(-8.*mf*C2plus_l(mu,l,mf,mf,Mz,s)
            +12.*mf*C1plus_l(mu,l,mf,mf,Mz,s)-4.*mf* C0)
            +lVprime*(-8.*mf*C1plus_l(mu,l,mf,mf,Mz,s)+4.*mf*C0) ));
    
}

complex EWSMOneLoopLEP2::FI_Mjq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double vq2 = vq(q,Mw_i)*vq(q,Mw_i);
    double aq2 = aq(q,Mw_i)*aq(q,Mw_i);
    double lVp = vq(q,Mw_i)*(vq2+3.*aq2);
    double lVprime = vq(q,Mw_i)*(vq2-aq2);
    double mf=SM.getQuarks(q).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(-8.*mf*C2plus_q(mu,q,mf,mf,Mz,s)
            +12.*mf*C1plus_q(mu,q,mf,mf,Mz,s)-4.*mf* C0)
            +lVprime*(-8.*mf*C1plus_q(mu,q,mf,mf,Mz,s)+4.*mf*C0)));
    
}


complex EWSMOneLoopLEP2::FI_Mkl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3fprime = 1./2.;
    double lVp = I3fprime/sW/cW/4./sW2;
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    complex C0 = -PV.C0(s,0.,Mw,0.);
    
    return (SM.getAle()/4./M_PI *(lVp*(-8.*mf*C2plus_l(mu,l,mfprime,mfprime,Mw,s)
            +12.*mf*C1plus_l(mu,l,mfprime,mfprime,Mw,s)-4.*mf* C0) ));
    
}


complex EWSMOneLoopLEP2::FI_Mkq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double lVp = (I3qprime(q)-Qqprime(q)*sW2)/sW/cW/4./sW2;
    double mf = SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    return (SM.getAle()/4./M_PI *(lVp*(-8.*mf*C2plus_q(mu,q,mfprime,mfprime,Mw,s)
            +12.*mf*C1plus_q(mu,q,mfprime,mfprime,Mw,s)-4.*mf* C0)));
    
}

complex EWSMOneLoopLEP2::FII_Mll(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3 = SM.getLeptons(l).getIsospin();
    double lV = cW*I3/2./sW2/sW;
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    complex C0 = -PV.C0(s,Mw,0.,Mw);
    
    return (SM.getAle()/4./M_PI *(lV* (8.*mf*C2plus_l(mu,l,Mw,Mw,mfprime,s)+
            2.*mf*C1plus_l(mu,l,Mw,Mw,mfprime,s) ) ) );
    
}

complex EWSMOneLoopLEP2::FII_Mlq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = cW*I3/2./sW2/sW;
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    
    return (SM.getAle()/4./M_PI *(lV* (8.*mf*C2plus_q(mu,q,Mw,Mw,mfprime,s)+
            2.*mf*C1plus_q(mu,q,Mw,Mw,mfprime,s) ) ) );
    
}


///////////////////


complex EWSMOneLoopLEP2::FIII_Mml(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double lV = vl(l,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double MH = SM.getMHl();

    complex C0 = -PV.C0(s,mf,MH,mf);
    
    return (SM.getAle()/4./M_PI *lV*(-4.*mf*C2plus_l(mu,l,mf,mf,MH,s) + 
            4.*mf*C1plus_l(mu,l,mf,mf,MH,s))) ;
    
}

complex EWSMOneLoopLEP2::FIII_Mmq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getQuarks(q).getMass();
    double MH = SM.getMHl();
    double lVp =  vq(q,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    
    complex C0 = -PV.C0(s,mf,MH,mf);
   
    return (SM.getAle()/4./M_PI *lVp*(-4.*mf*C2plus_q(mu,q,mf,mf,MH,s) 
            +4.*mf*C1plus_q(mu,q,mf,mf,MH,s)));
    
}


complex EWSMOneLoopLEP2::FIII_Mnl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double lVp = vl(l,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lVprime = -lVp;

    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lVp*(-4.*mf*C2plus_l(mu,l,mf,mf,Mz,s) 
            +2.*mf*C1plus_l(mu,l,mf,mf,Mz,s))+lVprime*2.*mf*C1plus_l(mu,l,mf,mf,Mz,s)));
    
}

complex EWSMOneLoopLEP2::FIII_Mnq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getQuarks(q).getMass();
    double lVp =  vq(q,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lVprime = -lVp;
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
   
    return (SM.getAle()/4./M_PI *(lVp*(-4.*mf*C2plus_q(mu,q,mf,mf,Mz,s) 
            +2.*mf*C1plus_q(mu,q,mf,mf,Mz,s))+lVprime*2.*mf*C1plus_q(mu,q,mf,mf,Mz,s) ));
    
}

complex EWSMOneLoopLEP2::FIII_Mol(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double I3fprime = 1./2.;
    
    double lVp = muf(Mw_i,mf)*(I3fprime/sW/cW/4./sW2);
    double lVprime = 0.;

    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    return (SM.getAle()/4./M_PI *(lVp*(-4.*mf*C2plus_l(mu,l,mfprime,mfprime,Mw,s) 
            +2.*mf*C1plus_l(mu,l,mfprime,mfprime,Mw,s))+lVprime*2.*mf*C1plus_l(mu,l,mfprime,mfprime,Mw,s)));
    
}

complex EWSMOneLoopLEP2::FIII_Moq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lVp = muf(Mw_i,mf)*muf(Mw_i,mf)*(I3qprime(q)-sW2*Qqprime(q))/sW/cW 
                 +muf(Mw_i,mfprime)*muf(Mw_i,mfprime)* (-Qqprime(q)*sW/cW);
    double lVprime = -2.* 0.5*((I3qprime(q)-sW2*Qqprime(q))/sW/cW
                     -Qqprime(q)*sW/cW)*muf(Mw_i,mf)*muf(Mw_i,mfprime);
  
    return (SM.getAle()/4./M_PI *(lVp*(-4.*mf*C2plus_q(mu,q,mfprime,mfprime,Mw,s) 
            +2.*mf*C1plus_q(mu,q,mfprime,mfprime,Mw,s))+lVprime*2.*mf*C1plus_q(mu,q,mfprime,mfprime,Mw,s) ));
    
}

complex EWSMOneLoopLEP2::FIV_Mpl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lV = -2.*alphaf_l(Mw_i,l)*(muf(Mw_i,mf)*muf(Mw_i,mf));

    return (SM.getAle()/4./M_PI * (lV*(4.*mf*C2plus_l(mu,l,Mw,Mw,mfprime,s)
            -2.*mf*C1plus_l(mu,l,Mw,Mw,mfprime,s))) );
    
}

complex EWSMOneLoopLEP2::FIV_Mpq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lV = -2.*alphaf_q(Mw_i,q)*(muf(Mw_i,mf)*muf(Mw_i,mf)
                +muf(Mw_i,mfprime)*muf(Mw_i,mfprime));
    double lVprime= 4.*alphaf_q(Mw_i,q)*muf(Mw_i,mf)*muf(Mw_i,mfprime);
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
    
    return (SM.getAle()/4./M_PI * (lV*(4.*mf*C2plus_q(mu,q,Mw,Mw,mfprime,s)
            -2.*mf*C1plus_q(mu,q,Mw,Mw,mfprime,s))
            +mfprime*lVprime*(2.*C1plus_q(mu,q,Mw,Mw,mfprime,s) - C0)  ) );
    
}

complex EWSMOneLoopLEP2::FV_Msl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lV = -al(l,Mw_i)*mf;

    return (SM.getAle()/4./M_PI * ( lV*(C1minus_l(mu,l,Mw,Mw,mfprime,s)
            + C1plus_l(mu,l,Mw,Mw,mfprime,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FV_Msq(const QCD::quark q, const double Mw_i, 
                                 const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lV = -aq(q,Mw_i)*mf;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( lV*(C1minus_q(mu,q,Mw,Mw,mfprime,s)
            + C1plus_q(mu,q,Mw,Mw,mfprime,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FV_Mtl(const StandardModel::lepton l, const double Mw_i, 
                                 const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getLeptons(l).getMass();
    double lV = -vl(l,Mw_i)*mf/2./sW2/cW2;
    //double mf = 0.;MAYBE NOW IT IS OK

    complex C0 = -PV.C0(s,Mz,mf,MH);
   
    return (SM.getAle()/4./M_PI * ( lV*(C1minus_l(mu,l,Mz,MH,mf,s)
            + C1plus_l(mu,l,Mz,MH,mf,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FV_Mtq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getQuarks(q).getMass();
    double lV = -vq(q,Mw_i)*mf/2./sW2/cW2;
    complex C0 = -PV.C0(s,Mz,mf,MH);
    //mf = 0.;MAYBE NOW IT IS OK
  
     return (SM.getAle()/4./M_PI * ( lV*(C1minus_q(mu,q,Mz,MH,mf,s)
            + C1plus_q(mu,q,Mz,MH,mf,s) ) )) ;
    
}


complex EWSMOneLoopLEP2::FVI_Mul(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lV = -al(l,Mw_i)*mf;

    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    return (SM.getAle()/4./M_PI * ( lV*(C1minus_l(mu,l,Mw,Mw,mfprime,s)
            - C1plus_l(mu,l,Mw,Mw,mfprime,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Muq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lV = -aq(q,Mw_i)*mf;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( lV*(C1minus_q(mu,q,Mw,Mw,mfprime,s)
            - C1plus_q(mu,q,Mw,Mw,mfprime,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Mvl(const StandardModel::lepton l, const double Mw_i, 
                                 const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getLeptons(l).getMass();
    double lV = -vl(l,Mw_i)*mf/2./sW2/cW2;

    return (SM.getAle()/4./M_PI * ( lV*(C1minus_l(mu,l,MH,Mz,mf,s)
            - C1plus_l(mu,l,MH,Mz,mf,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Mvq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getQuarks(q).getMass();
    double lV = -vq(q,Mw_i)*mf/2./sW2/cW2;
    complex C0 = -PV.C0(s,MH,mf,Mz);
    //mf = 0.;//MAYBE NOW IT IS OK
  
    return (SM.getAle()/4./M_PI * ( lV*(C1minus_q(mu,q,MH,Mz,mf,s)
            - C1plus_q(mu,q,MH,Mz,mf,s) ) )) ;
    
}

complex EWSMOneLoopLEP2::FMgammal_weak(const StandardModel::lepton l, 
                     const double Mw_i, const double mu, const double s) const{
    
    return (FI_Mal(l,Mw_i,mu,s)+FI_Mbl(l,Mw_i,mu,s)
            +FII_Mcl(l,Mw_i,mu,s)+FIII_Mdl(l,Mw_i,mu,s)+FIII_Mel(l,Mw_i,mu,s)
            +FIII_Mfl(l,Mw_i,mu,s)+FIV_Mgl(l,Mw_i,mu,s)+FV_Mhl(l,Mw_i,mu,s)
            +FVI_Mil(l,Mw_i,mu,s));
}


complex EWSMOneLoopLEP2::FMZl_weak(const StandardModel::lepton l, 
                     const double Mw_i, const double mu, const double s) const{

return (FI_Mjl(l,Mw_i,mu,s)+FI_Mkl(l,Mw_i,mu,s)
            +FII_Mll(l,Mw_i,mu,s)+FIII_Mml(l,Mw_i,mu,s)+FIII_Mnl(l,Mw_i,mu,s)
            +FIII_Mol(l,Mw_i,mu,s)+FIV_Mpl(l,Mw_i,mu,s)+FV_Msl(l,Mw_i,mu,s)
            +FV_Mtl(l,Mw_i,mu,s)+FVI_Mul(l,Mw_i,mu,s)+FVI_Mvl(l,Mw_i,mu,s));

}

complex EWSMOneLoopLEP2::FMgammaq_weak(const QCD::quark q, const double Mw_i,
                      const double mu, const double s) const{
    
    return (FI_Maq(q,Mw_i,mu,s)+FI_Mbq(q,Mw_i,mu,s)
            +FII_Mcq(q,Mw_i,mu,s)+FIII_Mdq(q,Mw_i,mu,s)+FIII_Meq(q,Mw_i,mu,s)
            +FIII_Mfq(q,Mw_i,mu,s)+FIV_Mgq(q,Mw_i,mu,s)+FV_Mhq(q,Mw_i,mu,s)
            +FVI_Miq(q,Mw_i,mu,s));
}


complex EWSMOneLoopLEP2::FMZq_weak(const QCD::quark q, const double Mw_i,
                      const double mu, const double s) const{

return (FI_Mjq(q,Mw_i,mu,s)+FI_Mkq(q,Mw_i,mu,s)
            +FII_Mlq(q,Mw_i,mu,s)+FIII_Mmq(q,Mw_i,mu,s)+FIII_Mnq(q,Mw_i,mu,s)
            +FIII_Moq(q,Mw_i,mu,s)+FIV_Mpq(q,Mw_i,mu,s)+FV_Msq(q,Mw_i,mu,s)
            +FV_Mtq(q,Mw_i,mu,s)+FVI_Muq(q,Mw_i,mu,s)+FVI_Mvq(q,Mw_i,mu,s));

}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////   AXIAL FORM FACTOR      /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////




complex EWSMOneLoopLEP2::FI_Aal(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double Qf = SM.getLeptons(l).getCharge();
    double lAp = -2.*Qf*vl(l,Mw_i)*al(l,Mw_i);
    double mf=SM.getLeptons(l).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lAp*(4.*C20_l(mu,l,mf,mf,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_l(mu,l,mf,mf,Mz,s)+2.*s*C2minus_l(mu,l,mf,mf,Mz,s)
            -4.*(2.*mf*mf-s)*C1plus_l(mu,l,mf,mf,Mz,s)+(2.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lAp*C0) );
    
}

complex EWSMOneLoopLEP2::FI_Aaq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double Qf = SM.getQuarks(q).getCharge();
    double lAp = -2.*Qf*vq(q,Mw_i)*aq(q,Mw_i);
    double mf=SM.getQuarks(q).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lAp*(4.*C20_q(mu,q,mf,mf,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_q(mu,q,mf,mf,Mz,s)+2.*s*C2minus_q(mu,q,mf,mf,Mz,s)
            -4.*(2.*mf*mf-s)*C1plus_q(mu,q,mf,mf,Mz,s)+(2.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lAp*C0) );
    
}


complex EWSMOneLoopLEP2::FI_Abl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double lAp = 0.;
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    complex C0 = -PV.C0(s,0.,Mw,0.);
    
    
    return (SM.getAle()/4./M_PI *(lAp*(4.*C20_l(mu,l,mfprime,mfprime,Mw,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_l(mu,l,mfprime,mfprime,Mw,s)+2.*s*C2minus_l(mu,l,mfprime,mfprime,Mw,s)
            -4.*(2.*mf*mf-s)*C1plus_l(mu,l,mfprime,mfprime,Mw,s)+(2.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lAp*C0) );
    
}


complex EWSMOneLoopLEP2::FI_Abq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double lAp = -Qqprime(q)/4./sW2;
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    return (SM.getAle()/4./M_PI *(lAp*(4.*C20_q(mu,q,mfprime,mfprime,Mw,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_q(mu,q,mfprime,mfprime,Mw,s)+2.*s*C2minus_q(mu,q,mfprime,mfprime,Mw,s)
            -4.*(2.*mf*mf-s)*C1plus_q(mu,q,mfprime,mfprime,Mw,s)+(2.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lAp*C0) );
    
}

complex EWSMOneLoopLEP2::FII_Acl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double I3 = SM.getLeptons(l).getIsospin();
    double lA = -I3/2./sW2;
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    complex C0 = -PV.C0(s,Mw,0.,Mw);
   
    
    return (SM.getAle()/4./M_PI *(lA* (12.*C20_l(mu,l,Mw,Mw,mfprime,s)-2.
            +2.*(4.*mf*mf-s)*C2plus_l(mu,l,Mw,Mw,mfprime,s)
            +2.*s*C2minus_l(mu,l,Mw,Mw,mfprime,s)-4.*(mf*mf-s)*C1plus_l(mu,l,Mw,Mw,mfprime,s))      
            ) );
    
}

complex EWSMOneLoopLEP2::FII_Acq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double I3 = SM.getQuarks(q).getIsospin();
    double lV = -I3/2./sW2;
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    
    return (SM.getAle()/4./M_PI *(lV* (12.*C20_q(mu, q, Mw,Mw,mfprime,s)-2.
            +2.*(4.*mf*mf-s)*C2plus_q(mu,q,Mw,Mw,mfprime,s)
            +2.*s*C2minus_q(mu,q,Mw,Mw,mfprime,s)-4.*(mf*mf-s)*C1plus_q(mu,q,Mw,Mw,mfprime,s) )      
            ) );
    
}



complex EWSMOneLoopLEP2::FIII_Afl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    
    return (0.) ;
    
}

complex EWSMOneLoopLEP2::FIII_Afq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double lAp =  Qqprime(q)*(muf(Mw_i,mf)*muf(Mw_i,mf)
                 -muf(Mw_i,mfprime)*muf(Mw_i,mfprime));
    double lAprime = 0.;
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
  
    return (SM.getAle()/4./M_PI *(lAp*(2.*C20_q(mu,q,mfprime,mfprime,Mw,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mfprime,mfprime,Mw,s) + s*C2minus_q(mu,q,mfprime,mfprime,Mw,s)
            +mf*mf*C0-4.*mf*mf*C1plus_q(mu,q,mfprime,mfprime,Mw,s))
            -mfprime*mfprime*lAp*C0+2.*mf*mfprime*lAprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIV_Agl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    double I3 = SM.getLeptons(l).getIsospin();
    double lA = 2.*I3*(muf(Mw_i,mf)*muf(Mw_i,mf));

    return (SM.getAle()/4./M_PI * (2.*lA*C20_l(mu,l,Mw,Mw,mfprime,s))) ;
    
}

complex EWSMOneLoopLEP2::FIV_Agq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double I3 = SM.getQuarks(q).getIsospin();
    double lA = 2.*I3*(muf(Mw_i,mf)*muf(Mw_i,mf)
                -muf(Mw_i,mfprime)*muf(Mw_i,mfprime));
  
    return (SM.getAle()/4./M_PI * (2.*lA*C20_q(mu,q,Mw,Mw,mfprime,s))) ;
    
}

complex EWSMOneLoopLEP2::FV_Ahl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double I3 = SM.getLeptons(l).getIsospin();
    double lA = -I3*mf/sW2/2.;
    double lAprime = 0.;

    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    return (SM.getAle()/4./M_PI * (2.*mf*lA*C1minus_l(mu,l,Mw,Mw,mfprime,s)
            +mfprime*lAprime*C0));
    
}

complex EWSMOneLoopLEP2::FV_Ahq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double I3 = SM.getQuarks(q).getIsospin();
    double lA = -I3*mf/sW2/2.;
    double lAprime= I3*mfprime/sW2/2.;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( 2.*mf*lA*C1minus_q(mu,q,Mw,Mw,mfprime,s)
            +mfprime*lAprime*C0  )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Ail(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double I3 = SM.getLeptons(l).getIsospin();
    double lA = I3*mf/sW2/2.;
    double lAprime = 0.;

    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    return (SM.getAle()/4./M_PI * ( -2.*mf*lA*C1minus_l(mu,l,Mw,Mw,mfprime,s) 
            + mfprime*lAprime*C0)) ;
    
}

complex EWSMOneLoopLEP2::FVI_Aiq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double I3 = SM.getQuarks(q).getIsospin();
    double lA = I3*mf/sW2/2.;
    double lAprime= I3*mfprime/sW2/2.;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( 2.*mf*lA*C1minus_q(mu,q,Mw,Mw,mfprime,s)
            +mfprime*lAprime*C0  )) ;
    
}



complex EWSMOneLoopLEP2::FI_Ajl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double vl2 = vl(l,Mw_i)*vl(l,Mw_i);
    double al2 = al(l,Mw_i)*al(l,Mw_i);
    double lAp = al(l,Mw_i)*(3.*vl2+al2);
    double lAm = al(l,Mw_i)*(-vl2+al2);
    
    double mf=SM.getLeptons(l).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lAp*(4.*C20_l(mu,l,mf,mf,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_l(mu,l,mf,mf,Mz,s)+2.*s*C2minus_l(mu,l,mf,mf,Mz,s)
            -4.*(2.*mf*mf-s)*C1plus_l(mu,l,mf,mf,Mz,s)+(2.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lAm*C0) );
    
}

complex EWSMOneLoopLEP2::FI_Ajq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double vq2 = vq(q,Mw_i)*vq(q,Mw_i);
    double aq2 = aq(q,Mw_i)*aq(q,Mw_i);
    double lAp = aq(q,Mw_i)*(3.*vq2+aq2);
    double lAm = aq(q,Mw_i)*(-vq2+aq2);
    double mf=SM.getQuarks(q).getMass();
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lAp*(4.*C20_q(mu,q,mf,mf,Mz,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_q(mu,q,mf,mf,Mz,s)+2.*s*C2minus_q(mu,q,mf,mf,Mz,s)
            -4.*(2.*mf*mf-s)*C1plus_q(mu,q,mf,mf,Mz,s)+(2.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lAm*C0) );
    
}


complex EWSMOneLoopLEP2::FI_Akl(const StandardModel::lepton l, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3fprime = 1./2.;
    double lAp = I3fprime/sW/cW/4./sW2;
    double lAm = 0.;
    double mf=SM.getLeptons(l).getMass();
    complex C0 = -PV.C0(s,0.,Mw,0.);
    double mfprime =0.;
    
    
    return (SM.getAle()/4./M_PI *(lAp*(4.*C20_l(mu,l,mfprime,mfprime,Mw,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_l(mu,l,mfprime,mfprime,Mw,s)+2.*s*C2minus_l(mu,l,mfprime,mfprime,Mw,s)
            -4.*(2.*mf*mf-s)*C1plus_l(mu,l,mfprime,mfprime,Mw,s)+(2.*mf*mf-2.*s)* C0)-
            2.*mf*mf*lAm*C0) );
    
}


complex EWSMOneLoopLEP2::FI_Akq(const QCD::quark q, const double Mw_i, const double mu,
                  const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double lAp = (I3qprime(q)-Qqprime(q)*sW2)/sW/cW/4./sW2;
    double lAm = -Qqprime(q)*sW/cW;
    double mf = SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    return (SM.getAle()/4./M_PI *(lAp*(4.*C20_q(mu,q,mfprime,mfprime,Mw,s)-2.
            +(8.*mf*mf-2.*s)*C2plus_q(mu,q,mfprime,mfprime,Mw,s)
            +2.*s*C2minus_q(mu,q,mfprime,mfprime,Mw,s)
            -4.*(2.*mf*mf-s)*C1plus_q(mu,q,mfprime,mfprime,Mw,s)
            +(2.*mf*mf-2.*s)* C0)-2.*mf*mf*lAm*C0) );
    
}

complex EWSMOneLoopLEP2::FII_All(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3 = SM.getLeptons(l).getIsospin();
    double lA = cW*I3/2./sW2/sW;
    double mf=SM.getLeptons(l).getMass();
    double mfprime = 0.;
    complex C0 = -PV.C0(s,Mw,0.,Mw);
    
    return (SM.getAle()/4./M_PI *(lA* (12.*C20_l(mu,l,Mw,Mw,mfprime,s)-2.
            +2.*(4.*mf*mf-s)*C2plus_l(mu,l,Mw,Mw,mfprime,s)
            +2.*s*C2minus_l(mu,l,Mw,Mw,mfprime,s)-4.*(mf*mf-s)*C1plus_l(mu,l,Mw,Mw,mfprime,s))      
            ) );
    
}

complex EWSMOneLoopLEP2::FII_Alq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double I3 = SM.getQuarks(q).getIsospin();
    double lA = cW*I3/2./sW2/sW;
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    
    return (SM.getAle()/4./M_PI *(lA* (12.*C20_q(mu, q, Mw,Mw,mfprime,s)-2.
            +2.*(4.*mf*mf-s)*C2plus_q(mu,q,Mw,Mw,mfprime,s)
            +2.*s*C2minus_q(mu,q,Mw,Mw,mfprime,s)-4.*(mf*mf-s)*C1plus_q(mu,q,Mw,Mw,mfprime,s) )      
            ) );
    
}


complex EWSMOneLoopLEP2::FIII_Aml(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double lAp = -al(l,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lAm = -lAp;
    double lAprime = lAm;
    double MH = SM.getMHl();

    complex C0 = -PV.C0(s,mf,MH,mf);
    
    return (SM.getAle()/4./M_PI *(lAp*(2.*C20_l(mu,l,mf,mf,MH,s)-  0.5+
            (4.*mf*mf-s)*C2plus_l(mu,l,mf,mf,MH,s) + s*C2minus_l(mu,l,mf,mf,MH,s)
            +mf*mf*C0-4.*mf*mf*C1plus_l(mu,l,mf,mf,MH,s))
            -mf*mf*lAm*C0+2.*mf*mf*lAprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Amq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getQuarks(q).getMass();
    double MH = SM.getMHl();
    double lAp = -aq(q,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lAm = -lAp;
    double lAprime = lAm;
    
    complex C0 = -PV.C0(s,mf,MH,mf);
   
    return (SM.getAle()/4./M_PI *(lAp*(2.*C20_q(mu,q,mf,mf,MH,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mf,mf,MH,s) + s*C2minus_q(mu,q,mf,mf,MH,s)
            +mf*mf*C0-4.*mf*mf*C1plus_q(mu,q,mf,mf,MH,s))
            -mf*mf*lAm*C0+2.*mf*mf*lAprime*C0) ) ;
    
}


complex EWSMOneLoopLEP2::FIII_Anl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getLeptons(l).getMass();
    double lAp = -al(l,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lAm = -lAp;
    double lAprime = lAp;

    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
    
    return (SM.getAle()/4./M_PI *(lAp*(2.*C20_l(mu,l,mf,mf,Mz,s)-  0.5+
            (4.*mf*mf-s)*C2plus_l(mu,l,mf,mf,Mz,s) + s*C2minus_l(mu,l,mf,mf,Mz,s)
            +mf*mf*C0-4.*mf*mf*C1plus_l(mu,l,mf,mf,Mz,s))
            -mf*mf*lAm*C0+2.*mf*mf*lAprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Anq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double mf=SM.getQuarks(q).getMass();
    double lAp = -aq(q,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    double lAm = -lAp;
    double lAprime = lAp;
    double Mz = SM.getMz();
    complex C0 = -PV.C0(s,mf,Mz,mf);
   
    return (SM.getAle()/4./M_PI *(lAp*(2.*C20_q(mu,q,mf,mf,Mz,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mf,mf,Mz,s) + s*C2minus_q(mu,q,mf,mf,Mz,s)
            +mf*mf*C0-4.*mf*mf*C1plus_q(mu,q,mf,mf,Mz,s))
            -mf*mf*lAm*C0+2.*mf*mf*lAprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Aol(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double I3fprime = 1./2.;
    
    double lAp = -muf(Mw_i,mf)*(I3fprime/sW/cW/4./sW2);
    double lAm = 0.;
    double lAprime = 0.;

    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
    
    return (SM.getAle()/4./M_PI *(lAp*(2.*C20_l(mu,l,mfprime,mfprime,Mw,s)-  0.5+
            (4.*mf*mf-s)*C2plus_l(mu,l,mfprime,mfprime,Mw,s) 
            + s*C2minus_l(mu,l,mfprime,mfprime,Mw,s)
            +mf*mf*C0-4.*mf*mf*C1plus_l(mu,l,mfprime,mfprime,Mw,s))
            -mfprime*mfprime*lAm*C0+2.*mf*mfprime*lAprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIII_Aoq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double sW = sqrt(sW2);
    double cW = sqrt(cW2);
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double lAp = -muf(Mw_i,mf)*muf(Mw_i,mf)*(I3qprime(q)-sW2*Qqprime(q))/sW/cW 
                 +muf(Mw_i,mfprime)*muf(Mw_i,mfprime)* (-Qqprime(q)*sW/cW);
    double lAm = muf(Mw_i,mf)*muf(Mw_i,mf)*(-Qqprime(q)*sW/cW)
            +muf(Mw_i,mfprime)*muf(Mw_i,mfprime)*(I3qprime(q)-sW2*Qqprime(q))/sW/cW;
    double lAprime = -2.* 0.5*((I3qprime(q)-sW2*Qqprime(q))/sW/cW
                     +Qqprime(q)*sW/cW)*muf(Mw_i,mf)*muf(Mw_i,mfprime);
    complex C0 = -PV.C0(s,mfprime,Mw,mfprime);
  
    return (SM.getAle()/4./M_PI *(lAp*(2.*C20_q(mu,q,mfprime,mfprime,Mw,s)-  0.5+
            (4.*mf*mf-s)*C2plus_q(mu,q,mfprime,mfprime,Mw,s) 
            + s*C2minus_q(mu,q,mfprime,mfprime,Mw,s)
            +mf*mf*C0-4.*mf*mf*C1plus_q(mu,q,mfprime,mfprime,Mw,s))
            -mfprime*mfprime*lAm*C0+2.*mfprime*mf*lAprime*C0) ) ;
    
}

complex EWSMOneLoopLEP2::FIV_Apl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lA = 2.*alphaf_l(Mw_i,l)*(muf(Mw_i,mf)*muf(Mw_i,mf));

    return (SM.getAle()/4./M_PI * (2.*lA*C20_l(mu,l,Mw,Mw,mfprime,s)));
    
}

complex EWSMOneLoopLEP2::FIV_Apq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double lA = 2.*alphaf_q(Mw_i,q)*(muf(Mw_i,mf)*muf(Mw_i,mf)
                -muf(Mw_i,mfprime)*muf(Mw_i,mfprime));
    
    return (SM.getAle()/4./M_PI * (2.*lA*C20_q(mu,q,Mw,Mw,mfprime,s)));
    
}

complex EWSMOneLoopLEP2::FIV_Aql(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double mf=SM.getLeptons(l).getMass();
    double lA = 2.*al(l,Mw_i)*(muf(Mw_i,mf)*muf(Mw_i,mf));

    return (SM.getAle()/4./M_PI * (2.*lA*C20_l(mu,l,Mz,MH,mf,s))) ;
    
}

complex EWSMOneLoopLEP2::FIV_Aqq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double mf=SM.getQuarks(q).getMass();
    double lA = 2.*aq(q,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    
    return (SM.getAle()/4./M_PI * (2.*lA*C20_q(mu,q,Mz,MH,mf,s))) ;
    
}
complex EWSMOneLoopLEP2::FIV_Arl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double mf=SM.getLeptons(l).getMass();
    double lA = 2.*al(l,Mw_i)*(muf(Mw_i,mf)*muf(Mw_i,mf));

    return (SM.getAle()/4./M_PI * (2.*lA*C20_l(mu,l,MH,Mz,mf,s))) ;
    
}

complex EWSMOneLoopLEP2::FIV_Arq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double mf=SM.getQuarks(q).getMass();
    double lA = 2.*aq(q,Mw_i)*muf(Mw_i,mf)*muf(Mw_i,mf);
    
    return (SM.getAle()/4./M_PI * (2.*lA*C20_q(mu,q,MH,Mz,mf,s))) ;
    
}

complex EWSMOneLoopLEP2::FV_Asl(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lA = -al(l,Mw_i)*mf;
    double lAprime = 0.;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);

    return (SM.getAle()/4./M_PI * (2.*mf*lA*C1minus_l(mu,l,Mw,Mw,mfprime,s)
            +mfprime*lAprime*C0));
    
}

complex EWSMOneLoopLEP2::FV_Asq(const QCD::quark q, const double Mw_i, 
                                 const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lA = -aq(q,Mw_i)*mf;
    double lAprime= -aq(qprime(q),Mw_i)*mfprime;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( 2.*mf*lA*C1minus_q(mu,q,Mw,Mw,mfprime,s)
            +mfprime*lAprime*C0  )) ;
    
}

complex EWSMOneLoopLEP2::FV_Atl(const StandardModel::lepton l, const double Mw_i, 
                                 const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lA = -al(l,Mw_i)*mf/2./sW2/cW2;
    double lAprime = lA;
    complex C0 = -PV.C0(s,Mz,mf,MH);
    //mf = 0.;//MAYBE NOW IT IS OK
   
    return (SM.getAle()/4./M_PI * (2.*mf*lA*C1minus_l(mu,l,Mz,MH,mf,s)
            +mfprime*lAprime*C0));
    
}

complex EWSMOneLoopLEP2::FV_Atq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getQuarks(q).getMass();
    double mfprime =SM.getQuarks(qprime(q)).getMass();
    double lA = -vq(q,Mw_i)*mf/2./sW2/cW2;
    double lAprime = lA;
    complex C0 = -PV.C0(s,Mz,mf,MH);
    //mf = 0.;//MAYBE NOW IT IS OK
  
    return (SM.getAle()/4./M_PI * ( 2.*mf*lA*C1minus_q(mu,q,Mz,MH,mf,s)
            +mfprime*lAprime*C0  )) ;
    
}


complex EWSMOneLoopLEP2::FVI_Aul(const StandardModel::lepton l, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mfprime = 0.;
    double mf=SM.getLeptons(l).getMass();
    double lA = al(l,Mw_i)*mf;
    double lAprime = 0.;

    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
   
    return (SM.getAle()/4./M_PI * ( -2.*mf*lA*C1minus_l(mu,l,Mw,Mw,mfprime,s) 
            + mfprime*lAprime*C0)) ;
    
}

complex EWSMOneLoopLEP2::FVI_Auq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double mf=SM.getQuarks(q).getMass();
    double mfprime=SM.getQuarks(qprime(q)).getMass();
    double lA = aq(q,Mw_i)*mf;
    double lAprime= -aq(qprime(q),Mw_i)*mfprime;
    complex C0 = -PV.C0(s,Mw,mfprime,Mw);
  
    return (SM.getAle()/4./M_PI * ( 2.*mf*lA*C1minus_q(mu,q,Mw,Mw,mfprime,s)
            +mfprime*lAprime*C0  )) ;
    
}

complex EWSMOneLoopLEP2::FVI_Avl(const StandardModel::lepton l, const double Mw_i, 
                                 const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getLeptons(l).getMass();
    double lA = al(l,Mw_i)*mf/2./sW2/cW2;
    double lAprime = -lA;
    //mf = 0.;//MAYBE NOW IT IS OK

    complex C0 = -PV.C0(s,MH,mf,Mz);
   
    return (SM.getAle()/4./M_PI * ( -2.*mf*lA*C1minus_l(mu,l,MH,Mz,mf,s) 
            + mf*lAprime*C0)) ;
    
}

complex EWSMOneLoopLEP2::FVI_Avq(const QCD::quark q, const double Mw_i, 
                                  const double mu, const double s) const{
    double Mw = SM.Mw();
    double Mz = SM.getMz();
    double MH = SM.getMHl();
    double sW2 = SM.sW2();
    double cW2 = SM.cW2();
    double mf=SM.getQuarks(q).getMass();
    double lA = -aq(q,Mw_i)*mf/2./sW2/cW2;
    double lAprime = lA;
    complex C0 = -PV.C0(s,MH,mf,Mz);
    //mf = 0.;//MAYBE NOW IT IS OK
  
    return (SM.getAle()/4./M_PI * ( 2.*mf*lA*C1minus_q(mu,q,MH,Mz,mf,s)
            +mf*lAprime*C0  )) ;
    
}

complex EWSMOneLoopLEP2::FAgammal_weak(const StandardModel::lepton l, 
                     const double Mw_i, const double mu, const double s) const{
    
    return (CAgammal(l,Mw_i,mu)+ FI_Aal(l,Mw_i,mu,s)+FI_Abl(l,Mw_i,mu,s)
            +FII_Acl(l,Mw_i,mu,s)+FIII_Afl(l,Mw_i,mu,s)+FIV_Agl(l,Mw_i,mu,s)
            +FV_Ahl(l,Mw_i,mu,s)+FVI_Ail(l,Mw_i,mu,s));
            
}


complex EWSMOneLoopLEP2::FAZl_weak(const StandardModel::lepton l, 
                     const double Mw_i, const double mu, const double s) const{

    return (CAZl(l,Mw_i,mu)+FI_Ajl(l,Mw_i,mu,s)+FI_Akl(l,Mw_i,mu,s)
            +FII_All(l,Mw_i,mu,s)+FIII_Aml(l,Mw_i,mu,s)+FIII_Anl(l,Mw_i,mu,s)
            +FIII_Aol(l,Mw_i,mu,s)+FIV_Apl(l,Mw_i,mu,s)+FIV_Aql(l,Mw_i,mu,s)
            +FIV_Arl(l,Mw_i,mu,s)+FV_Asl(l,Mw_i,mu,s)+FV_Atl(l,Mw_i,mu,s)
            +FVI_Aul(l,Mw_i,mu,s)+FVI_Avl(l,Mw_i,mu,s));

}

complex EWSMOneLoopLEP2::FAgammaq_weak(const QCD::quark q, const double Mw_i,
                      const double mu, const double s) const{
    
    return (CAgammaq(q,Mw_i,mu)+ FI_Aaq(q,Mw_i,mu,s)+FI_Abq(q,Mw_i,mu,s)
            +FII_Acq(q,Mw_i,mu,s)+FIII_Afq(q,Mw_i,mu,s)+FIV_Agq(q,Mw_i,mu,s)
            +FV_Ahq(q,Mw_i,mu,s) +FVI_Aiq(q,Mw_i,mu,s));
           
}

complex EWSMOneLoopLEP2::FAZq_weak(const QCD::quark q, const double Mw_i,
                      const double mu, const double s) const{

return (CAZq(q,Mw_i,mu)+FI_Ajq(q,Mw_i,mu,s)+FI_Akq(q,Mw_i,mu,s)
            +FII_Alq(q,Mw_i,mu,s)+FIII_Amq(q,Mw_i,mu,s)+FIII_Anq(q,Mw_i,mu,s)
            +FIII_Aoq(q,Mw_i,mu,s)+FIV_Apq(q,Mw_i,mu,s)+FIV_Aqq(q,Mw_i,mu,s)
            +FIV_Arq(q,Mw_i,mu,s)+FV_Asq(q,Mw_i,mu,s)+FV_Atq(q,Mw_i,mu,s)
            +FVI_Auq(q,Mw_i,mu,s)+FVI_Avq(q,Mw_i,mu,s));

}


complex EWSMOneLoopLEP2::E1(const double mu,const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const{
    
    return (Chi_Z(mu,s,Mw_i,W,X,Y)/s);
    
}

complex EWSMOneLoopLEP2::E2(const double mu,const double k, const double s, const double Mw_i,
                      const double W, const double X, const double Y) const{
    
    return (g_rhoe(k,Mw_i)*Chi_Z(mu,s,Mw_i,W,X,Y)/s);
    
}


complex EWSMOneLoopLEP2::F1_l(const double mu,const double rho, const double s,
                              const double Mw_i,const StandardModel::lepton l) const{
    
    return (FVgammal_weak(l,Mw_i,mu,s)-rho*FAgammal_weak(l,Mw_i,mu,s));
    
}

complex EWSMOneLoopLEP2::F2_l(const double mu,const double rho, const double s,
                              const double Mw_i,const StandardModel::lepton l) const{
    
    return (FVZl_weak(l,Mw_i,mu,s)-rho*FAZl_weak(l,Mw_i,mu,s));
    
}

complex EWSMOneLoopLEP2::F1_q(const double mu,const double rho, const double s,
                              const double Mw_i,const QCD::quark q) const{
    
    return (FVgammaq_weak(q,Mw_i,mu,s)-rho*FAgammaq_weak(q,Mw_i,mu,s));
    
}

complex EWSMOneLoopLEP2::F2_q(const double mu,const double rho, const double s,
                              const double Mw_i,const QCD::quark q) const{
    
    return (FVZq_weak(q,Mw_i,mu,s)-rho*FAZq_weak(q,Mw_i,mu,s));
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////        WEAK BOX CONTRIBUTION        /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////



complex EWSMOneLoopLEP2::A_CC(const double k, const double Mw_i) const{
    double Mw = SM.Mw();
    double sW2 = SM.sW2();
    double A;
    
    if(k == -0.5){
        A = 1./2./sW2/sW2;
    } else {
        A=0.;
    }
     
    return (A);
}

complex EWSMOneLoopLEP2::B_CCq(const double s,const double cos_theta,
                               const QCD::quark q, const double Mw_i) const{
    
     double I3 = SM.getQuarks(q).getIsospin();
     double Mw = SM.Mw();
     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
    
     if(q == SM.BOTTOM){
         mfprime = SM.getQuarks(SM.TOP).getMass();
     } else {
         mfprime = 0.;
     }
     
//     std::cout <<"PV.D27(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta) in B_CCl()= \t" 
//             << PV.D27(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta) << "\n\n"<<std::endl;
     
     complex x = PV.D27(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta);
     complex y = PV.D0(p12,p22,p32,p42,s,-t(0.,s,cos_theta),Mw,0.,Mw,mfprime);
     complex z = PV.C0(s,0.,Mw,Mw);
     complex a = PV.C0(s,Mw,mfprime,Mw);
     
     
     return (SM.getAle()*0.5/M_PI*(0.5-I3)*x-
             (0.5+I3)*SM.getAle()/4./M_PI*((Mw*Mw-t(0.,s,cos_theta))*y
             -z  //OTHER PROBLEM
             -a));

    
}


complex EWSMOneLoopLEP2::B_CCl(const double s,const double cos_theta,
                               const StandardModel::lepton l, const double Mw_i) const{
                               
    
     double I3 = SM.getLeptons(l).getIsospin();
     double Mw = SM.Mw();
     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
     
//     std::cout <<"PV.D27(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta) in B_CCl()= \t" 
//             << PV.D27(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta) << "\n\n"<<std::endl;
     
     complex x = PV.D27(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta);
     complex y = PV.D27(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta);
     
     return (SM.getAle()*0.5/M_PI*((0.5-I3)*x-
             (0.5+I3)*y
             ));
    
}

complex EWSMOneLoopLEP2::C_CCq(const double mu,const double s,const double cos_theta,
                               const QCD::quark q, const double Mw_i) const{
    
     double I3 = SM.getQuarks(q).getIsospin();
     double Mw = SM.Mw();
     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
    
     if(q == SM.BOTTOM){
         mfprime = SM.getQuarks(SM.TOP).getMass();
     } else {
         mfprime = 0.;
     }
     
//     std::cout <<"D24(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta) in C_CCq()= \t" 
//             << PV.D24(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta) << "\n\n"<<std::endl;
     
     complex x = PV.D11(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta);
     complex y = PV.D24(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta);
     complex z = PV.D25(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta);
     
     
     return (SM.getAle()*0.5/M_PI*(0.5-I3)*(x+
             y-
             z));
     
     return 0.;

    
}


complex EWSMOneLoopLEP2::C_CCl(const double mu,const double s,const double cos_theta,
                               const StandardModel::lepton l, const double Mw_i) const{
                               
    
     double I3 = SM.getLeptons(l).getIsospin();
     double Mw = SM.Mw();
     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
     
//     std::cout <<"D24(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta) in C_CCl()= \t" 
//             << PV.D24(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta) << "\n\n"<<std::endl;
     complex x = PV.D11(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta);
     complex y = PV.D24(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta);
     complex z = PV.D25(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta);
     
     
     return (SM.getAle()*0.5/M_PI*(0.5-I3)*(x+
             y-
             z));
     
     //return 0.;
    
}

complex EWSMOneLoopLEP2::D_CCq(const double mu,const double s,const double cos_theta,
                               const QCD::quark q, const double Mw_i) const{
    
     double I3 = SM.getQuarks(q).getIsospin();
     double Mw = SM.Mw();
     double mf = 0.;
     double mfprime;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
    
     if(q == SM.BOTTOM){
         mfprime = SM.getQuarks(SM.TOP).getMass();
     } else {
         mfprime = 0.;
     }
     
     return 0.;//because mf = 0.;
     
//     (SM.getAle()*mf*0.5/M_PI*(0.5-I3)*(PV.D12(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)+
//             PV.D22(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)-
//             PV.D26(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)));

    
}


complex EWSMOneLoopLEP2::D_CCl(const double mu,const double s,const double cos_theta,
                               const StandardModel::lepton l, const double Mw_i) const{
                               
    
     double I3 = SM.getLeptons(l).getIsospin();
     double Mw = SM.Mw();
     double mf = 0.;
     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
     
     return 0.;//because mf = 0.;
     
//     return (SM.getAle()*mf*0.5/M_PI*(0.5-I3)*(PV.D12(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)+
//             PV.D22(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)-
//             PV.D26(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)));
    
}

complex EWSMOneLoopLEP2::E_CCq(const double mu,const double s,const double cos_theta,
                               const QCD::quark q, const double Mw_i) const{
    
     double I3 = SM.getQuarks(q).getIsospin();
     double Mw = SM.Mw();
     double mf = 0.;
     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
    
     if(q == SM.BOTTOM){
         mfprime = SM.getQuarks(SM.TOP).getMass();
     } else {
         mfprime = 0.;
     }
     
     return 0.;//because mf = 0.;
     
//     (-SM.getAle()*mf*0.5/M_PI*(0.5-I3)*(PV.D13(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)+
//             PV.D26(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)-
//             PV.D23(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)));

    
}


complex EWSMOneLoopLEP2::E_CCl(const double mu,const double s,const double cos_theta,
                               const StandardModel::lepton l, const double Mw_i) const{
                               
    
     double I3 = SM.getLeptons(l).getIsospin();
     double Mw = SM.Mw();
     double mfprime = 0.;
     double mf = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
     
     return 0.;//because mf = 0.
     
//     (-SM.getAle()*mf*0.5/M_PI*(0.5-I3)*(PV.D13(p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)+
//             PV.D26(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)-
//             PV.D23(mu,p12,p22,p32,p42,s,Mw,0.,Mw,mfprime,cos_theta)));
    
}



complex EWSMOneLoopLEP2::A1_NCq(const double mu,const double s,const double cos_theta,
                               const QCD::quark q, const double Mw_i, const double rho, const double k) const{
    
//     double I3 = SM.getQuarks(q).getIsospin();
     double Mz = SM.getMz();
     double mf = 0.;
//     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
     complex F1;
//std::cout <<"F1_A1_NCq (prima di calcolarlo)= \t" << F1 << "\n\n"<<std::endl;

    F1=SM.getAle()/4./M_PI*((Mz*Mz-t(0.,s,cos_theta))*PV.D0(p12,p22,p32,p42,s,-t(0.,s,cos_theta),Mz,0.,Mz,mf)
               -PV.C0(s,0.,Mz,Mz) //OTHER PROBLEM
               -PV.C0(s,Mz,mf,Mz));
    
//    std::cout <<"F1_A1_NCq (dopo)= \t" << F1 << "\n\n"<<std::endl;
    //remember the sign - before C0 because of the different notation between 
    //Beenakker-vanderMarck-Hollik paperand Bardin-Passarino book
    complex F2=SM.getAle()*0.5/M_PI*PV.D27(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    
    complex F7=0.;//because mf = 0.
            
//            -2.*Mz*mf*SM.getAle()/4./M_PI*PV.D12(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    complex F8=0.;//because mf = 0.
//            -2.*Mz*mf*SM.getAle()/4./M_PI*PV.D13(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    
    double lkplus_q=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofq(q,k,Mw_i)*g_rhofq(q,k,Mw_i);
    double lkminus_q=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofq(q,-k,Mw_i)*g_rhofq(q,-k,Mw_i);
    double lbarkplus_q=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofq(q,k,Mw_i)*g_rhofq(q,-k,Mw_i);
    double lbarkminus_q=lbarkplus_q;
    
    complex B11kq = lkminus_q*(2.*F1-2.*F2)+lbarkplus_q*(F7+F7)-lbarkminus_q*(F8+F8);
    complex B12kq = lkplus_q*(2.*F2-2.*F1)+lbarkplus_q*(F8+F8)-lbarkminus_q*(F7+F7);
   
    if(rho == k){
        return (B11kq);
    } else if (rho == -k){
        return (B12kq);
    } else {
        throw "Invalid rho or k in EWSMOneLoopLEP2::A1_NCq()!!";
    }
}

complex EWSMOneLoopLEP2::A1_NCl(const double mu,const double s,const double cos_theta,
                               const StandardModel::lepton l, const double Mw_i, const double rho, const double k) const{
    
//     double I3 = SM.getQuarks(q).getIsospin();
     double Mz = SM.getMz();
     double mf = 0.;
//     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
     complex F1;

    complex x = PV.D0(p12,p22,p32,p42,s,-t(mf,s,cos_theta),Mz,0.,Mz,mf);
    complex y = PV.C0(s,0.,Mz,Mz);
    complex z = PV.C0(s,Mz,mf,Mz);
    complex a = PV.D27(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);




//std::cout <<"PV.D0(p12,p22,p32,p42,s,-t(mf,s,cos_theta),Mz,0.,Mz,mf) in A1_NCl (prima di calcolarlo)= \t" 
//          << PV.D0(p12,p22,p32,p42,s,-t(mf,s,cos_theta),Mz,0.,Mz,mf) << "\n\n"<<std::endl;
    if(t(mf,s,cos_theta) == 0.){
        F1 = SM.getAle()/4./M_PI*((mf*mf)*x
               -y //OTHER PROBLEM
               -z);
    } else {
        F1=SM.getAle()/4./M_PI*((mf*mf-t(mf,s,cos_theta))*x
               -y //OTHER PROBLEM
               -z);
    }
    
//    std::cout <<"F1_A1_NCl (dopo)= \t" << F1 << "\n\n"<<std::endl;
//    std::cout <<"PV.D0(p12,p22,p32,p42,s,-t(mf,s,cos_theta),Mz,0.,Mz,mf) in A1_NCl (dopo)= \t" 
//          << PV.D0(p12,p22,p32,p42,s,-t(mf,s,cos_theta),Mz,0.,Mz,mf) << "\n\n"<<std::endl;
    //remember the sign - before C0 because of the different notation between 
    //Beenakker-vanderMarck-Hollik paperand Bardin-Passarino book
    complex F2=SM.getAle()*0.5/M_PI*a;
    
    complex F7=0.;//because mf = 0.
//            -2.*Mz*mf*SM.getAle()/4./M_PI*PV.D12(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    complex F8=0.;//because mf = 0.
//            -2.*Mz*mf*SM.getAle()/4./M_PI*PV.D13(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    
    double lkplus_l=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofl(l,k,Mw_i)*g_rhofl(l,k,Mw_i);
    double lkminus_l=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofl(l,-k,Mw_i)*g_rhofl(l,-k,Mw_i);
    double lbarkplus_l=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofl(l,k,Mw_i)*g_rhofl(l,-k,Mw_i);
    double lbarkminus_l=lbarkplus_l;
    
    complex B11kl = lkminus_l*(2.*F1-2.*F2)+lbarkplus_l*(F7+F7)-lbarkminus_l*(F8+F8);
    complex B12kl = lkplus_l*(2.*F2-2.*F1)+lbarkplus_l*(F8+F8)-lbarkminus_l*(F7+F7);
    
    if(rho == k){
        return (B11kl);
    } else if (rho == -k){
        return (B12kl);
    } else {
        throw "Invalid rho or k in EWSMOneLoopLEP2::A1_NCl()!!";
    }
    
    //return 0.;
}

complex EWSMOneLoopLEP2::A2_NCq(const double mu,const double s,const double cos_theta,
                               const QCD::quark q, const double Mw_i, const double rho, const double k) const{
    
//     double I3 = SM.getQuarks(q).getIsospin();
     double Mz = SM.getMz();
     double mf = 0.;
//     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
     complex F3;
     //std::cout <<"F3_A2_NCq (prima di calcolarlo)= \t" << F3 << "\n\n"<<std::endl;
    
     complex x = PV.D11(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    complex y = PV.D24(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    complex z =  PV.D25(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    
    F3=SM.getAle()*0.5/M_PI*(x+y- z);
    
    //std::cout <<"F3_A2_NCq = \t" << F3 << "\n\n"<<std::endl;
    
    
    double lkplus_q=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofq(q,k,Mw_i)*g_rhofq(q,k,Mw_i);
    double lkminus_q=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofq(q,-k,Mw_i)*g_rhofq(q,-k,Mw_i);
    
    complex B21kq = 2.*lkminus_q*F3;
    complex B22kq = 2.*lkplus_q*F3;
   
    if(rho == k){
        return (B21kq);
    } else if (rho == -k){
        return (B22kq);
    } else {
        throw "Invalid rho or k in EWSMOneLoopLEP2::A2_NCq()!!";
    }
    
    
    //return 0.;
    
}

complex EWSMOneLoopLEP2::A2_NCl(const double mu,const double s,const double cos_theta,
                                const StandardModel::lepton l, const double Mw_i, const double rho, const double k) const{
    
//     double I3 = SM.getQuarks(q).getIsospin();
     double Mz = SM.getMz();
     double mf = 0.;
//     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;
     complex F3;
     //std::cout <<"F3_A2_NCl (prima di calcolarlo)= \t" << F3 << "\n\n"<<std::endl;
     complex x = PV.D11(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
     complex y = PV.D24(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
     complex z =  PV.D25(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);       
 
 
 
     F3=SM.getAle()*0.5/M_PI*(x + y - z);
    
    //std::cout <<"F3_A2_NCl (dopo averlo calcolato) = \t" << F3 << "\n\n"<<std::endl;
    
    double lkplus_l=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofl(l,k,Mw_i)*g_rhofl(l,k,Mw_i);
    double lkminus_l=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofl(l,-k,Mw_i)*g_rhofl(l,-k,Mw_i);
    
    complex B21kl = 2.*lkminus_l*F3;
    complex B22kl = 2.*lkplus_l*F3;
   
    if(rho == k){
        return (B21kl);
    } else if (rho == -k){
        return (B22kl);
    } else {
        throw "Invalid rho or k in EWSMOneLoopLEP2::A2_NCl()!!";
    }
    
    //return 0.;
}

complex EWSMOneLoopLEP2::A3_NCq(const double mu,const double s,const double cos_theta,
                               const QCD::quark q, const double Mw_i, const double rho, const double k) const{
    
//     double I3 = SM.getQuarks(q).getIsospin();
     double Mz = SM.getMz();
     double mf = 0.;
//     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;


    
    complex F4=0.;//because mf = 0.
//            SM.getAle()*mf*0.5/M_PI*(PV.D12(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)+
//             PV.D22(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)-
//             PV.D26(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta));
    complex F5=0.;//because mf = 0.
//            -SM.getAle()*mf*0.5/M_PI*(PV.D13(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)+
//               PV.D26(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)-
//               PV.D23(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta));
    complex F9=0.;//because mf = 0.
//            -2.*mf*SM.getAle()/4./M_PI*(PV.D12(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)-
//                PV.D13(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta));
    complex F10=0.;//because mf =0.
//            2.*mf*SM.getAle()/4./M_PI*(PV.D12(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)+
//                PV.D13(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta));
    
    
    double lkplus_q=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofq(q,k,Mw_i)*g_rhofq(q,k,Mw_i);
    double lkminus_q=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofq(q,-k,Mw_i)*g_rhofq(q,-k,Mw_i);
    double lbarkplus_q=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofq(q,k,Mw_i)*g_rhofq(q,-k,Mw_i);
    double lbarkminus_q=lbarkplus_q;
    
    complex B31kq = 2.*lkplus_q*F4-2.*lkminus_q*F4+lbarkplus_q*(F10-F9)-lbarkminus_q*(F10-F9);
    complex B32kq = 2.*lkplus_q*F5-2.*lkminus_q*F5+lbarkplus_q*(F10+F9)-lbarkminus_q*(F10+F9);
   
//    if(rho == k){
//        return (B31kq);
//    } else if (rho == -k){
//        return (B32kq);
//    } else {
//        throw "Invalid rho or k in EWSMOneLoopLEP2::A3_NCq()!!";
//    }
    
    return 0.;
    
}

complex EWSMOneLoopLEP2::A3_NCl(const double mu,const double s,const double cos_theta,
                                const StandardModel::lepton l, const double Mw_i, const double rho, const double k) const{
    
//     double I3 = SM.getQuarks(q).getIsospin();
     double Mz = SM.getMz();
     double mf = 0.;
//     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;


    
    complex F4=0.;//because mf = 0.
//            SM.getAle()*mf*0.5/M_PI*(PV.D12(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)+
//             PV.D22(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)-
//             PV.D26(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta));
    complex F5=0.;//because mf =0.
//            -SM.getAle()*mf*0.5/M_PI*(PV.D13(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)+
//               PV.D26(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)-
//               PV.D23(mu,p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta));
    complex F9=0.;//because mf =0.
//            -2.*mf*SM.getAle()/4./M_PI*(PV.D12(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)-
//                PV.D13(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta));
    complex F10=0.;//because mf =0.
//            2.*mf*SM.getAle()/4./M_PI*(PV.D12(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta)+
//                PV.D13(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta));
    
    double lkplus_l=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofl(l,k,Mw_i)*g_rhofl(l,k,Mw_i);
    double lkminus_l=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofl(l,-k,Mw_i)*g_rhofl(l,-k,Mw_i);
    double lbarkplus_l=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofl(l,k,Mw_i)*g_rhofl(l,-k,Mw_i);
    double lbarkminus_l=lbarkplus_l;
    
    complex B31kl = 2.*lkplus_l*F4-2.*lkminus_l*F4+lbarkplus_l*(F10-F9)-lbarkminus_l*(F10-F9);
    complex B32kl = 2.*lkplus_l*F5-2.*lkminus_l*F5+lbarkplus_l*(F10+F9)-lbarkminus_l*(F10+F9);
//   
//    if(rho == k){
//        return (B31kl);
//    } else if (rho == -k){
//        return (B32kl);
//    } else {
//        throw "Invalid rho or k in EWSMOneLoopLEP2::A3_NCl()!!";
//    }
    
    return 0.;
}

complex EWSMOneLoopLEP2::A4_NCq(const double mu,const double s,const double cos_theta,
                               const QCD::quark q, const double Mw_i, const double rho, const double k) const{
    
//     double I3 = SM.getQuarks(q).getIsospin();
     double Mz = SM.getMz();
     double mf = 0.;
//     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;

    complex F6=0.;//because mf =0.
            //-2.*Mz*mf*SM.getAle()/4./M_PI*PV.D11(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    
    double lbarkplus_q=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofq(q,k,Mw_i)*g_rhofq(q,-k,Mw_i);
    double lbarkminus_q=lbarkplus_q;
    
    complex B41kq = lbarkplus_q*F6-lbarkminus_q*F6;
    complex B42kq = lbarkplus_q*F6-lbarkminus_q*F6;
    
//    
//    if(rho == k){
//        return (B41kq);
//    } else if (rho == -k){
//        return (B42kq);
//    } else {
//        throw "Invalid rho or k in EWSMOneLoopLEP2::A4_NCq()!!";
//    }
    
    return 0.;
}

complex EWSMOneLoopLEP2::A4_NCl(const double mu,const double s,const double cos_theta,
                                const StandardModel::lepton l, const double Mw_i, const double rho, const double k) const{
    
//     double I3 = SM.getQuarks(q).getIsospin();
     double Mz = SM.getMz();
     double mf = 0.;
//     double mfprime = 0.;
     double p12=0.;
     double p22=0.;
     double p32=0.;
     double p42=0.;


    
    complex F6=0.;   //because mf =0.

            //-2.*Mz*mf*SM.getAle()/4./M_PI*PV.D11(p12,p22,p32,p42,s,Mz,0.,Mz,mf,cos_theta);
    
    double lbarkplus_l=g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*g_rhofl(l,k,Mw_i)*g_rhofl(l,-k,Mw_i);
    double lbarkminus_l=lbarkplus_l;
    
    complex B41kl = lbarkplus_l*F6-lbarkminus_l*F6;
    complex B42kl = lbarkplus_l*F6-lbarkminus_l*F6;
   
//    if(rho == k){
//        return (B41kl);
//    } else if (rho == -k){
//        return (B42kl);
//    } else {
//        throw "Invalid rho or k in EWSMOneLoopLEP2::A4_NCl()!!";
//    }
    
    return 0.;
}




/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////        QED VERTEX AND BOX CONTRIBUTION        ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////



complex EWSMOneLoopLEP2::Lambda1(const double s) const{
    
    double me = SM.getLeptons(1).getMass();
    complex i = complex::i();
    double phmass = 10.;
    complex Lambda1_tmp;
    
    Lambda1_tmp = (2.*log(phmass*phmass/(-s-i*EPSILON)) + 
            log((-s-i*EPSILON)/me/me))*(log((-s-i*EPSILON)/me/me) - 1.)+
            2.*(log((-s-i*EPSILON)/(me*me)))+ 4.*(M_PI*M_PI/12.-1.);
                         
    return (Lambda1_tmp);
}


complex EWSMOneLoopLEP2::E3(const double s) const{
    
    double Qe = SM.getLeptons(1).getCharge();
    
    return (SM.getAle()/4./M_PI*Qe*Qe*Lambda1(s));
    
}


complex EWSMOneLoopLEP2::Lambdaq(const double s, const QCD::quark q, const double mu) const{
    
    double mf = SM.getQuarks(q).getMass();
    double Qf = SM.getQuarks(q).getCharge();
    complex i = complex::i();
    complex C0 = -PV.C0(s,mf,0.,mf);//OTHER PROBLEM
    double phmass = 10.;
    
    complex Lambdaq_tmp = (SM.getAle()/4./M_PI*Qf*Qf*(-2.*(s-2.*mf*mf)*C0-4.-
            3.*(PV.B0(mu,s,mf,mf)-2.)-2.*log(phmass*phmass/mf/mf)));
    
    return (Lambdaq_tmp);
    
}

complex EWSMOneLoopLEP2::Lambdal(const double s,const StandardModel::lepton l, const double mu) const{
    
    double mf = SM.getLeptons(l).getMass();
    complex i = complex::i();
    double Qf = SM.getLeptons(l).getCharge();
    complex C0 = -PV.C0(s,mf,0.,mf);//OTHER PROBLEM
    double phmass = 10.;
    
    complex Lambdal_tmp = (SM.getAle()/4./M_PI*Qf*Qf*(-2.*(s-2.*mf*mf)*C0-4.-
            3.*(PV.B0(mu,s,mf,mf)-2.)-2.*log(phmass*phmass/mf/mf)));
    
    return (Lambdal_tmp);
    
}


complex EWSMOneLoopLEP2::LambdaMq(const double s, const QCD::quark q, const double mu) const{
    
    double mf = SM.getQuarks(q).getMass();
    double Qf = SM.getQuarks(q).getCharge();
    
    complex LambdaMq_tmp = (SM.getAle()/4./M_PI*Qf*Qf*2.*mf/(-s+4.*mf*mf)*
            (PV.B0(mu,s,mf,mf)-2.));
    
    return (LambdaMq_tmp);
    
}


complex EWSMOneLoopLEP2::LambdaMl(const double s, const StandardModel::lepton l, const double mu) const{
    
    double mf = SM.getLeptons(l).getMass();
    double Qf = SM.getLeptons(l).getCharge();
    
    complex LambdaMl_tmp = (SM.getAle()/4./M_PI*Qf*Qf*2.*mf/(-s+4.*mf*mf)*
            (PV.B0(mu,s,mf,mf)-2.));
    
    return (LambdaMl_tmp);
    
}

complex EWSMOneLoopLEP2::E4l(const double s,const StandardModel::lepton l,
                              const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const{
    
    double Qf = SM.getLeptons(l).getCharge();
    double Qe = SM.getLeptons(1).getCharge();
    
    return (Qf*Qe*Chi_Z(mu,s,Mw_i,W,X,Y)/s);
    
}


complex EWSMOneLoopLEP2::E4q(const double s,const QCD::quark q, 
                              const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const{
    
    double Qf = SM.getQuarks(q).getCharge();
    double Qe = SM.getLeptons(1).getCharge();
    
    return (Qf*Qe*Chi_Z(mu,s,Mw_i,W,X,Y)/s);
    
}


complex EWSMOneLoopLEP2::E5(const double s, const double k, const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const{
    
    return (g_rhoe(k,Mw_i)*Chi_Z(mu,s,Mw_i,W,X,Y)/s);
    
}


complex EWSMOneLoopLEP2::F5l(const double s, const double rho, const double mu,
                             const double Mw_i,const StandardModel::lepton l) const{
    
    double mf = SM.getLeptons(1).getMass();
    
    return (Lambdal(s,l,mu)*g_rhofl(l,rho,Mw_i)-4.*mf*rho*LambdaMl(s,l,mu)*al(l,Mw_i));
    
}

complex EWSMOneLoopLEP2::F5q(const double s, const double rho, const double mu,
                             const double Mw_i,const QCD::quark q) const{
    
    double mf = SM.getQuarks(q).getMass();
    
    return (Lambdaq(s,q,mu)*g_rhofq(q,rho,Mw_i)-4.*mf*rho*LambdaMq(s,q,mu)*aq(q,Mw_i));
    
}


complex EWSMOneLoopLEP2::G5l(const double s, const double mu,
                             const double Mw_i,const StandardModel::lepton l) const{
       
    return (2.*LambdaMl(s,l,mu)*vl(l,Mw_i));
    
}

complex EWSMOneLoopLEP2::G5q(const double s, const double mu,
                             const double Mw_i,const QCD::quark q) const{
    
    return (2.*LambdaMq(s,q,mu)*vq(q,Mw_i));
    
}


complex EWSMOneLoopLEP2::Gfunc(const double s,const double t) const{
 
    complex i = complex::i();
    
    return (s/(2.*(s+t))*(log(t/(s+i*EPSILON))-
           s*(s+2.*t)/(4.*(s+t)*(s+t))*(log(t/(s+i*EPSILON))*log(t/(s+i*EPSILON))+M_PI*M_PI)));
    
}

complex EWSMOneLoopLEP2::A_gammagammaq(const double s, const QCD::quark q, const double cos_theta) const{
    
    complex i = complex::i();
    double mf = SM.getQuarks(q).getMass();
    
    return (SM.getAle()/2./M_PI*(Gfunc(s,t(mf,s,cos_theta))+Gfunc(s,u(mf,s,cos_theta))));
}

complex EWSMOneLoopLEP2::A_gammagammal(const double s, const StandardModel::lepton l, const double cos_theta) const{
    
    complex i = complex::i();
    double mf = SM.getLeptons(l).getMass();
    
    return (SM.getAle()/2./M_PI*(Gfunc(s,t(mf,s,cos_theta))+Gfunc(s,u(mf,s,cos_theta))));
}


complex EWSMOneLoopLEP2::V_gammagammaq(const double s, const QCD::quark q, const double cos_theta) const{
    
    complex i = complex::i();
    double mf = SM.getQuarks(q).getMass();
    double phmass = 10.;
    
    return (SM.getAle()/2./M_PI*(Gfunc(s,t(mf,s,cos_theta))-Gfunc(s,u(mf,s,cos_theta))
            +2.*log(phmass*phmass/(-s-i*EPSILON))*log(t(mf,s,cos_theta)/u(mf,s,cos_theta))));    
}

complex EWSMOneLoopLEP2::V_gammagammal(const double s,const StandardModel::lepton l,const double cos_theta) const{
    
    complex i = complex::i();
    double mf = SM.getLeptons(l).getMass();
    double phmass = 10.;
    
    return (SM.getAle()/2./M_PI*(Gfunc(s,t(mf,s,cos_theta))-Gfunc(s,u(mf,s,cos_theta))
            +2.*log(phmass*phmass/(-s-i*EPSILON))*log(t(mf,s,cos_theta)/u(mf,s,cos_theta))));    
}


complex EWSMOneLoopLEP2::Afunc(const double s,const double t,const double GammaZ) const{
 
    complex i = complex::i();
    complex M = SM.getMz()*SM.getMz()-i*SM.getMz()*GammaZ;
    complex arg[2];
            arg[0] = s/M/M;
            arg[1] = -t/M/M;
    complex Li2[2];
    for (int i=0; i<2; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
    }
    
    return ((s-M*M)/(s+t)*(log(t/(s-M*M))+M*M/s*log(1-s/M/M)+(s+2.*t+M*M)/(s+t)*
            (log(-t/M/M)*log((M*M-s)/(M*M+t))+Li2[0]-Li2[1])));
    
}

complex EWSMOneLoopLEP2::A_gammaZq(const double s,const QCD::quark q, const double cos_theta,const double GammaZ) const{
    
  double mf = SM.getQuarks(q).getMass();
    
    
    return (SM.getAle()/2./M_PI*(Afunc(s,t(mf,s,cos_theta),GammaZ)+Afunc(s,u(mf,s,cos_theta),GammaZ)));
}

complex EWSMOneLoopLEP2::A_gammaZl(const double s,const StandardModel::lepton l, const double cos_theta,const double GammaZ) const{
    
    double mf = SM.getLeptons(l).getMass();
    
    return (SM.getAle()/2./M_PI*(Afunc(s,t(mf,s,cos_theta),GammaZ)+Afunc(s,u(mf,s,cos_theta),GammaZ)));
}


complex EWSMOneLoopLEP2::V_gammaZq(const double s,const QCD::quark q, const double cos_theta,const double GammaZ) const{
    
    
    double mf = SM.getQuarks(q).getMass();
    complex i = complex::i();
    complex M = SM.getMz()*SM.getMz()-i*SM.getMz()*GammaZ;
    double phmass = 10.;
    complex arg[2];
            arg[0] = 1+M*M/t(mf,s,cos_theta);
            arg[1] = 1+M*M/u(mf,s,cos_theta);
    complex Li2[2];
    for (int i=0; i<2; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
    }
    
    return (SM.getAle()/2./M_PI*(Afunc(s,t(mf,s,cos_theta),GammaZ)-Afunc(s,u(mf,s,cos_theta),GammaZ)
            +2.*Li2[0]-2.*Li2[1]+4.*log(M*phmass/(M*M-s))*log(t(mf,s,cos_theta)/u(mf,s,cos_theta))));  
    
}

complex EWSMOneLoopLEP2::V_gammaZl(const double s,const StandardModel::lepton l, const double cos_theta,
        const double GammaZ) const{
    
    
    complex i = complex::i();
    double mf = SM.getLeptons(l).getMass();
    double phmass = 10.;
    complex M = SM.getMz()*SM.getMz()-i*SM.getMz()*GammaZ;
    complex arg[2];
            arg[0] = 1+M*M/t(mf,s,cos_theta);
            arg[1] = 1+M*M/u(mf,s,cos_theta);
    complex Li2[2];
    for (int i=0; i<2; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
    }
    
    return (SM.getAle()/2./M_PI*(Afunc(s,t(mf,s,cos_theta),GammaZ)-Afunc(s,u(mf,s,cos_theta),GammaZ)
            +2.*Li2[0]-2.*Li2[1]+4.*log(M*phmass/(M*M-s))*log(t(mf,s,cos_theta)/u(mf,s,cos_theta))));  
    
}


complex EWSMOneLoopLEP2::E6l(const double s,const StandardModel::lepton l,
                              const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const{
    
    double Qf = SM.getLeptons(l).getCharge();
    double Qe = SM.getLeptons(1).getCharge();
    
    return (Qf*Qf*Qe*Qe*Chi_gamma(mu,s,Mw_i,W,X,Y)/s);
    
}


complex EWSMOneLoopLEP2::E6q(const double s,const QCD::quark q, 
                              const double mu, const double Mw_i,
                      const double W, const double X, const double Y) const{
    
    double Qf = SM.getQuarks(q).getCharge();
    double Qe = SM.getLeptons(1).getCharge();
    
    return (Qf*Qf*Qe*Qe*Chi_gamma(mu,s,Mw_i,W,X,Y)/s);
    
}

complex EWSMOneLoopLEP2::F6rhoq(const double s,const double rho, const double k,const QCD::quark q, const double cos_theta) const{
     
    return (V_gammagammaq(s,q,cos_theta)+rho*k*A_gammagammaq(s,q,cos_theta));
    
}

complex EWSMOneLoopLEP2::F6rhol(const double s,const double rho, const double k,const StandardModel::lepton l, const double cos_theta) const{
     
    return (V_gammagammal(s,l,cos_theta)+rho*k*A_gammagammal(s,l,cos_theta));
    
}

complex EWSMOneLoopLEP2::F7rhol(const double s,const double rho, const double k,
                                const double Mw_i,const StandardModel::lepton l, const double cos_theta,const double GammaZ) const{
     
    return (g_rhofl(l,rho,Mw_i)*(V_gammaZl(s,l,cos_theta,GammaZ)+rho*k*A_gammaZl(s,l,cos_theta,GammaZ)));
    
}


complex EWSMOneLoopLEP2::F7rhoq(const double s,const double rho, const double k,
        const double Mw_i,const QCD::quark q, const double cos_theta,const double GammaZ) const{
     
    return (g_rhofq(q,rho,Mw_i)*(V_gammaZq(s,q,cos_theta,GammaZ)+rho*k*A_gammaZq(s,q,cos_theta,GammaZ)));
    
}




/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////       SOFT PHOTON APPROXIMATION OF BREHMSTRAHLUNG CROSS SECTION     ///////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

    
    
    
double EWSMOneLoopLEP2::delta(const double s) const{
    
    double omega = 0.15*sqrt(s);//cut
    
    return (2.*omega/sqrt(s));
    
}

double EWSMOneLoopLEP2::xf(const double s, const double mf) const{
    
    double beta = sqrt(1.-4.*mf*mf/s);
    
    return ((beta-1.)/(beta+1.));
    
}
    
    

double EWSMOneLoopLEP2::Be(const double s) const{
    
    double me = SM.getLeptons(1).getMass();
    
    return (log(s/me/me)-1.);
    
}


double EWSMOneLoopLEP2::Bf(const double s, const double mf) const{
    
    double beta = sqrt(1.-4.*mf*mf/s);
    
    return (-(s-2.*mf*mf)/(s*beta)*log(-xf(s,mf))-1.);
    
}

double EWSMOneLoopLEP2::Bint(const double cos_theta, const double mf, const double s) const{
  
    double beta = sqrt(1.-4.*mf*mf/s);
    
    return (2.*log((1.-beta*cos_theta)/(1+beta*cos_theta)));
    
}
    
    double EWSMOneLoopLEP2::gammaIR_q(const double s, const QCD::quark q, const double cos_theta) const{
        
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getQuarks(q).getCharge();
        double mf = SM.getQuarks(q).getMass();
        int A = 0.;//if 0 we subtract the interference term
        double phmass = 10.;
        
        return (-SM.getAle()/M_PI*log(phmass*phmass/s)*(Qe*Qe*Be(s)+Qf*Qf*Bf(s,mf)+A*Qe*Qf*Bint(cos_theta,mf,s)));
        
    }
    
    
   double EWSMOneLoopLEP2::gammaIR_l(const double s, const StandardModel::lepton l, const double cos_theta) const{
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getLeptons(l).getCharge();
        double mf = SM.getLeptons(l).getMass();
        int A = 0.;//if 0 we subtract the interference term
        double phmass = 10.;
        
        return (-SM.getAle()/M_PI*log(phmass*phmass/s)*(Qe*Qe*Be(s)+Qf*Qf*Bf(s,mf)+A*Qe*Qf*Bint(cos_theta,mf,s)));
    }
   
   
   double EWSMOneLoopLEP2::gammadelta_q(const double s, const QCD::quark q, const double cos_theta) const{
        
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getQuarks(q).getCharge();
        double mf = SM.getQuarks(q).getMass();
        int A = 0.;//if 0 we subtract the interference term
        
        return (2.*SM.getAle()/M_PI*log(delta(s))*(Qe*Qe*Be(s)+Qf*Qf*Bf(s,mf)+A*Qe*Qf*Bint(cos_theta,mf,s)));
        
    }
    
    
   double EWSMOneLoopLEP2::gammadelta_l(const double s, const StandardModel::lepton l, const double cos_theta) const{
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getLeptons(l).getCharge();
        double mf = SM.getLeptons(l).getMass();
        int A = 0.;//if 0 we subtract the interference term
        
        return (2.*SM.getAle()/M_PI*log(delta(s))*(Qe*Qe*Be(s)+Qf*Qf*Bf(s,mf)+A*Qe*Qf*Bint(cos_theta,mf,s)));
    }
   
   
   complex EWSMOneLoopLEP2::gammadeltaINT_q(const double s, const QCD::quark q, const double cos_theta,const double GammaZ) const{
        
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getQuarks(q).getCharge();
        double mf = SM.getQuarks(q).getMass();
        double d = delta(s);
        complex i = complex::i();
        complex M = SM.getMz()*SM.getMz()-i*SM.getMz()*GammaZ;
        int A = 0.;//if 0 we subtract the interference term
        
        return (2.*SM.getAle()/M_PI*(Qe*Qe*Be(s)*log(d*(s-M*M)/(s-s*d-M*M))
                +Qf*Qf*Bf(s,mf)*log(d)
                +0.5*A*Qe*Qf*Bint(cos_theta,mf,s)*log(d*d*(s-M*M)/(s-s*d-M*M))));
        
    }
    
    
   complex EWSMOneLoopLEP2::gammadeltaINT_l(const double s, const StandardModel::lepton l, const double cos_theta,const double GammaZ) const{
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getLeptons(l).getCharge();
        double mf = SM.getLeptons(l).getMass();
        complex i = complex::i();
        double d = delta(s);
        complex M = SM.getMz()*SM.getMz()-i*SM.getMz()*GammaZ;
        int A = 0.;//if 0 we subtract the interference term
        
        return (2.*SM.getAle()/M_PI*(Qe*Qe*Be(s)*log(d*(s-M*M)/(s-s*d-M*M))
                +Qf*Qf*Bf(s,mf)*log(d)
                +0.5*A*Qe*Qf*Bint(cos_theta,mf,s)*log(d*d*(s-M*M)/(s-s*d-M*M))));
    }
   
   
   
   double EWSMOneLoopLEP2::gammadeltaRES_q(const double s, const QCD::quark q, const double cos_theta,const double GammaZ) const{
        
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getQuarks(q).getCharge();
        double mf = SM.getQuarks(q).getMass();
        double d = delta(s);
        complex i = complex::i();
        complex M = SM.getMz()*SM.getMz()-i*SM.getMz()*GammaZ;
        complex x = d*(s-M*M)/(s-s*d-M*M);
        double x_mod = x.abs();
        complex y = d*d*(s-M*M)/(s-s*d-M*M);
        double y_mod = y.abs();
        int A = 0.;//if 0 we subtract the interference term
        
        
        return (2.*SM.getAle()/M_PI*(Qe*Qe*Be(s)*log(x_mod)
                +Qf*Qf*Bf(s,mf)*log(d)
                +A*Qe*Qf*Bint(cos_theta,mf,s)*log(y_mod)));
        
    }
    
    
   double EWSMOneLoopLEP2::gammadeltaRES_l(const double s, const StandardModel::lepton l, const double cos_theta,const double GammaZ) const{
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getLeptons(l).getCharge();
        double mf = SM.getLeptons(l).getMass();
        complex i = complex::i();
        double d = delta(s);
        complex M = SM.getMz()*SM.getMz()-i*SM.getMz()*GammaZ;
        complex x = d*(s-M*M)/(s-s*d-M*M)*SM.getMz();
        double x_mod = x.abs();
        complex y = d*d*(s-M*M)/(s-s*d-M*M);
        double y_mod = y.abs();
        int A = 0.;//if 0 we subtract the interference term
        
        return (2.*SM.getAle()/M_PI*(Qe*Qe*Be(s)*log(x_mod)
                +Qf*Qf*Bf(s,mf)*log(d)
                +A*Qe*Qf*Bint(cos_theta,mf,s)*log(y_mod)));
    }
    
    
   
   double EWSMOneLoopLEP2::gammatail_q(const double s, const QCD::quark q,const double GammaZ) const{
        
        double Qe = SM.getLeptons(1).getCharge();
        double Mz=SM.getMz();
        double x = (Mz*Mz-s+s*delta(s))/(Mz*GammaZ);
        double y = (Mz*Mz-s)/(Mz*GammaZ);
        
        return (2.*SM.getAle()/M_PI*Qe*Qe*Be(s)*((s-Mz*Mz)/Mz/GammaZ)*(
                atan(x)-atan(y)));
        
    }
    
    
   double EWSMOneLoopLEP2::gammatail_l(const double s, const StandardModel::lepton l,const double GammaZ) const{
        double Qe = SM.getLeptons(1).getCharge();
     
        double Mz=SM.getMz();
        double d = delta(s);
        double x = (Mz*Mz-s+s*d)/(Mz*GammaZ);
        double y = (Mz*Mz-s)/(Mz*GammaZ);
        
        return (2.*SM.getAle()/M_PI*Qe*Qe*Be(s)*((s-Mz*Mz)/Mz/GammaZ)*(
                atan(x)-atan(y)));
    }
   
   
   double EWSMOneLoopLEP2::gammafin_q(const double s, const QCD::quark q,const double cos_theta) const{
        
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getQuarks(q).getCharge();
        double mf = SM.getQuarks(q).getMass();
        double beta = sqrt(1.-4.*mf*mf/s);
        int A = 0.;//if 0 we subtract the interference term
        
        complex arg[5];
            arg[0] = 1.+xf(s,mf);
            arg[1] = 1.-(1-beta)/(1-beta*cos_theta);
            arg[2] = 1.-(1+beta)/(1-beta*cos_theta);
            arg[3] = 1.-(1-beta)/(1+beta*cos_theta);
            arg[4] = 1.-(1+beta)/(1+beta*cos_theta);
        complex Li2[5];
        for (int i=0; i<5; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
        }
        
        
        return (-SM.getAle()/M_PI*(Qe*Qe*(Be(s)*Be(s)*0.5+M_PI*M_PI/3.-0.5)
                +Qf*Qf*(1./beta*log(-xf(s,mf))+(s-2.*mf*mf)/(s*beta)*(2.*Li2[0].real()
                +0.5*log(-xf(s,mf))*log(-xf(s,mf))))+2.*A*Qe*Qf*(Li2[1].real()+
                Li2[2].real()-Li2[3].real()-Li2[4].real())));
        
    }
    
    
   double EWSMOneLoopLEP2::gammafin_l(const double s, const StandardModel::lepton l,const double cos_theta) const{
        double Qe = SM.getLeptons(1).getCharge();
        double Qf = SM.getLeptons(l).getCharge();
        double mf = SM.getLeptons(l).getMass();
        double beta = sqrt(1.-4.*mf*mf/s);
        int A = 0.;//if 0 we subtract the interference term
        
        complex arg[5];
            arg[0] = 1.+xf(s,mf);
            arg[1] = 1.-(1.-beta)/(1.-beta*cos_theta);
            arg[2] = 1.-(1.+beta)/(1.-beta*cos_theta);
            arg[3] = 1.-(1.-beta)/(1.+beta*cos_theta);
            arg[4] = 1.-(1.+beta)/(1.+beta*cos_theta);
        complex Li2[5];
        for (int i=0; i<5; i++) {
                gsl_sf_result re, im;
                gsl_sf_complex_dilog_xy_e(arg[i].real(), arg[i].imag(), &re, &im);
                Li2[i].real() = re.val;
                Li2[i].imag() = im.val;
       }
        
        
        return (-SM.getAle()/M_PI*(Qe*Qe*(Be(s)*Be(s)*0.5+M_PI*M_PI/3.-0.5)
                +Qf*Qf*(1./beta*log(-xf(s,mf))+(s-2.*mf*mf)/(s*beta)*(2.*Li2[0].real()+0.5*log(-xf(s,mf))*log(-xf(s,mf))))
                +2.*A*Qe*Qf*(Li2[1].real()+Li2[2].real()-Li2[3].real()-Li2[4].real())));
        
}
   
   
    
    double EWSMOneLoopLEP2::deltagammagamma_softl(const double s, const StandardModel::lepton l, const double cos_theta) const {
        
        return (gammaIR_l(s,l,cos_theta)+gammadelta_l(s,l,cos_theta)+gammafin_l(s,l,cos_theta));
        
    }
   
   
   double EWSMOneLoopLEP2::deltagammagamma_softq(const double s, const QCD::quark q, const double cos_theta) const {
        
        return (gammaIR_q(s,q,cos_theta)+gammadelta_q(s,q,cos_theta)+gammafin_q(s,q,cos_theta));
        
    }
   
   complex EWSMOneLoopLEP2::deltagammaZ_softl(const double s, const StandardModel::lepton l, const double cos_theta,const double GammaZ) const {
        
        return (gammaIR_l(s,l,cos_theta)+gammadeltaINT_l(s,l,cos_theta,GammaZ)+gammafin_l(s,l,cos_theta));
        
    }
   
   
   complex EWSMOneLoopLEP2::deltagammaZ_softq(const double s, const QCD::quark q, const double cos_theta,const double GammaZ) const {
        
        return (gammaIR_q(s,q,cos_theta)+gammadeltaINT_q(s,q,cos_theta,GammaZ)+gammafin_q(s,q,cos_theta));
        
    }
   
   
   double EWSMOneLoopLEP2::deltaZZ_softl(const double s, const StandardModel::lepton l, const double cos_theta,const double GammaZ) const {
        
        return (gammaIR_l(s,l,cos_theta)+gammadeltaRES_l(s,l,cos_theta,GammaZ)+gammafin_l(s,l,cos_theta)+gammatail_l(s,l,GammaZ));
        
    }
   
   
   double EWSMOneLoopLEP2::deltaZZ_softq(const double s, const QCD::quark q, const double cos_theta,const double GammaZ) const {
        
        return (gammaIR_q(s,q,cos_theta)+gammadeltaRES_q(s,q,cos_theta,GammaZ)+gammafin_q(s,q,cos_theta)+gammatail_q(s,q,GammaZ));
        
    }
   
    
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////       CHIRALITY AMPLITUDE     //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////



complex EWSMOneLoopLEP2::M1rhok_M1rhopk_l(const double s, const double rho,
                                        const double rhoprime,const double k, StandardModel::lepton l,const double cos_theta) const{
    
    double mf = SM.getLeptons(l).getMass(); 
    complex x;
    
    if(rho==rhoprime){
        if(rho==k){
            x = 4.*(u(mf,s,cos_theta)-mf*mf)*(u(mf,s,cos_theta)-mf*mf);
            //x = 4.*u(mf,s,cos_theta)*u(mf,s,cos_theta);
        } else if(rho==-k){
            x = 4.*(t(mf,s,cos_theta)-mf*mf)*(t(mf,s,cos_theta)-mf*mf);
            //x = 4.*t(mf,s,cos_theta)*t(mf,s,cos_theta);
        }
    } else if(rho==-rhoprime){
        x = 4.*mf*mf*s;
        //x = 0.;
    }
    
    return x;
    
}


complex EWSMOneLoopLEP2::M1rhok_M1rhopk_q(const double s, const double rho
                                        ,const double rhoprime,const double k, const QCD::quark q,const double cos_theta) const{
    
    double mf = SM.getQuarks(q).getMass(); 
    double x;
    
    if(rho==rhoprime){
        if(rho==k){
           x = 4.*(u(mf,s,cos_theta)-mf*mf)*(u(mf,s,cos_theta)-mf*mf);
        } else if(rho==-k){
            x = 4.*(t(mf,s,cos_theta)-mf*mf)*(t(mf,s,cos_theta)-mf*mf);
        }
    } else if(rho==-rhoprime){
        x = mf*mf*s;
    }
    
    return (x);
}

complex EWSMOneLoopLEP2::M2rhok_M1rhopk_l(const double s, const double rho
                                        ,const double rhoprime,const double k, StandardModel::lepton l,const double cos_theta) const{
    
    double mf = SM.getLeptons(l).getMass(); 
    double x;
    
    if(rho==rhoprime){
        if(rho==k){
            x = 2.*(u(mf,s,cos_theta)*t(mf,s,cos_theta)-mf*mf*mf*mf)*(u(mf,s,cos_theta)-mf*mf);
        } else if(rho==-k){
            x = -2.*(u(mf,s,cos_theta)*t(mf,s,cos_theta)-mf*mf*mf*mf)*(t(mf,s,cos_theta)-mf*mf);
        }
    } else if(rho==-rhoprime){
        return (0.);
    }
    return (x);
}



complex EWSMOneLoopLEP2::M2rhok_M1rhopk_q(const double s, const double rho
                                        ,const double rhoprime,const double k, const QCD::quark q,const double cos_theta) const{
    
    double mf = SM.getQuarks(q).getMass(); 
    double x;
    
    if(rho==rhoprime){
        if(rho==k){
            x = 2.*(u(mf,s,cos_theta)*t(mf,s,cos_theta)-mf*mf*mf*mf)*(u(mf,s,cos_theta)-mf*mf);
        } else if(rho==-k){
            x = -2.*(u(mf,s,cos_theta)*t(mf,s,cos_theta)-mf*mf*mf*mf)*(t(mf,s,cos_theta)-mf*mf);
        }
    } else if(rho==-rhoprime){
        x = 0.;
    }
    
    return (x);        
            
}

complex EWSMOneLoopLEP2::M3rhok_M1rhopk_l(const double s, const double rho
                                        ,const double rhoprime,const double k, StandardModel::lepton l,const double cos_theta) const{
    
    double mf = SM.getLeptons(l).getMass(); 
    
    
    return (-2.*mf*(u(mf,s,cos_theta)*t(mf,s,cos_theta)-mf*mf*mf*mf));
}



complex EWSMOneLoopLEP2::M3rhok_M1rhopk_q(const double s, const double rho
                                        ,const double rhoprime,const double k, const QCD::quark q,const double cos_theta) const{
    
    double mf = SM.getQuarks(q).getMass(); 
    
    return (-2.*mf*(u(mf,s,cos_theta)*t(mf,s,cos_theta)-mf*mf*mf*mf));
}

complex EWSMOneLoopLEP2::M4rhok_M1rhopk_l(const double s, const double rho
                                        ,const double rhoprime,const double k, StandardModel::lepton l,const double cos_theta) const{
    
    double mf = SM.getLeptons(l).getMass(); 
    
    double x;
    
    if(rho==-k){
        if(rhoprime==-k){
            x = 4.*mf*s*(t(mf,s,cos_theta)-mf*mf);
        } else if(rhoprime==k){
            x = 4.*(u(mf,s,cos_theta)-mf*mf);
        }
    } else if(rho==k){
        x = 0.;
    }
    
    return (x);

}

complex EWSMOneLoopLEP2::M4rhok_M1rhopk_q(const double s, const double rho
                                        ,const double rhoprime,const double k, const QCD::quark q,const double cos_theta) const{
    
    double mf = SM.getQuarks(q).getMass(); 
    double x;
    
    if(rho==-k){
        if(rhoprime==-k){
           x = 4.*mf*s*(t(mf,s,cos_theta)-mf*mf);
        } else if(rhoprime==k){
            x = 4.*(u(mf,s,cos_theta)-mf*mf);
        }
    } else if(rho==k){
        x = 0.;
    } 
    return (x);
}


complex EWSMOneLoopLEP2::ATOTq(const double mu,const QCD::quark q, const double rho,
                               const double k, const double s, const double Mw_i,const double cos_theta,
                               const double W,const double X,const double Y,const double GammaZ) const{
    
    
    return (1./s*Bq(mu,q,rho,k,s,Mw_i,W,X,Y)+Cq(mu,q,rho,k,s,Mw_i,W,X,Y)+Dq_rho(mu,q,rho,k,s,Mw_i,W,X,Y)+
            E1(mu,k,s,Mw_i,W,X,Y)*F1_q(mu,rho,s,Mw_i,q)+E2(mu,k,s,Mw_i,W,X,Y)*F2_q(mu,rho,s,Mw_i,q)
            +1./s*E3(s)*Aq(mu,q,rho,k,s,Mw_i,W,X,Y)+E4q(s,q,mu,Mw_i,W,X,Y)*Lambdaq(s,q,mu)
            +E5(s,k,mu,Mw_i,W,X,Y)*F5q(s,rho,mu,Mw_i,q)+E6q(s,q,mu,Mw_i,W,X,Y)*F6rhoq(s,rho,k,q,cos_theta)
            +g_rhoe(k,Mw_i)*E5(s,k,mu,Mw_i,W,X,Y)*F7rhoq(s,rho,k,Mw_i,q,cos_theta,GammaZ)+A1_NCq(mu,s,cos_theta,q,Mw_i,rho,k));
    
    
}


complex EWSMOneLoopLEP2::ATOTl(const double mu,const StandardModel::lepton l, const double rho,
                  const double k, const double s, const double Mw_i,const double cos_theta,
                               const double W,const double X,const double Y,const double GammaZ) const{
    
    return (1./s*Bl(mu,l,rho,k,s,Mw_i,W,X,Y)
            +Cl(mu,l,rho,k,s,Mw_i,W,X,Y)
            +Dl_rho(mu,l,rho,k,s,Mw_i,W,X,Y)
            +E1(mu,k,s,Mw_i,W,X,Y)*F1_l(mu,rho,s,Mw_i,l)
            +E2(mu,k,s,Mw_i,W,X,Y)*F2_l(mu,rho,s,Mw_i,l)
            +1./s*E3(s)*Al(mu,l,rho,k,s,Mw_i,W,X,Y)
            +E4l(s,l,mu,Mw_i,W,X,Y)*Lambdal(s,l,mu)
            +E5(s,k,mu,Mw_i,W,X,Y)*F5l(s,rho,mu,Mw_i,l)
            +E6l(s,l,mu,Mw_i,W,X,Y)*F6rhol(s,rho,k,l,cos_theta)
            +g_rhoe(k,Mw_i)*E5(s,k,mu,Mw_i,W,X,Y)*F7rhol(s,rho,k,Mw_i,l,cos_theta,GammaZ));
//            +A1_NCl(mu,s,cos_theta,l,Mw_i,rho,k));
}


complex EWSMOneLoopLEP2::BTOTq(const double mu,const QCD::quark q, const double rho,
                               const double k, const double s, const double Mw_i,const double cos_theta,
                               const double W,const double X,const double Y) const{
    
    return (
            //A3_NCq(mu,s,cos_theta,q,Mw_i,rho,k)+
            E1(mu,k,s,Mw_i,W,X,Y)*FMgammaq_weak(q,Mw_i,mu,s)
            +E2(mu,k,s,Mw_i,W,X,Y)*FMZq_weak(q,Mw_i,mu,s)+2.*E4q(s,q,mu,Mw_i,W,X,Y)*LambdaMq(s,q,mu)
            +E5(s,k,mu,Mw_i,W,X,Y)*G5q(s,mu,Mw_i,q));
    
}


complex EWSMOneLoopLEP2::BTOTl(const double mu,const StandardModel::lepton l, const double rho,
                  const double k, const double s, const double Mw_i,const double cos_theta,
                               const double W,const double X,const double Y) const{
    
    return (
            //A3_NCl(mu,s,cos_theta,l,Mw_i,rho,k)+
            E1(mu,k,s,Mw_i,W,X,Y)*FMgammal_weak(l,Mw_i,mu,s)
            +E2(mu,k,s,Mw_i,W,X,Y)*FMZl_weak(l,Mw_i,mu,s)+2.*E4l(s,l,mu,Mw_i,W,X,Y)*LambdaMl(s,l,mu)
            +E5(s,k,mu,Mw_i,W,X,Y)*G5l(s,mu,Mw_i,l));
    
}



double EWSMOneLoopLEP2::MTOTq_sq(const QCD::quark q,const double k,const double s,
                                  const double Mw_i,const double cos_theta,
                               const double W,const double X,const double Y,const double GammaZ) const{   
    double Mq;
    
    double mu = 10.;
    complex a1 = 1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*ATOTq(mu,q,0.5,k,s,Mw_i,cos_theta,W,X,Y,GammaZ).conjugate()*M1rhok_M1rhopk_q(s,0.5,0.5,k,q,cos_theta);
    complex a2 = 1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*ATOTq(mu,q,-0.5,k,s,Mw_i,cos_theta,W,X,Y,GammaZ).conjugate()*M1rhok_M1rhopk_q(s,0.5,-0.5,k,q,cos_theta);
    complex a3 = 1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*ATOTq(mu,q,0.5,k,s,Mw_i,cos_theta,W,X,Y,GammaZ).conjugate()*M1rhok_M1rhopk_q(s,-0.5,0.5,k,q,cos_theta);
    complex a4 =1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*ATOTq(mu,q,-0.5,k,s,Mw_i,cos_theta,W,X,Y,GammaZ).conjugate()*M1rhok_M1rhopk_q(s,-0.5,-0.5,k,q,cos_theta);
    complex a5 =1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*A2_NCq(mu,s,cos_theta,q,Mw_i,0.5,k).conjugate()*M2rhok_M1rhopk_q(s,0.5,0.5,k,q,cos_theta);
    complex a6 =1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*A2_NCq(mu,s,cos_theta,q,Mw_i,-0.5,k).conjugate()*M2rhok_M1rhopk_q(s,0.5,-0.5,k,q,cos_theta);
    complex a7 =1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*A2_NCq(mu,s,cos_theta,q,Mw_i,0.5,k).conjugate()*M2rhok_M1rhopk_q(s,-0.5,0.5,k,q,cos_theta);
    complex a8 =1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*A2_NCq(mu,s,cos_theta,q,Mw_i,-0.5,k).conjugate()*M2rhok_M1rhopk_q(s,-0.5,-0.5,k,q,cos_theta);
    complex a9 =1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*BTOTq(mu,q,0.5,k,s,Mw_i,cos_theta,W,X,Y).conjugate()*M3rhok_M1rhopk_q(s,0.5,0.5,k,q,cos_theta);
    complex a10 =1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*BTOTq(mu,q,-0.5,k,s,Mw_i,cos_theta,W,X,Y).conjugate()*M3rhok_M1rhopk_q(s,0.5,-0.5,k,q,cos_theta);
    complex a11 =1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*BTOTq(mu,q,0.5,k,s,Mw_i,cos_theta,W,X,Y).conjugate()*M3rhok_M1rhopk_q(s,-0.5,0.5,k,q,cos_theta);
    complex a12 =1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*BTOTq(mu,q,-0.5,k,s,Mw_i,cos_theta,W,X,Y).conjugate()*M3rhok_M1rhopk_q(s,-0.5,-0.5,k,q,cos_theta);
    complex a13 =1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*A4_NCq(mu,s,cos_theta,q,Mw_i,0.5,k).conjugate()*M4rhok_M1rhopk_q(s,0.5,0.5,k,q,cos_theta);
    complex a14 =1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*A4_NCq(mu,s,cos_theta,q,Mw_i,-0.5,k).conjugate()*M4rhok_M1rhopk_q(s,0.5,-0.5,k,q,cos_theta);
    complex a15 =1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*A4_NCq(mu,s,cos_theta,q,Mw_i,0.5,k).conjugate()*M4rhok_M1rhopk_q(s,-0.5,0.5,k,q,cos_theta);
    complex a16 =1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*A4_NCq(mu,s,cos_theta,q,Mw_i,-0.5,k).conjugate()*M4rhok_M1rhopk_q(s,-0.5,-0.5,k,q,cos_theta);
    complex a17=1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*B_CCq(s,cos_theta,q,Mw_i)).conjugate()*M1rhok_M1rhopk_q(s,0.5,-0.5,-0.5,q,cos_theta);
    complex a18=1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*B_CCq(s,cos_theta,q,Mw_i)).conjugate()*M1rhok_M1rhopk_q(s,-0.5,-0.5,-0.5,q,cos_theta);
    complex a19=1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*C_CCq(mu,s,cos_theta,q,Mw_i)).conjugate()*M2rhok_M1rhopk_q(s,0.5,-0.5,-0.5,q,cos_theta);
    complex a20=1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*C_CCq(mu,s,cos_theta,q,Mw_i)).conjugate()*M2rhok_M1rhopk_q(s,-0.5,-0.5,-0.5,q,cos_theta);
    complex a21 = 1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*D_CCq(mu,s,cos_theta,q,Mw_i)).conjugate()*M3rhok_M1rhopk_q(s,0.5,0.5,-0.5,q,cos_theta);
    complex a22 = 1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*D_CCq(mu,s,cos_theta,q,Mw_i)).conjugate()*M3rhok_M1rhopk_q(s,-0.5,0.5,-0.5,q,cos_theta);
    complex a23 =1./s*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*E_CCq(mu,s,cos_theta,q,Mw_i)).conjugate()*M3rhok_M1rhopk_q(s,0.5,-0.5,-0.5,q,cos_theta) ;
    complex a24 = 1./s*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*E_CCq(mu,s,cos_theta,q,Mw_i)).conjugate()*M3rhok_M1rhopk_q(s,-0.5,-0.5,-0.5,q,cos_theta);
            
    Mq = 1./s/s*(M1rhok_M1rhopk_q(s,0.5,0.5,k,q,cos_theta).real()*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y).abs2()+
         M1rhok_M1rhopk_q(s,-0.5,-0.5,k,q,cos_theta).real()*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y).abs2() +
         2.*(M1rhok_M1rhopk_q(s,0.5,-0.5,k,q,cos_theta)*Aq(mu,q,0.5,k,s,Mw_i,W,X,Y)*Aq(mu,q,-0.5,k,s,Mw_i,W,X,Y).conjugate()).real())+
         2.*(a1).real()+
         2.*(a2).real()+
         2.*(a3).real()+
         2.*(a4).real()+
         2.*(a5).real()+
         2.*(a6).real()+
         2.*(a7).real()+
         2.*(a8).real()+
         2.*(a9).real()+
         2.*(a10).real()+
         2.*(a11).real()+
         2.*(a12).real()+
         2.*(a13).real()+
         2.*(a14).real()+
         2.*(a15).real()+
         2.*(a16).real()+ 
         2.*(a17).real()+
         2.*(a18).real()+      
         2.*(a19).real()+
         2.*(a20).real()+    
         2.*(a21).real()+
         2.*(a22).real()+
         2.*(a23).real()+
         2.*(a24).real();
    
    return (Mq);
    
}


double EWSMOneLoopLEP2::MTOTl_sq(const StandardModel::lepton l, const double k,
                                const double s, const double Mw_i,const double cos_theta,
                               const double W,const double X,const double Y,const double GammaZ)const{   
    
    double mu = 10.;
    double Ml = 0.;
//    complex a1 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*ATOTl(mu,l,0.5,k,s,Mw_i,cos_theta,W,X,Y,GammaZ).conjugate()*M1rhok_M1rhopk_l(s,0.5,0.5,k,l,cos_theta);
//    complex a2 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*ATOTl(mu,l,-0.5,k,s,Mw_i,cos_theta,W,X,Y,GammaZ).conjugate()*M1rhok_M1rhopk_l(s,0.5,-0.5,k,l,cos_theta);
//    complex a3 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*ATOTl(mu,l,0.5,k,s,Mw_i,cos_theta,W,X,Y,GammaZ).conjugate()*M1rhok_M1rhopk_l(s,-0.5,0.5,k,l,cos_theta);
//    complex a4 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*ATOTl(mu,l,-0.5,k,s,Mw_i,cos_theta,W,X,Y,GammaZ).conjugate()*M1rhok_M1rhopk_l(s,-0.5,-0.5,k,l,cos_theta);
//    complex a5 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*A2_NCl(mu,s,cos_theta,l,Mw_i,0.5,k).conjugate()*M2rhok_M1rhopk_l(s,0.5,0.5,k,l,cos_theta);
//    complex a6 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*A2_NCl(mu,s,cos_theta,l,Mw_i,-0.5,k).conjugate()*M2rhok_M1rhopk_l(s,0.5,-0.5,k,l,cos_theta);
//    complex a7 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*A2_NCl(mu,s,cos_theta,l,Mw_i,0.5,k).conjugate()*M2rhok_M1rhopk_l(s,-0.5,0.5,k,l,cos_theta);
//    complex a8 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*A2_NCl(mu,s,cos_theta,l,Mw_i,-0.5,k).conjugate()*M2rhok_M1rhopk_l(s,-0.5,-0.5,k,l,cos_theta);
//    complex a9 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*BTOTl(mu,l,0.5,k,s,Mw_i,cos_theta,W,X,Y).conjugate()*M3rhok_M1rhopk_l(s,0.5,0.5,k,l,cos_theta);
//    complex a10 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*BTOTl(mu,l,-0.5,k,s,Mw_i,cos_theta,W,X,Y).conjugate()*M3rhok_M1rhopk_l(s,0.5,-0.5,k,l,cos_theta);
//    complex a11 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*BTOTl(mu,l,0.5,k,s,Mw_i,cos_theta,W,X,Y).conjugate()*M3rhok_M1rhopk_l(s,-0.5,0.5,k,l,cos_theta);
//    complex a12 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*BTOTl(mu,l,-0.5,k,s,Mw_i,cos_theta,W,X,Y).conjugate()*M3rhok_M1rhopk_l(s,-0.5,-0.5,k,l,cos_theta);
//    complex a13 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*A4_NCl(mu,s,cos_theta,l,Mw_i,0.5,k).conjugate()*M4rhok_M1rhopk_l(s,0.5,0.5,k,l,cos_theta);
//    complex a14 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*A4_NCl(mu,s,cos_theta,l,Mw_i,-0.5,k).conjugate()*M4rhok_M1rhopk_l(s,0.5,-0.5,k,l,cos_theta);
//    complex a15 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*A4_NCl(mu,s,cos_theta,l,Mw_i,0.5,k).conjugate()*M4rhok_M1rhopk_l(s,-0.5,0.5,k,l,cos_theta);
//    complex a16 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*A4_NCl(mu,s,cos_theta,l,Mw_i,-0.5,k).conjugate()*M4rhok_M1rhopk_l(s,-0.5,-0.5,k,l,cos_theta);
//    complex a17 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*B_CCl(s,cos_theta,l,Mw_i)).conjugate()*M1rhok_M1rhopk_l(s,0.5,-0.5,-0.5,l,cos_theta);
//    complex a18 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*B_CCl(s,cos_theta,l,Mw_i)).conjugate()*M1rhok_M1rhopk_l(s,-0.5,-0.5,-0.5,l,cos_theta);
//    complex a19 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*C_CCl(mu,s,cos_theta,l,Mw_i)).conjugate()*M2rhok_M1rhopk_l(s,0.5,-0.5,-0.5,l,cos_theta);
//    complex a20 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*C_CCl(mu,s,cos_theta,l,Mw_i)).conjugate()*M2rhok_M1rhopk_l(s,-0.5,-0.5,-0.5,l,cos_theta);
//    complex a21 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*D_CCl(mu,s,cos_theta,l,Mw_i)).conjugate()*M3rhok_M1rhopk_l(s,0.5,0.5,-0.5,l,cos_theta);
//    complex a22 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*D_CCl(mu,s,cos_theta,l,Mw_i)).conjugate()*M3rhok_M1rhopk_l(s,-0.5,0.5,-0.5,l,cos_theta);
//    complex a23 = 1./s*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*E_CCl(mu,s,cos_theta,l,Mw_i)).conjugate()*M3rhok_M1rhopk_l(s,0.5,-0.5,-0.5,l,cos_theta);
//    complex a24 = 1./s*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y)*(A_CC(k,Mw_i)*E_CCl(mu,s,cos_theta,l,Mw_i)).conjugate()*M3rhok_M1rhopk_l(s,-0.5,-0.5,-0.5,l,cos_theta);
    
    Ml = 1./s/s*(M1rhok_M1rhopk_l(s,0.5,0.5,k,l,cos_theta).real()*Al(mu,l,0.5,k,s,Mw_i,W,X,Y).abs2()+
         M1rhok_M1rhopk_l(s,-0.5,-0.5,k,l,cos_theta).real()*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y).abs2() +
         2.*(M1rhok_M1rhopk_l(s,0.5,-0.5,k,l,cos_theta)*Al(mu,l,0.5,k,s,Mw_i,W,X,Y)*Al(mu,l,-0.5,k,s,Mw_i,W,X,Y).conjugate()).real());
//            +
//         2.*(a1).real()+
//         2.*(a2).real()+
//         2.*(a3).real()+
//         2.*(a4).real()+
//         2.*(a5).real()+
//         2.*(a6).real()+
//         2.*(a7).real()+
//         2.*(a8).real()+
//         2.*(a9).real()+
//         2.*(a10).real()+
//         2.*(a11).real()+
//         2.*(a12).real();
//         2.*(a13).real()+
//         2.*(a14).real()+
//         2.*(a15).real()+
//         2.*(a16).real()+ 
//         2.*(a17).real()+
//         2.*(a18).real()+      
//         2.*(a19).real()+
//         2.*(a20).real()+    
//         2.*(a21).real()+
//         2.*(a22).real()+
//         2.*(a23).real()+
//         2.*(a24).real();
            
    return (Ml);
    
}



double EWSMOneLoopLEP2::dsigmaBrem_l(const StandardModel::lepton l, const double k,
                                 const double s, const double Mw_i, const double cos_theta,
                      const double W, const double X, const double Y,const double GammaZ) const{
    
    double mf = SM.getLeptons(l).getMass();
    double mu = 10.;
    double beta = sqrt(1.-4.*mf*mf/s);
    double Qe = SM.getLeptons(1).getCharge();
    double Qf = SM.getLeptons(l).getCharge();
    double M1 = M1rhok_M1rhopk_l(s,0.5,0.5,k,l,cos_theta).real()+M1rhok_M1rhopk_l(s,-0.5,-0.5,k,l,cos_theta).real()
                +2.*(M1rhok_M1rhopk_l(s,-0.5,0.5,k,l,cos_theta)).real();
    double M2 = g_rhofl(l,0.5,Mw_i)*g_rhofl(l,0.5,Mw_i)*M1rhok_M1rhopk_l(s,0.5,0.5,k,l,cos_theta).real()
                +g_rhofl(l,-0.5,Mw_i)*g_rhofl(l,-0.5,Mw_i)*M1rhok_M1rhopk_l(s,-0.5,-0.5,k,l,cos_theta).real()
                +2.*g_rhofl(l,0.5,Mw_i)*g_rhofl(l,-0.5,Mw_i)*(M1rhok_M1rhopk_l(s,-0.5,0.5,k,l,cos_theta)).real();
    complex M3 = g_rhofl(l,0.5,Mw_i)*M1rhok_M1rhopk_l(s,0.5,0.5,k,l,cos_theta).real()
                +g_rhofl(l,-0.5,Mw_i)*M1rhok_M1rhopk_l(s,-0.5,-0.5,k,l,cos_theta).real()
                +g_rhofl(l,0.5,Mw_i)*M1rhok_M1rhopk_l(s,0.5,-0.5,k,l,cos_theta)+
                 g_rhofl(l,-0.5,Mw_i)*M1rhok_M1rhopk_l(s,-0.5,0.5,k,l,cos_theta);
    double MgammaZerosquared = Qe*Qe*Qf*Qf*Chi_gamma(mu,s,Mw_i,W,X,Y).abs2()/s/s*M1;
    double MZetaZerosqared = g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*Chi_Z(mu,s,Mw_i,W,X,Y).abs2()/s/s*M2;
    
    complex MgammaZ = Qe*Qf*g_rhoe(k,Mw_i)*Chi_gamma(mu,s,Mw_i,W,X,Y).conjugate()*Chi_Z(mu,s,Mw_i,W,X,Y)/s/s*M3;
    
    return (SM.getAle()*SM.getAle()/4./s*beta* (deltagammagamma_softl(s,l,cos_theta)*MgammaZerosquared
            + deltaZZ_softl(s,l,cos_theta,GammaZ)*MZetaZerosqared
            +2.*(deltagammaZ_softl(s,l,cos_theta,GammaZ)*MgammaZ).real()));
    
}


double EWSMOneLoopLEP2::dsigmaBrem_q(const QCD::quark q, const double k,
                                 const double s, const double Mw_i, const double cos_theta,
                      const double W, const double X, const double Y,const double GammaZ) const{
    
    double mf = SM.getQuarks(q).getMass();
    double mu = 10.;
    double beta = sqrt(1.-4.*mf*mf/s);
    double Qe = SM.getLeptons(1).getCharge();
    double Qf = SM.getQuarks(q).getCharge();
    double M1 = M1rhok_M1rhopk_q(s,0.5,0.5,k,q,cos_theta).real()+M1rhok_M1rhopk_q(s,-0.5,-0.5,k,q,cos_theta).real()
                +2.*(M1rhok_M1rhopk_q(s,-0.5,0.5,k,q,cos_theta)).real();
    double M2 = g_rhofq(q,0.5,Mw_i)*g_rhofq(q,0.5,Mw_i)*M1rhok_M1rhopk_q(s,0.5,0.5,k,q,cos_theta).real()
                +g_rhofq(q,-0.5,Mw_i)*g_rhofq(q,-0.5,Mw_i)*M1rhok_M1rhopk_q(s,-0.5,-0.5,k,q,cos_theta).real()
                +2.*g_rhofq(q,0.5,Mw_i)*g_rhofq(q,-0.5,Mw_i)*(M1rhok_M1rhopk_q(s,-0.5,0.5,k,q,cos_theta)).real();
    complex M3 = g_rhofq(q,0.5,Mw_i)*M1rhok_M1rhopk_q(s,0.5,0.5,k,q,cos_theta).real()
                +g_rhofq(q,-0.5,Mw_i)*M1rhok_M1rhopk_q(s,-0.5,-0.5,k,q,cos_theta).real()
                +g_rhofq(q,0.5,Mw_i)*M1rhok_M1rhopk_q(s,0.5,-0.5,k,q,cos_theta)+
                 g_rhofq(q,-0.5,Mw_i)*M1rhok_M1rhopk_q(s,-0.5,0.5,k,q,cos_theta);
    double MgammaZerosquared = Qe*Qe*Qf*Qf*Chi_gamma(mu,s,Mw_i,W,X,Y).abs2()/s/s*M1;
    double MZetaZerosquared = g_rhoe(k,Mw_i)*g_rhoe(k,Mw_i)*Chi_Z(mu,s,Mw_i,W,X,Y).abs2()/s/s*M2;
    
    complex MgammaZ = Qe*Qf*g_rhoe(k,Mw_i)*Chi_gamma(mu,s,Mw_i,W,X,Y).conjugate()*Chi_Z(mu,s,Mw_i,W,X,Y)/s/s*M3;
    
    return (SM.getAle()*SM.getAle()/4./s*beta* (deltagammagamma_softq(s,q,cos_theta)*MgammaZerosquared
            + deltaZZ_softq(s,q,cos_theta,GammaZ)*MZetaZerosquared
            +2.*(deltagammaZ_softq(s,q,cos_theta,GammaZ) * MgammaZ).real()));
    
}


double EWSMOneLoopLEP2::dsigma_l(const StandardModel::lepton l,const double s,
                                  const double Mw_i, const double cos_theta,
                               const double W,const double X,const double Y,const double GammaZ) const{
    
    double mf = SM.getLeptons(l).getMass();
    double beta = sqrt(1.-4.*mf*mf/s);
    //double beta = 1.;
    
    return (0.25*(SM.getAle()*SM.getAle()/4./s*beta*MTOTl_sq(l,0.5,s,Mw_i,cos_theta,W,X,Y,GammaZ)+
            //+dsigmaBrem_l(l,0.5,s,Mw_i,cos_theta,W,X,Y,GammaZ)+
            SM.getAle()*SM.getAle()/4./s*beta*MTOTl_sq(l,-0.5,s,Mw_i,cos_theta,W,X,Y,GammaZ)));
            //+dsigmaBrem_l(l,-0.5,s,Mw_i,cos_theta,W,X,Y,GammaZ));
}


double EWSMOneLoopLEP2::dsigma_q(const QCD::quark q,const double s,const double Mw_i,
                                   const double cos_theta,const double W,const double X,
                               const double Y,const double GammaZ) const{
    
    double mf = SM.getQuarks(q).getMass();
    double beta = sqrt(1.-4.*mf*mf/s);
    
    return (0.25*(SM.getAle()*SM.getAle()/4./s*beta*MTOTq_sq(q,0.5,s,Mw_i,cos_theta,W,X,Y,GammaZ)
            +dsigmaBrem_q(q,0.5,s,Mw_i,cos_theta,W,X,Y,GammaZ)+
            SM.getAle()*SM.getAle()/4./s*beta*MTOTq_sq(q,-0.5,s,Mw_i,cos_theta,W,X,Y,GammaZ)
            +dsigmaBrem_q(q,-0.5,s,Mw_i,cos_theta,W,X,Y,GammaZ)));
    
}
