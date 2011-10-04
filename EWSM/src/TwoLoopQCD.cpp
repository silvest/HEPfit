/* 
 * File:   TwoLoopQCD.cpp
 * Author: mishima
 */

#include <gsl/gsl_sf.h>
#include "TwoLoopQCD.h"


TwoLoopQCD::TwoLoopQCD(const EWSMcommon& EWSMC_i) : EWSMC(EWSMC_i) {
}

//TwoLoopQCD::TwoLoopQCD(const TwoLoopQCD& orig) {
//}

TwoLoopQCD::~TwoLoopQCD() {
}


////////////////////////////////////////////////////////////////////////

double TwoLoopQCD::DeltaAlpha_l() const {
    return (0.0);
}    

double TwoLoopQCD::DeltaAlpha_t() const {   
    double xt = pow(EWSMC.GetSM().getMz()
                /EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass(), 2.0);
    double tmp = (5.062 + xt*0.8315)*EWSMC.GetSM().getAlsMz()/M_PI;
    tmp *= -4.0/45.0*EWSMC.GetSM().getAle()/M_PI*xt;
    return tmp;
} 

double TwoLoopQCD::DeltaRho() const {
    return ( 3.0*EWSMC.GetXt_alpha()*EWSMC.GetAlsMt()/M_PI*deltaQCD_2() );     
}

double TwoLoopQCD::DeltaR_rem() const {
    return ( (2.0*DeltaR_ud() + DeltaR_tb())
              + EWSMC.GetCW2()/EWSMC.GetSW2()/EWSMC.GetF_AlphaToGF()*DeltaRho() );     
}

complex TwoLoopQCD::deltaRho_rem_l(const StandardModel::lepton l) const {
    return ( (2.0*DeltaRho_ud() + DeltaRho_tb()) - DeltaRho() );       
}

complex TwoLoopQCD::deltaRho_rem_q(const StandardModel::quark q) const {
    return ( (2.0*DeltaRho_ud() + DeltaRho_tb()) - DeltaRho() );    
}

complex TwoLoopQCD::deltaKappa_rem_l(const StandardModel::lepton l) const {
    return ( (2.0*DeltaKappa_ud() + DeltaKappa_tb())
              - EWSMC.GetCW2()/EWSMC.GetSW2()*DeltaRho() );  
}

complex TwoLoopQCD::deltaKappa_rem_q(const StandardModel::quark q) const {
    return ( (2.0*DeltaKappa_ud() + DeltaKappa_tb())
              - EWSMC.GetCW2()/EWSMC.GetSW2()*DeltaRho() );    
}


////////////////////////////////////////////////////////////////////////

double TwoLoopQCD::deltaQCD_2() const {
    return ( -2.0/3.0*(1.0+2.0*EWSMC.GetZeta2()) );
}

double TwoLoopQCD::F1(const double x) const {
    if (x < 0.0 || x >= 1.0) throw "x is out of range in TwoLoopQCD::F1";    
       
    /* Zeta functions */
    double zeta_2 = EWSMC.GetZeta2();
    double zeta_3 = EWSMC.GetZeta3();    

    if (x == 0.0) return (23.0/16.0 - zeta_2/2.0 - 3.0/2.0*zeta_3);     
    
    /* Dilogarithm and Trilogarithm */
    double Li2_x, Li3_x, Li3_mx_1mx;
    double Mw = EWSMC.GetMw();
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass(); 
    if (x == Mw*Mw/Mt/Mt) {
        Li2_x = EWSMC.GetLi2_MW2toMTOP2();
        Li3_x = EWSMC.GetLi3_MW2toMTOP2();
        Li3_mx_1mx = EWSMC.GetLi3_for_F1();
    } else { 
        Li2_x = gsl_sf_dilog(x);
        Polylogarithms* myPolyLog;
        myPolyLog = new Polylogarithms();
        Li3_x = myPolyLog->Li3(x);
        Li3_mx_1mx = myPolyLog->Li3(-x/(1.0 - x));
        delete myPolyLog;
    }
    
    double b = log(1.0 - x);    

    double F1;
    F1 = (x - 3.0/2.0 + 1.0/2.0/x/x)
           *(- Li3_x - Li3_mx_1mx + b/3.0*(2.0*Li2_x - zeta_2) + b*b*b/6.0)
         + 1.0/3.0*(x + 1.0/2.0 - 1.0/2.0/x)*Li2_x 
         + b*b/6.0*(x - 3.0/4.0 - 3.0/2.0/x + 5.0/4.0/x/x)
         - b/4.0*(x - 5.0/2.0 + 2.0/3.0/x + 5.0/6.0/x/x)
         + zeta_3*(x - 3.0/2.0) + zeta_2/3.0*(x - 7.0/4.0 - 1.0/2.0/x)
         + 13.0/12.0 - 5.0/24.0/x;
    return F1;
}

double TwoLoopQCD::V1(const double r) const {
    if (r < 0.0 || r >= 1.0) throw "r is out of range in TwoLoopQCD::V1";

    /* Zeta functions */
    double zeta_3 = EWSMC.GetZeta3(); 
    
    if (r == 0.0) return (0.0); 
    
    double Mz = EWSMC.GetSM().getMz(); 
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass(); 

    /* Logarithms etc */
    double Phi, gamma, h;
    if (r==Mz*Mz/4.0/Mt/Mt) {
        Phi = EWSMC.GetPhi_QCD2();
        gamma = EWSMC.GetGamma_QCD2();
        h = EWSMC.GetH_QCD2();
    } else {
        Phi= asin(sqrt(r));
        gamma = log(2.0*sqrt(r));
        h = log(2.0*sqrt(1.0-r));
    }
        
    /* Clausen functions */
    double Cl3_2Phi, Cl3_4Phi, Cl2_2Phi, Cl2_4Phi;
    if (r == Mz*Mz/4.0/Mt/Mt) {
        Cl3_2Phi = EWSMC.GetCl3_2Phi();
        Cl3_4Phi = EWSMC.GetCl3_4Phi();
        Cl2_2Phi = EWSMC.GetCl2_2Phi();
        Cl2_4Phi = EWSMC.GetCl2_4Phi();
    } else {
        ClausenFunctions* myClausen;
        myClausen = new ClausenFunctions();        
        Cl3_2Phi = myClausen->Cl3(2.0*Phi);
        Cl3_4Phi = myClausen->Cl3(4.0*Phi);    
        Cl2_2Phi = myClausen->Cl2(2.0*Phi); 
        Cl2_4Phi = myClausen->Cl2(4.0*Phi); 
        delete myClausen;
    }

    double V1;    
    V1 = 4.0*(r - 1.0/4.0/r)
           *(2.0*Cl3_2Phi - Cl3_4Phi  + 8.0/3.0*Phi*(Cl2_2Phi - Cl2_4Phi) 
             - 4.0/3.0*Phi*Phi*(gamma + 2.0*h))
         + sqrt(1.0/r - 1.0)*(8.0/3.0*(r + 1.0/2.0)
                               *(- Cl2_2Phi + Cl2_4Phi + 2.0*Phi*(gamma + 2.0*h)) 
                              - 2.0*Phi*(r + 3.0/2.0))
         + 8.0*Phi*Phi*(r - 1.0/6.0 - 7.0/48.0/r) + 13.0/6.0 + zeta_3/r;
    return V1;
}

double TwoLoopQCD::A1(const double r) const {
    if (r < 0.0 || r >= 1.0) throw "r is out of range in TwoLoopQCD::A1";
 
    /* Zeta functions */
    double zeta_2 = EWSMC.GetZeta2();
    double zeta_3 = EWSMC.GetZeta3(); 
        
    if (r == 0.0) return (3.0*(7.0/4.0 - zeta_2 - 2.0*zeta_3));         
     
    double Mz = EWSMC.GetSM().getMz();   
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();     

    /* Logarithms etc */
    double Phi, gamma, h;
    if (r==Mz*Mz/4.0/Mt/Mt) {
        Phi = EWSMC.GetPhi_QCD2();
        gamma = EWSMC.GetGamma_QCD2();
        h = EWSMC.GetH_QCD2();
    } else {
        Phi= asin(sqrt(r));
        gamma = log(2.0*sqrt(r));
        h = log(2.0*sqrt(1.0-r));
    }
    
    /* Clausen functions */
    double Cl3_2Phi, Cl3_4Phi, Cl2_2Phi, Cl2_4Phi;
    if (r == Mz*Mz/4.0/Mt/Mt) {
        Cl3_2Phi = EWSMC.GetCl3_2Phi();
        Cl3_4Phi = EWSMC.GetCl3_4Phi();
        Cl2_2Phi = EWSMC.GetCl2_2Phi();
        Cl2_4Phi = EWSMC.GetCl2_4Phi();
    } else {
        ClausenFunctions* myClausen;
        myClausen = new ClausenFunctions();        
        Cl3_2Phi = myClausen->Cl3(2.0*Phi);
        Cl3_4Phi = myClausen->Cl3(4.0*Phi);    
        Cl2_2Phi = myClausen->Cl2(2.0*Phi); 
        Cl2_4Phi = myClausen->Cl2(4.0*Phi); 
        delete myClausen;
    }   
    
    double A1;    
    A1 = 4.0*(r - 3.0/2.0 + 1.0/2.0/r)*(2.0*Cl3_2Phi - Cl3_4Phi 
              + 8.0/3.0*Phi*(Cl2_2Phi - Cl2_4Phi) - 4.0/3.0*Phi*Phi*(gamma + 2.0*h))
         + sqrt(1.0/r - 1.0)*(8.0/3.0*(r - 1.0)
                                *(- Cl2_2Phi + Cl2_4Phi + 2.0*Phi*(gamma + 2.0*h)) 
                              - 2.0*Phi*(r - 3.0 + 1.0/4.0/r))
         + 8.0*Phi*Phi*(r - 11.0/12.0 + 5.0/48.0/r + 1.0/32.0/r/r) 
         - 3.0*zeta_2 + 13.0/6.0 + (-2.0*zeta_3 + 1.0/4.0)/r;    
    return A1;
}

double TwoLoopQCD::V1prime(const double r) const {
    if (fabs(r) >= 1.0) throw "r is out of range in TwoLoopQCD::V1prime()";
    /* Zeta functions */
    double zeta_3 = EWSMC.GetZeta3();     
    return (4.0*zeta_3 - 5.0/6.0 + 656.0/81.0*r );
}

double TwoLoopQCD::A1prime(const double r) const {
    if (fabs(r) >= 1.0) throw "r is out of range in TwoLoopQCD::A1prime()";
    /* Zeta functions */
    double zeta_3 = EWSMC.GetZeta3();     
    return (4.0*zeta_3 - 49.0/18.0 + 1378.0/405.0*r );    
}

double TwoLoopQCD::DeltaR_ud() const {
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    
    /* Logarithm */
    double log_cW2 = EWSMC.GetLog_cW2();     
    
    double DeltaR;
    DeltaR = - log_cW2;
    DeltaR *= (cW2 - sW2)/4.0/sW2/sW2;
    DeltaR *= EWSMC.GetSM().getAle()*EWSMC.GetSM().getAlsMz()/M_PI/M_PI;
    return DeltaR;   
}

double TwoLoopQCD::DeltaR_tb() const {
    double Mw = EWSMC.GetMw();
    double Mz = EWSMC.GetSM().getMz();  
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();
    double wt = Mt*Mt/Mw/Mw;
    double zt = Mt*Mt/Mz/Mz;
    double rZ4t = Mz*Mz/4.0/Mt/Mt;
    double xWt = Mw*Mw/Mt/Mt;
    
    double vt = EWSMC.vf(EWSMC.GetSM().TOP);
    double at = EWSMC.af(EWSMC.GetSM().TOP);
    double vb = EWSMC.vf(EWSMC.GetSM().BOTTOM);
    double ab = EWSMC.af(EWSMC.GetSM().BOTTOM);
    
    /* Zeta functions */
    double zeta_2 = EWSMC.GetZeta2();
    
    /* Logarithm */
    double log_zt = - 2.0*EWSMC.GetLogMZtoMTOP();
    
    double DeltaR;
    DeltaR = pow(EWSMC.Qf(EWSMC.GetSM().TOP), 2.0)*V1prime(0.0)
             + cW2/sW2/sW2*wt/4.0*(zeta_2 + 1.0/2.0)
             - zt/sW2/sW2*( vt*vt*V1(rZ4t) + at*at*(A1(rZ4t) - A1(0.0)) )
             + (cW2 - sW2)/sW2/sW2*wt*(F1(xWt) - F1(0.0))
             - vb*ab/2.0/sW2/sW2*log_zt;
    DeltaR *= EWSMC.GetSM().getAle()*EWSMC.GetAlsMt()/M_PI/M_PI;
    return DeltaR;  
}

double TwoLoopQCD::DeltaRho_ud() const {
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    
    double DeltaRho;
    DeltaRho = pow(EWSMC.vf(EWSMC.GetSM().UP), 2.0) 
               + pow(EWSMC.vf(EWSMC.GetSM().DOWN), 2.0)
               + pow(EWSMC.af(EWSMC.GetSM().UP), 2.0) 
               + pow(EWSMC.af(EWSMC.GetSM().DOWN), 2.0); 
    DeltaRho /= 4.0*sW2*cW2;    
    DeltaRho *= EWSMC.GetSM().getAle()*EWSMC.GetSM().getAlsMz()/M_PI/M_PI;
    return DeltaRho;      
}

double TwoLoopQCD::DeltaRho_tb() const {
    double Mz = EWSMC.GetSM().getMz();  
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();
    double zt = Mt*Mt/Mz/Mz;
    double rZ4t = Mz*Mz/4.0/Mt/Mt;
    
    double vt = EWSMC.vf(EWSMC.GetSM().TOP);
    double at = EWSMC.af(EWSMC.GetSM().TOP);
    
    double DeltaRho;
    DeltaRho = - (vt*vt*V1prime(rZ4t) + at*at*A1prime(rZ4t))
               + 4.0*zt*(vt*vt*V1(rZ4t) + at*at*A1(rZ4t))
               - 4.0*zt*F1(0.0);
    DeltaRho /= 4.0*sW2*cW2;
    DeltaRho *= EWSMC.GetSM().getAle()*EWSMC.GetAlsMt()/M_PI/M_PI;
    return DeltaRho;   
}

complex TwoLoopQCD::DeltaKappa_ud() const {
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    
    /* Logarithm */
    double log_cW2 = EWSMC.GetLog_cW2();     
        
    complex DeltaKappa(0.0,0.0,false);
    DeltaKappa = cW2/4.0/sW2/sW2*log_cW2 
                 + M_PI/4.0/sW2*(1.0 - 20.0/9.0*sW2)*(complex::i());
    DeltaKappa *= EWSMC.GetSM().getAle()*EWSMC.GetSM().getAlsMz()/M_PI/M_PI;
    return DeltaKappa;     
}

complex TwoLoopQCD::DeltaKappa_tb() const {
    double Mw = EWSMC.GetMw();
    double Mz = EWSMC.GetSM().getMz();  
    double sW2 = EWSMC.GetSW2();
    double cW2 = EWSMC.GetCW2();
    double Mt = EWSMC.GetSM().getQuarks(EWSMC.GetSM().TOP).getMass();
    double wt = Mt*Mt/Mw/Mw;
    double zt = Mt*Mt/Mz/Mz;
    double rZ4t = Mz*Mz/4.0/Mt/Mt;
    double xWt = Mw*Mw/Mt/Mt;
    
    double vt = EWSMC.vf(EWSMC.GetSM().TOP);
    double at = EWSMC.af(EWSMC.GetSM().TOP);
    double Qt = EWSMC.Qf(EWSMC.GetSM().TOP);
    double vb = EWSMC.vf(EWSMC.GetSM().BOTTOM);
    double ab = EWSMC.af(EWSMC.GetSM().BOTTOM);
    double Qb = EWSMC.Qf(EWSMC.GetSM().BOTTOM);
    
    /* Logarithm */
    double log_zt = - 2.0*EWSMC.GetLogMZtoMTOP();    
    
    complex DeltaKappa(0.0,0.0,false);
    DeltaKappa = 4.0*cW2*wt*(vt*vt*V1(rZ4t) + at*at*A1(rZ4t) - F1(xWt))
                 + 4.0*sW2*(fabs(Qt) - 4.0*sW2*Qt*Qt)*zt*V1(rZ4t)
                 + (vb*vb + ab*ab + sW2*(fabs(Qb) - 4.0*sW2*Qb*Qb))*log_zt;
    DeltaKappa += M_PI*sW2*(1.0/3.0 - 4.0/9.0*sW2)*(complex::i());
    DeltaKappa /= 4.0*sW2*sW2;
    DeltaKappa *= EWSMC.GetSM().getAle()*EWSMC.GetAlsMt()/M_PI/M_PI;
    return DeltaKappa;   
}




