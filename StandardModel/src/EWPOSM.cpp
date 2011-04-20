/*
 * File:   EWPOSM.cpp
 * Author: aleksandr
 *
 * Created on February 23, 2011, 4:33 PM
 */

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_zeta.h>
#include "EWPOSM.h"



EWPOSM::EWPOSM(StandardModel& sm) {
    me = sm.getLeptons(sm.ELECTRON).getMass();
    mmu = sm.getLeptons(sm.MU).getMass();
    mtau = sm.getLeptons(sm.TAU).getMass();
    
    mu = sm.getQuarks(sm.UP).getMass();
    md = sm.getQuarks(sm.DOWN).getMass();
    ms = sm.getQuarks(sm.STRANGE).getMass();
    mt = sm.getQuarks(sm.TOP).getMass();

    mc = sm.mcMz();
    mb = sm.mbMz();

    alsMz = sm.getAlsMz();
    GF = sm.getGF();
    ale = sm.getAle();
    dAletotmz = sm.dAleTotalMz();
    mZ = sm.getMZ();
    mHl = sm.getMHl();
}

EWPOSM::~EWPOSM() {
}


///////////////////////////////////////////////////////////////////////////

void EWPOSM::Check_string(const std::string ferm) {
    if( ferm!="up" && ferm!="down" && ferm!="chlepton" ) {
        std::cout << "wrong input for the string 'ferm', "
                  << "which must be either up, down or chlepton"<< std::endl;
        exit(EXIT_FAILURE);
    }
}

double EWPOSM::Qf(const std::string ferm) {
    Check_string(ferm);
    
    if(ferm=="up") {
        return (2.0/3.0);
    } else if(ferm=="down") {
        return (-1.0/3.0);
    } else if(ferm=="chlepton") {
        return (-1.0);
    } else {
        exit(EXIT_FAILURE);        
    }
}


///////////////////////////////////////////////////////////////////////////
// Analytic formulas for the W boson mass and the effective weak mixing angle

double EWPOSM::mW() const {
    // Eqs. (6), (7) and (9) in hep-ph/0311148
    // applicable for 100 GeV <= mHl <= 1 TeV

    if(mHl<100.||mHl>1000.)
    {
        std::cout << "Higgs mass out of range in mW()" << std::endl;
        exit(EXIT_FAILURE);
    }

    const double Mw0 = 80.3800;
    const double c1 = 0.05253;
    const double c2 = 0.010345;
    const double c3 = 0.001021;
    const double c4 = -0.000070;
    const double c5 = 1.077;
    const double c6 = 0.5270;
    const double c7 = 0.0698;
    const double c8 = 0.004055;
    const double c9 = 0.000110;
    const double c10 = 0.0716;
    const double c11 = 115.0;

    // mt, mZ, dAletotalmZ() and alsMz have to be varied within their combined
    // 2 sigma region around their central values (year 2003) adopted below.
    double dH = log(mHl/100.0);
    double dh = pow((mHl/100.0), 2.0);
    double dt = pow((mt/174.3), 2.0) - 1.0;
    double dZ = mZ/91.1875 - 1.0;
    double dalphae = dAletotmz/0.05907 - 1.0;
    double dalphas = alsMz/0.119 - 1.0;

    return (Mw0 - c1*dH - c2*dH*dH + c3*pow(dH, 4.0)
                + c4*(dh - 1.0) - c5*dalphae + c6*dt - c7*dt*dt
                - c8*dH*dt + c9*dh*dt - c10*dalphas + c11*dZ);
}

double EWPOSM::sw2() const {
    return ( 1.0 - mW()*mW()/mZ/mZ );
}

double EWPOSM::cw2() const { 
    return ( mW()*mW()/mZ/mZ );
}

double EWPOSM::sin2thwall(const std::string ferm) const {
    //Effective mixing angle for leptons,and c,b quarks
    //http://arXiv.org/abs/0811.1364v2, http://arXiv.org/abs/hep-ph/0608099v2
    // applicable for 10 GeV <= mHl <= 1 TeV

    //order is charged lepton, c, b quark
    const double s0[4] ={0.2312527,0.2311395,0.2327580,0.2310286};
    const double d1[4] = {4.729*0.0001,4.726*0.0001,4.749*0.0001,4.720*0.0001};
    const double d2[4] = {2.07*0.00001,2.07*0.00001,2.03*0.00001,2.06*0.00001};
    const double d3[4] = {3.85*0.000001,3.85*0.000001 ,3.94*0.000001,3.85*0.000001};
    const double d4[4] = {-1.85*0.000001,-1.85*0.000001,-1.84*0.000001,-1.85*0.000001};
    const double d5[4] = {0.0207,0.0207,0.0208,0.0207};
    const double d6[4] = {-0.002851,-0.002853,-9.93*0.0001,-0.002848};
    const double d7[4] = {1.82*0.0001,1.83*0.0001,7.08*0.00001,1.81*0.0001};
    const double d8[4] = {-9.74*0.000001,-9.73*0.000001,-7.61*0.000001,-9.73*0.000001};
    const double d9[4] = {3.98*0.0001,3.98*0.0001,4.03*0.0001,3.97*0.0001};
    const double d10[4] = {-0.655,-0.655,0.661,-0.655};

    double L_H = log(mHl/100.0);
    double Delta_H = mHl/100.0;
    double Delta_alphae = dAletotmz/0.05907 - 1.0;
    double Delta_t = pow((mt/178.0), 2.0) - 1.0;
    double Delta_alphas = alsMz/0.117 - 1.0;
    double Delta_Z = mZ/91.1876 - 1.0;
    int i=5;
    if(ferm=="chlepton") {i=0;}
    if(ferm=="c"||"up"){i=1;}
    if(ferm=="b"){i=2;}
    if(ferm=="down"){i=1;}

    double sin2t=sw2();
    if(i!=5) {
        sin2t= s0[i] + d1[i]*L_H + d2[i]*L_H*L_H + d3[i]*L_H*L_H*L_H*L_H
                + d4[i]*(Delta_H*Delta_H - 1.0) + d5[i]*Delta_alphae + d6[i]*Delta_t
                + d7[i]*Delta_t*Delta_t + d8[i]*Delta_t*(Delta_H - 1.0)
                + d9[i]*Delta_alphas + d10[i]*Delta_Z;
    }
    return sin2t;
}

double EWPOSM::vf(const std::string ferm) {
    Check_string(ferm);
    
    return ( 1.0-4.0*std::abs(Qf(ferm))*sin2thwall(ferm) );
}


///////////////////////////////////////////////////////////////////////////
// Loop functions in Bardin et al., Z.Phys.C 44,493-502 (1989)

//A.9
double EWPOSM::lambda(double Q2, double M12, double M22) {
    return ( (Q2+M12+M22)*(Q2+M12+M22)-4.0*M12*M22 );
}

//C.6
double EWPOSM::fcurve(double Q2, double M12, double M22) {
    //double j;
    double sql=lambda(Q2,M12,M22);
    sql=sqrt(sql);
    return 1.0/sql*log((Q2+M12+M22+sql)/(Q2+M12+M22-sql));
    //return j;
}

//C.6
gslpp::complex EWPOSM::fcurve_c(double Q2, double M12, double M22) {
    gslpp::complex sql(lambda(Q2,M12,M22));
    gslpp::complex xx(Q2/2.0/M12/M12);
    double m1=sqrt(M12);
    double m2=sqrt(M22);
    if((Q2!=0)|| (M12!=M22)) {
        if((std::abs(Q2)>=-(m1-m2)*(m1-m2))) {
            sql=sqrt(sql);
            if((M12!=M22) || std::abs(M12/Q2)>0.000001 ) {
                sql=1.0/sql*log((sql+Q2+M12+M22))-1.0/sql*log((-sql+Q2+M12+M22));
            } else {
                sql=1.0/sql*log((sql+Q2+M12+M22))-1.0/sql*log(xx);
            }
        }
    
        if((Q2<=-(m1-m2)*(m1-m2))&&(Q2>=-(M12+M22))) {
            sql=2.0/sqrt(-lambda(Q2,M12,M22))*atan(sqrt(-lambda(Q2,M12,M22))/(Q2+M12+M22));
        }

        if((Q2>-(m1+m2)*(m1+m2))&&Q2<=-(M12+M22)) { 
            sql=2.0/sqrt(-lambda(Q2,M12,M22))*(atan(sqrt(-lambda(Q2,M12,M22))/(Q2+M12+M22))+M_PI);
        }

    } else { 
        sql=1.0/M12;
    }

    return sql;
}

//C.1
gslpp::complex EWPOSM::fcurve2(double s, double Mv2) {
    return ( fcurve_c(-s,Mv2,Mv2)*s );
}

//C.3
double EWPOSM::lambdas(double s, double Mv2) {
    return ( s*s-4.0*Mv2*s );
}

//C.2
gslpp::complex EWPOSM::fcurve3(double s, double Mv2) {
    gslpp::complex II(0.0,1.0,false);
    gslpp::complex ls(-lambdas(s,Mv2),0.0,false);
    ls=sqrt(ls);
    gslpp::complex jj3=log((s+II*ls)/(s-II*ls));
    return ( jj3*jj3 );
}

//A.8
gslpp::complex EWPOSM::LL(double Q2,double M12, double M22) {
    return ( fcurve_c(Q2,M12,M22)*lambda(Q2,M12,M22) );
}

//D.1
gslpp::complex EWPOSM::I0(double Q2, double M12, double M22) {
    gslpp::complex l(LL(Q2,M12,M22));
    gslpp::complex i;
    if(Q2!=0.0) {
        i=0.5*log(M12*M22/mW()/mW()/mW()/mW())-2.0-(M12-M22)/2.0/Q2*log(M12/M22)
                +1.0/2.0/Q2*l;
    } else { 
        if(M12==M22) {
            i=log(M12/mW()/mW());
        } else { 
            i= -log(mW()*mW())+1.0/(M12-M22)*log(M12/M22); 
        }   
    }
    return i;
}

//D.2
gslpp::complex EWPOSM::I1(double Q2, double M12, double M22) {
    gslpp::complex l(LL(Q2,M12,M22));
    gslpp::complex i;
    if(M22!=0.0) {
        if(std::abs(M12+M22)<=std::abs(Q2)) {
            if(Q2!=0.0) { 
                i=0.25*log(M12*M22/mW()/mW()/mW()/mW())-1.0-(M12-M22)/2.0/Q2
                        -(2.0*M12/Q2+(M12-M22)*(M12-M22)/Q2/Q2)*log(M12/M22)
                        +(1.0+(M12-M22)/Q2)/4.0/Q2*l;
            } else {
                if(M12!=M22) {
                    i=0.25*(-1.0+2.0*M22/(M12-M22)-2.0*M22*M22*log(M12/M22)
                            /(M12-M22)/(M12-M22) + 2.0*log(M12/mW()/mW())); 
                } else {
                    i=0.5*log(M12/mW()/mW());
                }
            }
        } else { 
            if(M12!=M22) {
                i=(1.0/(120.0*pow((M12-M22),8.0)))*
                        (-30.0*(M12 - 3*M22)*pow((M12 - M22),7.0) +
                        20.0* pow((M12 - M22),5.0)*(M12*M12- 5.0*M12*M22 - 2.0*M22*M22)*Q2 - 
                        5.0* pow((M12 - M22),3.0)*(pow(M12,3.0) - 11.0*M12*M12*M22 - 47.0*M12*M22*M22 -
                        3.0* M22*M22*M22)*Q2*Q2 + 
                        2.0* (pow(M12,5.0) - 20*pow(M12,4.0)*M22 - 220.0*pow(M12,3.0)*M22*M22 + 80.0*M12*M12*pow(M22,3.0) + 
                        155.0*M12 *pow(M22,4.0) + 4*pow(M22,5.0))* pow(Q2,3.0) + 
                        60.0* (-pow((M12 - M22),4.0)*M22*M22*(M12*M12 + M22*M22 - 
                        2.0* M12 *(M22 + Q2)) *log(M12/M22) +
                        M12*M22*M22*Q2*Q2 *(-3.0* pow(M12,3.0) + 2.0*M22*M22 *(-M22 + Q2) +
                        4.0* pow(M12,2.0) *(M22 + Q2) + M12 *M22 *(M22 + 8.0* Q2))* log(M12)  +
                        M12*M22*M22*Q2*Q2*((M12 - M22)*(M12-M22)* (3.0 *M12 + 2.0* M22) -
                        2.0* (2.0*M12*M12 + 4.0*M12*M22 + M22*M22)*Q2)*log(M22)
                        +pow((M12-M22),8.0)*log(M12/(mW()*mW())) )
                        );
            } else {
                i=(70.0*M12*M12*Q2 - 7.0*M12*Q2*Q2 + Q2*Q2*Q2 
                        + 420.0*M12*M12*M12*log(M12/mW()/mW()))/(840.0 *M12*M12*M12);
            }
        }
    } else {
        if(Q2!=0.0) {
            //auxiliary variable
            gslpp::complex aa(M12/(M12+Q2));
            i=0.5/Q2/Q2*(-Q2*(M12+2*Q2)-
                    (M12+Q2)*(M12+Q2)*log(aa)
                    +Q2*Q2*log(M12/mW()/mW())
                    );
        } else{
            i=-0.25+0.5*log(M12/mW()/mW());
        }
    }
    return i;
}

//D.3
gslpp::complex EWPOSM::I3(double Q2, double M12, double M22) {
    gslpp::complex i;
    double M12x=0.01;
    double M22x=0.01;
    
    if(M12>0.1) M12x=M12;
    if(M22>0.1) M22x=M22;
    gslpp::complex l(LL(Q2,M12x,M22x));
    if(Q2!=0.0) { 
        i=log(M12x*M22x/mW()/mW()/mW()/mW())/12.0 -5.0/18.0
                +(M12x+M22x)/3.0/Q2+(M12x-M22x)*(M12x-M22x)/3.0/Q2/Q2
                +((M12x*M12x-M22x*M22x)/Q2/Q2/4.0+(M12x-M22x)*(M12x-M22x)*(M12x-M22x)/6.0/Q2/Q2/Q2)*log(M12x/M22x)
                +(0.5-(M12x+M22x)/Q2/2.0-(M12x-M22x)*(M12x-M22x)/Q2/Q2)/6.0/Q2*l;
    } else i=0.0;//because I3(0.) is always multiplied by 0!!!
    /*else{if(std::abs((M12x-M22x)/M22x)>0.5){i=-1.0/36.0/(M12x-M22x)/(M12x-M22x)/(M12x-M22x)*
                    (6.0*(3*M12x*M22x*M22x-M22x*M22x*M22x)*log(M12x/M22x))-(M12x-M22x)*
                    (-5.0*(M12x*M12x+22.0*M12x*M22x-5*M22x*M22x+6.0*(M12x-M22x)*(M12x-M22x)*log(M12x/mW()/mW())));}
                 else i=1.0/6.0*log(M12x/mW()/mW())+(M22x-M12x)/12.0/M12x-(M22x-M12x)*(M22x-M12x)/40.0/M12x/M12x;}*/
    return i;
}


///////////////////////////////////////////////////////////////////////////
// Loop functions copied from the ZFitter codes

gslpp::complex EWPOSM::DL(double Q2, double Q2SBT,double M12, double M22) { 
    // another copy of ZFitter code
    //    DL(Q2,AMQ2,AM12,AM22)=(L(Q2,AM12,AM22)-L(Q2SBT,AM12,AM22))/(Q2-Q2SBT)
     double EPS=0.001;
     double Q2S=Q2SBT+M12+M22;
     double ALAM=Q2S*Q2S-4.0*M12*M22;
     double DSLAM=sqrt(std::abs(ALAM));
     double QD=Q2-Q2SBT;
     double RQD=std::abs(QD/DSLAM);
     gslpp::complex XDL;
     //double RR=4.0*M12/M22;
     if(RQD>=EPS) {
         XDL=(LL(Q2,M12,M22)-LL(Q2SBT,M12,M22))/QD;
     } else {
         double RR=4.0*M12/M22;
         if(RR >1.0) {
             gslpp::complex XJS=fcurve_c(Q2SBT,M12,M22);
             XDL=2.0+Q2S*XJS+QD/ALAM*(Q2S-2.0*M12*M22*XJS)
                     +(QD/ALAM)*(QD/ALAM)*(-Q2S*Q2S/3.0-8.0/3.0*M12*M22*Q2S*XJS);
         } else {
             double RAT=QD/M22;
             XDL=4.0+2.0/3.0*RAT-2.0/15.0*RAT*RAT;
         }
     }
     return XDL;
}

gslpp::complex EWPOSM::DI0(double Q2, double MV2, double M12, double M22) {    
    //double RAT=std::abs(M12/MV2);
    double AL12=log(M12/M22);
    return ( (LL(Q2,M12,M22)-(M12-M22)*AL12)/2.0/Q2
            -DL(Q2,-MV2,M12,M22)/2.0 );
}

gslpp::complex EWPOSM::DI3(double Q2, double AMZ2, double M12, double M22) {
    double EPS=0.0001;
    double AL12=log(M12/M22);
    double DM12=(M12-M22)/Q2;
    double SM12=(M12+M22)/Q2;
    double AQ=Q2/AMZ2;
    double SMV1=1.0-AQ;
    double SMV2=1.0-AQ+AQ*AQ;
    gslpp::complex XDI3(0.0);
    if(M22>0.00001) { 
        XDI3=SM12/3.0+DM12*DM12/3.0*SMV1
                +(SM12*DM12/4.0*SMV1+DM12*DM12*DM12/6.0*SMV2)*AL12
                +(0.50-SM12/2.0*SMV1-DM12*DM12*SMV2)/6.0/Q2*LL(Q2,M12,M22)
                -(0.50+0.50*(M12+M22)/AMZ2
                -((M12-M22)/AMZ2)*(M12-M22)/AMZ2/AMZ2)/6.0*DL(Q2,-AMZ2,M12,M22);
    } else {
    // CHAIN1 WILL BE USED ONLY FOR W-WIDTH
        double QV=(Q2+AMZ2)/AMZ2;
        if(std::abs(QV)>EPS) { 
            XDI3=(I3(Q2,M12,M22)-I3(-AMZ2,M12,M22))/QV;
        }
        double R1V=M12/AMZ2;
        double VQ=AMZ2/Q2;
        double VQ2=1.0-VQ+VQ*VQ;
        AQ=M12/Q2;
        double RDI3=-1.0/6.0+(VQ-0.50)/3.0*R1V+VQ2/3.0*R1V*R1V
                    -QV*(0.50+R1V)/6.0
                    -VQ*R1V*R1V*((VQ-1.0)/2.0+VQ2/3.0*R1V)*log(std::abs(1.0+1.0/AQ));
        XDI3=RDI3;
    }
    return XDI3;
}


///////////////////////////////////////////////////////////////////////////
// Self energies of the W and Z bosons 

//Eq. 70 of Gfitter
double EWPOSM::sigmaBF0() {
    double rw=mHl*mHl/mW()/mW();
    return ( 5.0*cw2()*(1.0+cw2())/8.0-17.0/4.0+5.0/8.0/cw2()-rw/8.0
            +(9.0/4.0+3.0/4.0/cw2()-3/sw2())*log(cw2())
            +3.0*rw/4.0/(1.0-rw)*log(rw) );
}

//Eq. 71
gslpp::complex EWPOSM::sigmaBFMw() {
    gslpp::complex lwh(LL(-(mW()*mW()),(mW()*mW()),(mHl*mHl)));
    gslpp::complex lwz(LL(-(mW()*mW()),(mW()*mW()),(mZ*mZ)));
    double rw=mHl*mHl/mW()/mW();
    lwh=lwh*(0.5-rw/6.0+rw*rw/24.0)/(mW()*mW());
    lwz=lwz*(-2.0*cw2()-17.0/6.0+2.0/3.0/cw2()+1.0/24.0/cw2()/cw2())/(mW()*mW());
    //total formula
    lwh=lwh+lwz+(-157.0/9.0+23.0/12.0/cw2()+1.0/12.0/cw2()/cw2()-rw/2.0
            +rw*rw/12.0+1.0/cw2()*(-7.0/2.0+7.0/12.0/cw2()+1.0/24.0/cw2()/cw2())*log(cw2())
            +rw*(-3.0/4.0+rw/4.0-rw*rw/24.0)*log(rw));
    return lwh;
}

//Eq. 72
gslpp::complex EWPOSM::sigmaBFMz() {
    //gsl_complex j;
    double rz=mHl*mHl/mZ/mZ;
    double rw=mHl*mHl/mW()/mW();

    gslpp::complex lzh(LL(-(mZ*mZ),(mZ*mZ),(mHl*mHl)));
    gslpp::complex lww(LL(-(mZ*mZ),(mW()*mW()),(mW()*mW())));
    lzh=lzh*(0.5-rz/6.0+rz*rz/24.0)/mW()/mW();
    lww=lww*(-2.0*cw2()*cw2()*cw2()-17.0/6.0*cw2()*cw2()+2.0/3.0*cw2()+1.0/24.0)/mW()/mW();
    //total formula
    lww=lww+lzh-8.0*cw2()*cw2()-34.0*cw2()/3.0+35.0/18.0*(1.0+1.0/cw2())-rw/2.0
            +rz*rz/12.0/cw2()+rw*(-3.0/4.0+rz/4.0-rz*rz/24.0)*log(rz)
            +5.0/6.0/cw2()*log(cw2());
    return lww;
}

//Eq. 73
gslpp::complex EWPOSM::sigmaBFMz1() {
    double rz=mHl*mHl/mZ/mZ;
    double rw=mHl*mHl/mW()/mW();
    gslpp::complex lzh(LL(-(mZ*mZ),(mZ*mZ),(mHl*mHl)));
    gslpp::complex lww(LL(-(mZ*mZ),(mW()*mW()),(mW()*mW())));
    lzh=lzh*(0.5-5.0*rz/24.0+rz*rz/12.0+1.0/2.0/(rz-4.0))/mW()/mW();
    lww=lww*(-cw2()*cw2()*cw2()+7.0/6.0*cw2()*cw2()-17.0/12.0*cw2()-1.0/8.0)/mW()/mW();
    //total formula
    lww=lww+lzh-4*cw2()*cw2()+17.0/3.0*cw2()-23.0/9.0+
            5.0/18.0/cw2()-rw/2.0+rw*rz/6.0
            +rw*(-3.0/4.0+3.0*rz/8.0-rz*rz/12.0)*log(rz)-1.0/12.0/cw2()*log(cw2())
            +log(rz)/2.0/cw2();
    return lww;
}

//Eq. 74
// at the end we will have to write it as a sum ,
//(261) of hep-ph 9709229
// but right now I will write a code only for the first generation;
gslpp::complex EWPOSM::sigmaFww(double s) {
    double muq[3]={mu,mc,mt};
    double mdq[3]={md,ms,mb};
    double ml[3]={me,mmu,mtau};
    gslpp::complex i3;
    gslpp::complex i1u;
    gslpp::complex i1d;
    gslpp::complex i3l;
    gslpp::complex i1cl;
    gslpp::complex i1nu;
    double RMLW;
    double ALLW,ALTW,ALBW;
    double RMTW;
    double RMBW;
    //not clear but I guess s=(mW()*mW()) in this case
    gslpp::complex sig(0.0);
    if(s!=0.0) {
        for(int i=0;i<3;i++) {
            i3=I3(-(s),muq[i]*muq[i],mdq[i]*mdq[i]);
            i1u=I1(-(s),muq[i]*muq[i],mdq[i]*mdq[i]);
            i1d=I1(-(s),mdq[i]*mdq[i],muq[i]*muq[i]);
            //Neutrino lepton contribution
            i3l=I3(-(s),ml[i]*ml[i],0.0);
            i1cl=I1(-(s),ml[i]*ml[i],0.0);
            sig=sig+3.0*(-2.0*s*i3/mW()/mW()+i1u/(mW()*mW())*(muq[i]*muq[i])
                    +i1d/(mW()*mW())*(mdq[i]*mdq[i]))
                    //lepton term
                    -2.0*s*i3l/mW()/mW()+i1cl/(mW()*mW())*(ml[i]*ml[i])
                    ;
        }
    } else {
        //actually you do not need this else your function works anyway...code is just copied from ZFitter code        
        for(int i=0;i<3;i++) {
            RMLW=ml[i]*ml[i]/mW()/mW();
            RMTW=muq[i]*muq[i]/mW()/mW();
            RMBW=mdq[i]*mdq[i]/mW()/mW();
            ALLW=log(RMLW);
            ALTW=log(RMTW);
            ALBW=log(RMBW);
            sig=sig+0.5*(RMLW*ALLW-RMLW/2.0);
            if(RMTW!=RMBW){
                sig=sig+3.0/2.0*((RMTW*RMTW*ALTW-RMBW*RMBW*ALBW)/(RMTW-RMBW)
                        -(RMTW+RMBW)/2.0);
            } else {
                sig=sig+3.0*(RMTW*ALTW-RMTW/2.0);
            }
        }
    }
    return sig;
}

//Eq. 75
// at the end we will have to write it as a sum ,
// but right now I will write a code only for the first generation;
gslpp::complex EWPOSM::sigmaFzz() {
    //gsl_complex j;
    //not clear but I guess s=(mZ*mZ) in this case
    double s=mZ*mZ;
    const double muq[3]={mu,mc,mt};
    const double mdq[3]={md,ms,mb};
    const double ml[3]={me,mmu,mtau};
    gslpp::complex i3u;
    gslpp::complex i0u;
    gslpp::complex i3d;
    gslpp::complex i0d;
    gslpp::complex i3l;
    gslpp::complex i0l;
    gslpp::complex sig(0.0);
    for(int i=0;i<3;i++) {
        gslpp::complex i3u=(I3(-s,muq[i]*muq[i],muq[i]*muq[i]));
        gslpp::complex i0u=(I0(-s,muq[i]*muq[i],muq[i]*muq[i]));
        gslpp::complex i3d=(I3(-s,mdq[i]*mdq[i],mdq[i]*mdq[i]));
        gslpp::complex i0d=(I0(-s,mdq[i]*mdq[i],mdq[i]*mdq[i]));
        gslpp::complex i3l=(I3(-s,ml[i]*ml[i],ml[i]*ml[i]));
        gslpp::complex i0l=(I0(-s,ml[i]*ml[i],ml[i]*ml[i]));

        sig=sig+1.0/2.0/cw2()*(3.0*(-s/mZ/mZ*(1.0+vf("up")*vf("up"))*i3u+i0u*muq[i]*muq[i]/(mZ*mZ))
                +3.0*(-s/mZ/mZ*(1.0+vf("down")*vf("down"))*i3d+i0d*mdq[i]*mdq[i]/(mZ*mZ))
                +(-s/mZ/mZ*(1.0+vf("chlepton")*vf("chlepton"))*i3l+i0l*ml[i]*ml[i]/(mZ*mZ)));
        //std::cout<<(1.0+vf("chlepton")*vf("chlepton"))<<"\n";
    }
    //neutrino part from ZFitter
    gslpp::complex XALR(2.0*log(mW()/mZ),M_PI,false);
    sig=sig+3.0/cw2()/6.0*(XALR+5.0/3.0);
    return sig;
}

//Eq. 76
// at the end we will have to write it as a sum ,
// but right now I will write a code only for the first generation;
gslpp::complex EWPOSM::sigmaFzz1() {
    //not clear but I guess s=(mW()*mW()) in this case
    gslpp::complex sig(0.0);
    double ru;
    double rd;
    double rl;

    gslpp::complex fu;
    gslpp::complex fd;
    gslpp::complex fl;
    const double muq[3]={mu,mc,mt};
    const double mdq[3]={md,ms,mb};
    const double ml[3]={me,mmu,mtau};

    double NG=3.0;
    double ALQ=0.0;
    double AIL1=-M_PI;
    double AIL2=0.0;
    double ALR=log(cw2());
    gslpp::complex XL1(ALQ,AIL1,false);
    gslpp::complex XL2(ALQ,AIL2,false);
    //double QD=0.0;
    double RQD=0.0;
    gslpp::complex XLQ=XL1+1.0+RQD/2.0+RQD*RQD/3.0;
    gslpp::complex XSNU=NG/cw2()/6.0*(-5.0/3.0-ALR+XLQ);
    
    
    for(int i=0;i<3;i++) {
        fu=fcurve_c(-(mZ*mZ),muq[i]*muq[i],muq[i]*muq[i]);
        fd=fcurve_c(-(mZ*mZ),mdq[i]*mdq[i],mdq[i]*mdq[i]);
        fl=fcurve_c(-(mZ*mZ),ml[i]*ml[i],ml[i]*ml[i]);
        ru=muq[i]*muq[i]/(mW()*mW());
        rd=mdq[i]*mdq[i]/(mW()*mW());
        rl=ml[i]*ml[i]/(mW()*mW());
        
        sig= sig
                -0.0*(ru/2.0*(1.0-ru*(mW()*mW())*fu)+1.0/6.0/cw2()*(1+vf("up")*vf("up"))*
                (0.5*log(ru*cw2())+ru*cw2()+(-1.0/4.0/cw2()+ru/2.0-ru*ru*cw2())*(mW()*mW())*fu))
                
                -0.0*(rd/2.0*(1.0-rd*(mW()*mW())*fd)+1.0/6.0/cw2()*(1+vf("down")*vf("down"))*
                (0.5*log(rd*cw2())+rd*cw2()+(-1.0/4.0/cw2()+rd/2.0-rd*rd*cw2())*(mW()*mW())*fd))
                
                -LEPTON_TMP*(rl/2.0*(1.0-rl*(mW()*mW())*fl)+1.0/6.0/cw2()*(1+vf("chlepton")*vf("chlepton"))*
                (0.5*log(rl*cw2())+rl*cw2()+(-1.0/4.0/cw2()+rl/2.0-rl*rl*cw2())*(mW()*mW())*fl))
                ;
    }
    sig=sig+0.0*XSNU;
    return sig;
}

gslpp::complex EWPOSM::sigmaFzz1a() {
    /*I will try to translate ZFitter code to c, the formulas they are using are different from   
     *the ones present in the papers, and I have no idea where they come from...*/
    /*Q2-s always equal to -mZ^2*/

    double AML2[3]={me*me,mmu*mmu,mtau*mtau};
    double AMu2[3]={mu*mu,mc*mc,mt*mt};
    double AMd2[3]={md*md,ms*ms,mb*mb};
    
    double NG=3.0;
    double ALQ=0.0;
    double AIL1=-M_PI;
    double AIL2=0.0;
    double ALR=log(cw2());
    gslpp::complex XL1(ALQ,AIL1,false);
    gslpp::complex XL2(ALQ,AIL2,false);
    //double QD=0.0;
    double RQD=0.0;
    gslpp::complex XLQ=XL1+1.0+RQD/2.0+RQD*RQD/3.0;
    gslpp::complex XSNU=NG/cw2()/6.0*(-5.0/3.0-ALR+XLQ);
    gslpp::complex XSI3(0.0,0.0,false);
    gslpp::complex XSQMI3(0.0,0.0,false);
    gslpp::complex XSQ2I3(0.0,0.0,false);
    gslpp::complex XSDI0(0.0,0.0,false);
    //LEPTONS
    for(int ii=0;ii<3;ii++) {
        XSI3=XSI3+I3(-mZ*mZ,AML2[ii],AML2[ii])-DI3(-mZ*mZ,mZ*mZ,AML2[ii],AML2[ii]);
        XSDI0=XSDI0+AML2[ii]/mW()/mW()*DI0(-mZ*mZ,mZ*mZ,AML2[ii],AML2[ii]);
    }
      XSQMI3=XSI3;
      XSQ2I3=XSI3;
      for(int I=0;I<3;I++) {
          //up quarks
          XSI3=XSI3+3.0*(I3(-mZ*mZ,AMu2[I],AMu2[I])-DI3(-mZ*mZ,mZ*mZ,AMu2[I],AMu2[I]));
          XSQMI3=XSQMI3+3.0*(2.0/3.0)*(I3(-mZ*mZ,AMu2[I],AMu2[I])
                  -DI3(-mZ*mZ,mZ*mZ,AMu2[I],AMu2[I]));
          XSQ2I3=XSQ2I3+3.0*(2.0/3.0)*(2.0/3.0)*(I3(-mZ*mZ,AMu2[I],AMu2[I])
                  -DI3(-mZ*mZ,mZ*mZ,AMu2[I],AMu2[I]));
          
          XSDI0=XSDI0+3.0*AMu2[I]/mW()/mW()*DI0(-mZ*mZ,mZ*mZ,AMu2[I],AMu2[I]);
          // std::cout<<DI0(-mZ*mZ,mZ*mZ,AMu2[I],AMu2[I])<<"\n";
          
          //std::cout<<XSDI0<<"\n";
          //down quarks
          XSI3=XSI3+3.0*(I3(-mZ*mZ,AMd2[I],AMd2[I])-DI3(-mZ*mZ,mZ*mZ,AMd2[I],AMd2[I]));
          
          XSQMI3=XSQMI3+3.0*(1.0/3.0)*(I3(-mZ*mZ,AMd2[I],AMd2[I])
                  -DI3(-mZ*mZ,mZ*mZ,AMd2[I],AMd2[I]));
          XSQ2I3=XSQ2I3+3.0*(1.0/3.0)*(1.0/3.0)*(I3(-mZ*mZ,AMd2[I],AMd2[I])
                  -DI3(-mZ*mZ,mZ*mZ,AMd2[I],AMd2[I]));
          XSDI0=XSDI0+3.0*AMd2[I]/mW()/mW()*DI0(-mZ*mZ,mZ*mZ,AMd2[I],AMd2[I]);
          // std::cout<<XSDI0<<"\n";
      }
      
      return ( (8.0*sw2()*sw2()*XSQ2I3-4.0*sw2()*XSQMI3+XSI3)/cw2()+XSDI0/2.0+XSNU );
}

//Eq (63) of Gfitter 
gslpp::complex EWPOSM::delta_rho_W_F() {
    gslpp::complex d(sigmaBFMw());
    d=-d.real()+sigmaBF0()-sigmaFww(mW()*mW()).real()+sigmaFww(0.0);
    return d;
}

//unnumbered equation after (65)
gslpp::complex EWPOSM::delta_rho_Z_F() {
    gslpp::complex d(sigmaBFMz());
    d=-d+sigmaBF0()-sigmaFzz()+sigmaFww(0.0);
    return d;
}

//Eq A.3.41 of ZFitter
gslpp::complex EWPOSM::Delta_rho() {
    gslpp::complex x;    
    x=sigmaBFMw()+sigmaFww(mW()*mW())-sigmaFzz()-sigmaBFMz();
    return x;
}


///////////////////////////////////////////////////////////////////////////
// Z-gamma mixing 

//Eq. 68 this is only fermionic contribution
gslpp::complex EWPOSM::piZg() {
    gslpp::complex pizg(I3(-(mZ*mZ),mu*mu,mu*mu));
    
    pizg=(pizg+I3(-(mZ*mZ),mc*mc,mc*mc)+I3(-(mZ*mZ),mt*mt,mt*mt))*2.0*3.0*std::abs(Qf("up"))*vf("up");
    //std::cout<<pizg<<"\n";
    pizg=pizg+2.0*3.0*std::abs(Qf("down"))*vf("down")*( I3(-(mZ*mZ),md*md,md*md)+I3(-(mZ*mZ),ms*ms,ms*ms)+I3(-(mZ*mZ),mb*mb,mb*mb));
    //std::cout<<pizg<<"\n";
    pizg=pizg+LEPTON_TMP*2.0*std::abs(Qf("chlepton"))*vf("chlepton")*(I3(-(mZ*mZ),me*me,me*me)+I3(-(mZ*mZ),mmu*mmu,mmu*mmu)+I3(-(mZ*mZ),mtau*mtau,mtau*mtau));
    return pizg;
}

//I guess this is something like -\Pi_{Z\gamma} bosonic contribution I do not know exactly it is copied from ZFitter ...
gslpp::complex EWPOSM::XAMM1() {
    gslpp::complex XL4=LL(-mZ*mZ,mW()*mW(),mW()*mW())/mW()/mW();
    gslpp::complex X=2.0/9.0/cw2()+35.0/18.0-34.0/3.0*cw2()-8.0*cw2()*cw2()
            +(1.0/24.0+2.0/3.0*cw2()-17.0/6.0*cw2()*cw2()-2.0*cw2()*cw2()*cw2())*XL4;
    return X;
}


///////////////////////////////////////////////////////////////////////////
// Vertex functions

//Eq 66(for simplicity I will write it down for the (mW()*mW())...)
//I still do not understand why to put zero for the decay width...
gslpp::complex EWPOSM::v1W(double s) {

    const double GammaW0=2.085; // ?????
    
    int i;
    double Rw=mW()*mW()/s;
    gsl_sf_result re,im;
    gslpp::complex Rwt(Rw,-mW()*GammaW0*0.0/s,false);
    gslpp::complex rwt=1.0/Rwt+1.0;
    double z=rwt.abs();
    double arg=rwt.arg();
    i=gsl_sf_complex_dilog_e (z,arg,&re,&im);
    gslpp::complex lterm(re.val-1.0/6.0*M_PI*M_PI,M_PI*log(1+1/Rw),false);
    //std::cout<<-log(-Rwt)*(3.0+2.0*Rw)<<"\n";
    //std::cout<<2.0*(1+Rw)*(1+Rw)*lterm<<"\n";
    gslpp::complex vv=-log(-Rwt)*(3.0+2.0*Rw)-2.0*Rw-7.0/2.0+
            2.0*(1+Rw)*(1+Rw)*lterm;
    return vv;
}

//Eq 66(for simplicity I will write it down for the (mZ*mZ)...)
//I do not understand this code I just copied it from Z FItter it does not agree with the formulas they give in the paper...
gslpp::complex EWPOSM::v1Z(double s) {

    const double GammaZ0=2.4952; // ?????
    
    int i;
    double Rz=mZ*mZ/s;
    gsl_sf_result re,im;
    gslpp::complex Rzt(Rz,-0.0*mZ*GammaZ0/s,false);
    gslpp::complex rzt=1.0/Rzt+1.0;
    double z=rzt.abs();
    double arg=rzt.arg();
    i=gsl_sf_complex_dilog_e(z,arg,&re,&im);
    gslpp::complex lterm(re.val-1.0/6.0*M_PI*M_PI,M_PI*log(1.0+1.0/Rz),false);
    gslpp::complex vv=-log(-Rzt)*(3.0+2.0*Rz)-2.0*Rz-7.0/2.0+
            2.0*(1+Rz)*(1+Rz)*lterm;
    return vv;
}

//g++ -L/usr/local/lib Zdec.o -lgsl -lgslcblas -lm
//g++ -Wall -I/usr/local/include -c Zdec.cpp

//Eq 67
gslpp::complex EWPOSM::v2W(double s) {
    double Rw=(mW()*mW())/s;
    gslpp::complex vv(fcurve3(s,(mW()*mW())));
    gslpp::complex lww(LL(-s,(mW()*mW()),(mW()*mW())));
    vv=vv*2*Rw*(Rw+2.0)-(7.0/6.0+Rw)*lww/s-1.0/6.0-2.0*Rw;
    return vv;
}

//Eq 67, this is an alternative copy pasted from ZFitter code, I do not understand it...
gslpp::complex EWPOSM::v2Wa() {
    double SR=sqrt(4.0*cw2()-1.0);
    double AT=atan(SR/(2.0*cw2()-1.0));
    gslpp::complex V2ZWW=2.0/9.0/cw2()/cw2()+43.0/18.0/cw2()-1.0/6.0-2.0*cw2()
            +(-1.0/12.0/cw2()/cw2()-1.50/cw2()+7.0/3.0+2.0*cw2())*SR*AT
            -2.0*cw2()*(2.0+cw2())*AT*AT;
    return V2ZWW;
}
		
//unnumbered formula after (65)
gslpp::complex EWPOSM::uf(const std::string ferm) {
    Check_string(ferm);
    
    gslpp::complex uf(0.0);
    double qf = Qf(ferm);
    gslpp::complex v1z(v1Z(mZ*mZ));
    gslpp::complex v1w(v1W(mZ*mZ));
    gslpp::complex v2w(v2Wa());
    uf=1.0/4.0/cw2()*(1.0-6.0*std::abs(qf)*sw2()+12.0*qf*qf*sw2()*sw2())*(v1z)
            +(0.5-cw2()-std::abs(qf)*sw2())*(v1w)+cw2()*(v2w);
    // std::cout<<v1Z(mZ*mZ)<<"\n";
    //std::cout<<v1W(mZ*mZ)<<"\n";
    //std::cout<<v2W(mZ*mZ)<<"\n";
    return uf;
}


///////////////////////////////////////////////////////////////////////////
// delta rho

//Eq. 2.3.15 of ZFitter
gslpp::complex EWPOSM::delta_rho() {
    double SCALE = ale/4.0/M_PI/sw2()*(41.0/6.0-11.0/3.0*cw2())*log(cw2());    
    // std::cout<<"step 1"<<sw2()/cw2()*SCALE<<"\n";
    // std::cout<<"step 2"<<ale/4.0/M_PI/sw2()*Deltarho()<<"\n";
    return ( ale/4.0/M_PI/sw2()*Delta_rho()+sw2()/cw2()*SCALE );
}

double EWPOSM::delta_rho_QCD() {
    //const double zeta2=1.64493;
    double zeta2 = gsl_sf_zeta_int(2);
    //double ut=mt*mt/(mZ*mZ);
    //double ub=mb*mb/(mZ*mZ);
    return ( alsMz/M_PI*(-1.0/4.0)*(zeta2+0.5)*ale/M_PI*mt*mt/sw2()/(mW()*mW()) );
}


///////////////////////////////////////////////////////////////////////////
// delta r

// Eq (16) of ZFitter
double EWPOSM::delta_r() const {
    return ( 1.0-GF*mZ*mZ/sqrt(8.0)/M_PI/ale
             *(4.0*mW()*mW()/mZ/mZ-4.0*mW()*mW()*mW()*mW()/mZ/mZ/mZ/mZ) );
}


///////////////////////////////////////////////////////////////////////////
// Reminder "delta rho_rem"

//Eq. (64) of Gfitter
//Eq. A.4.70 of ZFitter
gslpp::complex EWPOSM::delta_rho_rem(const std::string ferm) {
    Check_string(ferm);
    
    gslpp::complex d(0.0);
    d=sigmaFzz1a();
    d=d+sigmaBFMz1()-delta_rho_Z_F()-11.0/2.0+5.0/8.0*cw2()*(1+cw2())
            -9.0/4.0*cw2()/sw2()*log(cw2())+2.0*uf(ferm);
    d=ale/4.0/M_PI/sw2()*d;
    return d;
}


///////////////////////////////////////////////////////////////////////////
// Reminder "delta kappa_rem"

//Eq (65) of Gfitter
//Eq A.4.71 of ZFitter
gslpp::complex EWPOSM::delta_k_rem(const std::string ferm) {
    Check_string(ferm);
    
    gslpp::complex k(0.0);
    double qf = Qf(ferm);
    //delta_rho_Z_F()-delta_rho_W_F() is total delta rho^F 
    k = delta_rho_Z_F()-delta_rho_W_F();
    gslpp::complex v1z(v1Z(mZ*mZ));
    
    // why here is k.real I do not know I just copied it from ZFitter code ...
    k=-cw2()/sw2()*k.real()-piZg()+XAMM1()+sw2()*sw2()/cw2()*qf*qf*v1z-uf(ferm);
    //std::cout<<-piZg()<<"\n";
    //std::cout<<"v1z"<<v1z*sw2()*sw2()/cw2()*qf*qf<<"\n";
    //std::cout<<"uf"<<uf(ferm)<<"\n";
    //std::cout<<k<<"\n";
    k=ale/4.0/M_PI/sw2()*k;
    return k;
}


///////////////////////////////////////////////////////////////////////////
// Reminder "Delta r_rem"

// EW corrections
//   Eq 62 of Gfitter
//   Eq 2.4.23 of ZFitter
gslpp::complex EWPOSM::Delta_r_rem() {
    gslpp::complex d(delta_rho_W_F());
    d=( 1.0/sw2()*(d+11.0/2.0-5.0/8.0*cw2()*(1.0+cw2())+9.0/4.0*cw2()/sw2()*log(cw2()))
            -2.0/3.0+1.0/sw2()*(1.0/6.0*24.0-1.0/6.0-7.0*cw2())*log(cw2()) )
            *sqrt(2.0)*GF*mZ*mZ*sw2()*cw2()/4.0/M_PI/M_PI;
    return d;
}

// F_1(x) in Eq. (30) of B.A.Kniehl, NPB347, 86 (1990)
double EWPOSM::F1(double x) {
    double b=log(1.0-x);
    //const double c1=0.2933;
    //const double c2=-0.0192;
    //const double c3=0.4417;
    //const double c4=-0.2837;
    const double f1=-1.1881;
    const double f2=-2.0979;
    const double f3=4.1157;
    const double f4=-2.2082;
    const double f5=3.6968;
    const double f6=-2.1815;
    //const double zeta2=1.64493;
    double zeta2 = gsl_sf_zeta_int(2);

    return ( (1.0-x)*(1.0-x)*b*(3.0/8.0*b-zeta2-9.0/8.0)+
            f1+f2*x+f3*x*x+f4*x*x*x+f5*x*x*x*x+f6*x*x*x*x*x );
}

// F_1^inf(x) below Eq. (32) of B.A.Kniehl, NPB347, 86 (1990)
double EWPOSM::F1inf(double x) {
    const double f1=-1.1881;
    const double f2=-2.0979;
    //const double zeta2=1.64493;
    double zeta2 = gsl_sf_zeta_int(2);
    
    return ( f1+(f2+zeta2+9.0/8.0)*x );
}

// EW+QCD corrections
double EWPOSM::Delta_r_rem_QCD() {  
    const double c1=0.2933;
    const double c2=-0.0192;
    const double c3=0.4417;
    //const double c4=-0.2837;
    //const double f1=-1.1881;
    //const double f2=-2.0979;
    //const double f3=4.1157;
    //const double f4=-2.2082;
    //const double f5=3.6968;
    //const double f6=-2.1815;
    //const double zeta2=1.64493;
    double zeta2 = gsl_sf_zeta_int(2);
    double ut=mt*mt/(mZ*mZ);
    double ub=mb*mb/(mZ*mZ);
    return ( alsMz/M_PI*(-cw2()/sw2()*(-1.0/4.0)*(zeta2+0.5)*ale/M_PI*mt*mt/sw2()/(mW()*mW())
            +ale/M_PI*(1.0/2.0/sw2()*(-1.0/4.0/sw2()+1.0/3.0)*log(ut)
            -1.0/9.0*log(ub)+c1/sw2()/sw2()+c2/sw2()+c3 +1.0/sw2()*(1.0/sw2()-2.0)*
            (F1(cw2()/ut)-F1inf(cw2()/ut)))) );
}

// Total corrections
//   Eq. 61 of Gfitter 
gslpp::complex EWPOSM::Delta_r_rem_tot() {
    return ( Delta_r_rem() + Delta_r_rem_QCD() );
}


///////////////////////////////////////////////////////////////////////////
// Electroweak form factors

//Eq. 59 in Gfitter, arXiv:0811.0009
gslpp::complex EWPOSM::rhoZf(const std::string ferm) {
    Check_string(ferm);
    
    gslpp::complex r(0.0);
    r=delta_rho_rem(ferm);
    gslpp::complex drmtot(Delta_r_rem_tot());
    r=(1.0+r)/(1.0+delta_rho()*(1.0-drmtot));
    return r;
}

//Eq.60 in Gfitter, arXiv:0811.0009
gslpp::complex EWPOSM::kZf(const std::string ferm) {
    Check_string(ferm);
    
    gslpp::complex k(0.0);
    k=delta_k_rem(ferm);
    gslpp::complex drho(delta_rho());
    gslpp::complex drrtot(Delta_r_rem_tot());
    //std::cout<<k<<" "<<drho<<"\n";
    k=(1.0+k)*(1.0-cw2()/sw2()*drho*(1.0-drrtot));
    return k;
}

//extracting k from sin^2theta_eff
double EWPOSM::rekZf(const std::string ferm) {
    Check_string(ferm);
    
    return ( sin2thwall(ferm)/sw2() );
}

//Eq 2.4.12 of ZFitter, hep-ph/9908433v3
gslpp::complex EWPOSM::gzf(const std::string ferm) { 
    Check_string(ferm);
    
    gslpp::complex k(kZf(ferm));
    double qf = Qf(ferm);
    k=1.0-4.0*std::abs(qf)*(k*sw2()+ale*ale*35.0/18.0*(1.0-8.0/3.0*k.real()*sw2()));      
    return k;
}

//Eq 2.4.12 of ZFitter, hep-ph/9908433v3
gslpp::complex EWPOSM::gzf1(const std::string ferm) { 
    Check_string(ferm);
    
    double qf = Qf(ferm);
    gslpp::complex xx( (1.0-4.0*std::abs(qf)*sin2thwall(ferm)),
                        -4.0*std::abs(qf)*kZf(ferm).imag()*sw2(), false );
    return xx;
}


///////////////////////////////////////////////////////////////////////////

//some functions from ZFitter code
gslpp::complex EWPOSM::spence(double X) {
    //double F1 = 1.64493406684822618;    

    return 0.0;
}
      
gslpp::complex EWPOSM::fspence(double X) {
    return 0.0;
}


// ROQCD Eq. 3.4.5 of ZFitter paper only light quarks
double EWPOSM::ROQCD() {
    // for simplicity I set everywhere s=M_z^2
    return ( -ale*alsMz/M_PI/M_PI/16.0/cw2()/sw2()
            *(2+vf("up")*vf("up")+vf("down")*vf("down"))*(-1.0) );
}

// AKQCD Eq. 3.4.5 of ZFitter paper only light quarks
double EWPOSM::AKQCD() {
    return ( ale*alsMz/M_PI/M_PI/4.0/sw2()/sw2()*(cw2()*log(cw2())) );
}

// \delta r_{ud} from Eq 3.4.4 of ZFitter paper only light quarks
double EWPOSM::CLQQCD() {
    return ( -ale*alsMz/M_PI/M_PI/4.0*(cw2()-sw2())/sw2()/sw2()*log(cw2()) );
}





