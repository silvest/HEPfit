/*
 * File:   EWPOSM.cpp
 * Author: aleksandr
 *
 * Created on February 23, 2011, 4:33 PM
 */

#include <gslpp_complex.h>
#include <gslpp_vector_double.h>
#include <gslpp_vector_complex.h>
#include <gslpp_matrix_double.h>
#include <gslpp_matrix_complex.h>
#include <iostream>
#include <math.h>
#include <TF1.h>
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "QCD.h"
#include "StandardModel.h"
#include "EWPOSM.h"
using namespace gslpp;
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

    // mt, mZ, dAle5Mz and alsMz have to be varied within their combined
    // 2 sigma region around their central values (year 2003) adopted below.
    double dH = log(mHl/100.0);
    double dh = pow((mHl/100.0), 2.0);
    double dt = pow((mt/174.3), 2.0) - 1.0;
    double dZ = mZ/91.1875 - 1.0;
    double dalphae = dAle5Mz/0.05907 - 1.0;
    double dalphas = alsMz/0.119 - 1.0;

    return (Mw0 - c1*dH - c2*dH*dH + c3*pow(dH, 4.0)
                + c4*(dh - 1.0) - c5*dalphae + c6*dt - c7*dt*dt
                - c8*dH*dt + c9*dh*dt - c10*dalphas + c11*dZ);
}

double EWPOSM::sin2thwall(const std::string& ferm) const {
    //Effective mixing angle for leptons,and c,b quarks
    //http://arXiv.org/abs/0811.1364v2, http://arXiv.org/abs/hep-ph/0608099v2
    // applicable for 10 GeV <= mHl <= 1 TeV

    //order is charged lepton, c, b quark
    const double s0[3] ={0.2312527,0.2311395,0.2327580};
    const double d1[3] = {4.729*0.0001,4.726*0.0001,4.749*0.0001};
    const double d2[3] = {2.07*0.00001,2.07*0.00001,2.03*0.00001};
    const double d3[3] = {3.85*0.000001,3.85*0.000001 ,3.94*0.000001};
    const double d4[3] = {-1.85*0.000001,-1.85*0.000001,-1.84*0.000001};
    const double d5[3] = {0.0207,0.0207,0.0208};
    const double d6[3] = {-0.002851,-0.002853,-9.93*0.0001};
    const double d7[3] = {1.82*0.0001,1.83*0.0001,7.08*0.00001};
    const double d8[3] = {-9.74*0.000001,-9.73*0.000001,-7.61*0.000001};
    const double d9[3] = {3.98*0.0001,3.98*0.0001,4.03*0.0001};
    const double d10[3] = {-0.655,-0.655,0.661};

    double L_H = log(mHl/100.0);
    double Delta_H = mHl/100.0;
    double Delta_alphae = dAle5Mz/0.05907 - 1.0;
    double Delta_t = pow((mt/178.0), 2.0) - 1.0;
    double Delta_alphas = alsMz/0.117 - 1.0;
    double Delta_Z = mZ/91.1876 - 1.0;
    int i=4;
    if(ferm=="l") {i=0;}
    if(ferm=="c"){i=1;}
    if(ferm=="b"){i=2;}


    double sin2t;
    if(i!=4) {sin2t= s0[i] + d1[i]*L_H + d2[i]*L_H*L_H + d3[i]*L_H*L_H*L_H*L_H
                 + d4[i]*(Delta_H*Delta_H - 1.0) + d5[i]*Delta_alphae + d6[i]*Delta_t
                 + d7[i]*Delta_t*Delta_t + d8[i]*Delta_t*(Delta_H - 1.0)
                 + d9[i]*Delta_alphas + d10[i]*Delta_Z;}
    return sin2t;
}

//here are some functions from Bardin et al  Z.Phys.C 44,493-502 (1989)

//A.9
double EWPOSM::lambda(double Q2, double M12, double M22 )
{return (Q2+M12+M22)*(Q2+M12+M22)-4*M12*M22;
	}

//C.6
double EWPOSM::fcurve(double Q2, double M12, double M22)
{double j;
double sql=lambda(Q2,M12,M22);
sql=sqrt(sql);
 return 1.0/sql*log((Q2+M12+M22+sql)/(Q2+M12+M22-sql));
return j ;}

//C.6
gsl_complex EWPOSM::fcurve_c(double Q2, double M12, double M22)
{gsl_complex j;
complex sql(lambda(Q2,M12,M22),0,0);
complex xx(Q2/2.0/M12/M12,0, false);
sql=sqrt(sql);
if((M12!=M22) || abs(M12/Q2)>0.000001 ) {sql=1.0/sql*log((sql+Q2+M12+M22))-1.0/sql*log((-sql+Q2+M12+M22));}
else
{sql=1.0/sql*log((sql+Q2+M12+M22))-1.0/sql*log(xx);}
GSL_SET_COMPLEX(&j,sql.real(),sql.imag());
return j ;}


//C.1
gsl_complex EWPOSM::fcurve2(double s,double Mv2)
{gsl_complex j2;
	j2=gsl_complex_mul_real(fcurve_c(-s,Mv2,Mv2),s);
	return j2;
	}

// C.
double EWPOSM::lambdas(double s, double Mv2)
{
return s*s-4*Mv2*s;}


//C.2
gsl_complex EWPOSM::fcurve3(double s, double Mv2)
{gsl_complex j3;
complex II(0.0,1.0,false);
complex ls(-lambdas(s,Mv2),0.0,false);
ls=sqrt(ls);
complex jj3=log((s+II*ls)/(s-II*ls));
jj3=jj3*jj3;
GSL_SET_COMPLEX(&j3, jj3.real(), jj3.imag());
return j3;
	}
//A.8

gsl_complex EWPOSM::LL(double Q2,double M12, double M22)
{
	return gsl_complex_mul_real(fcurve_c(Q2,M12,M22),lambda(Q2,M12,M22));
	}

//D.1
gsl_complex EWPOSM::I0(double Q2, double M12, double M22)
{gsl_complex i0;
complex l(LL(Q2,M12,M22));
complex i=0.5*log(M12*M22/mW()/mW())-2.0-(M12-M22)/2/Q2*log(M12/M22)
+1/2/Q2*l;
GSL_SET_COMPLEX(&i0,i.real(),i.imag());
return i0;
}



//D.2
gsl_complex EWPOSM::I1(double Q2, double M12, double M22)
{gsl_complex i1;
complex l(LL(Q2,M12,M22));
complex i=0.25*log(M12*M22/mW()/mW())-1.0-(M12-M22)/2/Q2
-(2*M12/Q2+(M12-M22)*(M12-M22)/Q2/Q2)*log(M12/M22)
+(1.0+(M12-M22)/Q2)/4.0/Q2*l;
GSL_SET_COMPLEX(&i1,i.real(),i.imag());
return i;
}


//D.3
gsl_complex EWPOSM::I3(double Q2, double M12, double M22)
{gsl_complex i1;
complex l(LL(Q2,M12,M22));
complex i=log(M12*M22/mW()/mW())/12.0 -5.0/18.0
+(M12+M22)/3.0/Q2+(M12-M22)*(M12-M22)/3.0/Q2/Q2+
((M12*M12-M22*M22)/Q2/Q2/4.0+(M12-M22)*(M12-M22)*(M12-M22)/6.0/Q2/Q2/Q2)*log(M12/M22)
+(0.5-(M12+M22)/Q2/2.0-(M12-M22)*(M12-M22)/Q2/Q2)/6.0/Q2*l;
GSL_SET_COMPLEX(&i1,i.real(),i.imag());
return i;
}






//Equation 62 of GFitter
gsl_complex EWPOSM::deltarrem()
{gsl_complex dr;
complex d(deltarhowF());
d=(1.0/sw2*(d+11.0/2.0-5.0/8.0*cw2*(1+cw2)+9.0/4.0*cw2/sw2*log(cw2))
-2.0/3.0+1.0/sw2*(1/6*3-1/6-7*cw2)*log(cw2))*sqrt(2)*GF*mZ*mZ*sw2*cw2/4.0/M_PI/M_PI;
GSL_SET_COMPLEX(&dr,d.real(),d.imag());
	return dr;}



//Equation (63)
gsl_complex EWPOSM::deltarhowF()
{gsl_complex dr;
complex d(sigmaBFMw());
d=-d+sigmaBF0()-sigmaFFMw();
GSL_SET_COMPLEX(&dr,d.real(),d.imag());
return dr;}

//Equation (64)
gsl_complex EWPOSM::deltarhorem(const std::string& ferm)
{gsl_complex dr;
if((ferm=="u")||(ferm=="d")||ferm=="l")
	{
	complex d(sigmaFFMz1());
        complex uff(uf(ferm));
	d=d+sigmaBFMz1()-deltarhozF()-11.0/2.0+5.0/8.0*cw2*(1+cw2)-9.0/4.0*cw2/sw2*log(cw2)+2*uff;
	d=ale/4.0/M_PI/sw2;
	GSL_SET_COMPLEX(&dr,d.real(),d.imag());}

else {GSL_SET_COMPLEX(&dr,0.0,0.0);std::cout<<"wrong input for fermion parameter,  it can be equal only to 'u','d', 'l' deltarhorem"<<"\n";}
return dr;}


//Equation (65)
gsl_complex EWPOSM::deltakrem(const std::string& ferm)
{
gsl_complex kr;
double qf;
if((ferm=="u")||(ferm=="d")||ferm=="l")
{
	if(ferm=="u") qf=Qf[0];
	if(ferm=="d") qf=Qf[1];
	if(ferm=="l") qf=Qf[2];
	complex k(deltarhozF());
        complex v1z(v1Z(mZ*mZ));
	k=-cw2/sw2*k+piZg()+sw2*sw2/cw2*qf*qf*v1z-uf(ferm);
	k=ale/4.0/M_PI/sw2;
	GSL_SET_COMPLEX(&kr,k.real(),k.imag());}
else {GSL_SET_COMPLEX(&kr,0.0,0.0);
		std::cout<<"wrong input for fermion parameter "<< ferm<<"   it can be equal only to 'u','d', 'l'   deltakrem"<<"\n";}
return kr;}

//unnumbered formula after (65)
gsl_complex EWPOSM::uf(const std::string& ferm)
{
gsl_complex u;
double qf;
if((ferm=="u")||(ferm=="d")||(ferm=="l"))
{
	if(ferm=="u") qf=Qf[0];
	if(ferm=="d") qf=Qf[1];
	if(ferm=="l") qf=Qf[2];
	complex v1z(v1Z(mZ*mZ));
	complex v1w(v1W(mZ*mZ));
	complex v2w(v2W(mZ*mZ));
	complex uf=1/4.0/cw2*(1-6.0*abs(qf)*sw2+12*qf*qf*sw2*sw2)*v1z
	+(0.5-cw2-abs(qf)*sw2)*v1w+cw2*v2w;
	GSL_SET_COMPLEX(&u,uf.real(),uf.imag());
}
else {
GSL_SET_COMPLEX(&u,0.0,0.0);std::cout<<"wrong input for fermion parameter,  it can be equal only to 'u','d', 'l'  uf "<<"\n";};
return u;
	}

//unnumbered equation after (65)
gsl_complex EWPOSM::deltarhozF()
{gsl_complex dr;
complex d(sigmaBFMz());
d=-d+sigmaBF0()-sigmaFFMz();
GSL_SET_COMPLEX(&dr,d.real(),d.imag());
return dr;}



//Equation 66(for simplicity I will write it down for the (mW()*mW())...)
gsl_complex EWPOSM::v1W(double s)
{int i;
gsl_complex v;
double Rw=mW()*mW()/s;
gsl_sf_result re,im;
complex Rwt(Rw,-mW()*GammaW()/s,false);
complex rwt=Rwt+1;
double z=rwt.abs();
double arg=rwt.arg();
i=gsl_sf_complex_dilog_e (z,arg,&re,&im);
complex lterm(re.val-gsl_sf_dilog(1.0),im.val,false);
complex vv=-log(-Rwt)*(3.0+2.0*Rw)-2.0*Rw-7.0/2.0+
2.0*(1+Rw)*(1+Rw)*lterm;
GSL_SET_COMPLEX(&v,vv.real(),vv.imag());
	return v;
		}
//Equation 66(for simplicity I will write it down for the (mZ*mZ)...)
gsl_complex EWPOSM::v1Z(double s)
{int i;
gsl_complex v;
double Rz=mZ*mZ/s;
gsl_sf_result re,im;
complex Rzt(Rz,-mZ*GammaZ()/s,false);
complex rzt=Rzt+1;
double z=rzt.abs();
double arg=rzt.arg();
i=gsl_sf_complex_dilog_e(z,arg,&re,&im);
complex lterm(re.val-gsl_sf_dilog(1.0),im.val,false);
complex vv=-log(-Rzt)*(3.0+2.0*Rz)-2.0*Rz-7.0/2.0+
2.0*(1+Rz)*(1+Rz)*lterm;
GSL_SET_COMPLEX(&v,vv.real(),vv.imag());
	return v;
		}
//g++ -L/usr/local/lib Zdec.o -lgsl -lgslcblas -lm
//g++ -Wall -I/usr/local/include -c Zdec.cpp

//Equation 67
gsl_complex EWPOSM::v2W(double s)
{gsl_complex v;
double Rw=(mW()*mW())/s;
complex vv(fcurve3(s,(mW()*mW())));
complex lww(LL(-s,(mW()*mW()),(mW()*mW())));
vv=vv*2*Rw*(Rw+2.0)-(7.0/6.0+Rw)*lww/s-1.0/6.0-2.0*Rw;
	GSL_SET_COMPLEX(&v,vv.real(),vv.imag());
	return v;
		}


//Eq.68
gsl_complex EWPOSM::piZg()
{gsl_complex a;
    complex zero(0,0,false);
	complex pizg(I3((mZ*mZ),mu*mu,mu*mu));
       
	pizg=(pizg+I3((mZ*mZ),mc*mc,mc*mc)+I3((mZ*mZ),mt*mt,mt*mt))*2.0*3.0*Qf[0]*vf("u");
	pizg=pizg+2.0*3.0*Qf[1]*vf("d")*(zero+ I3((mZ*mZ),md*md,md*md)+I3((mZ*mZ),ms*ms,ms*ms)+I3((mZ*mZ),mb*mb,mb*mb));
	pizg=pizg+2.0*Qf[2]*vf("l")*(zero+I3((mZ*mZ),me*me,me*me)+I3((mZ*mZ),mmu*mmu,mmu*mmu)+I3((mZ*mZ),mtau*mtau,mtau*mtau));
GSL_SET_COMPLEX(&a,pizg.real(),pizg.imag());
return a;}


//Equation 70 of GFitter
double EWPOSM::sigmaBF0()
{
  double rw=mHl*mHl/mW()/mW();
return 5.0*cw2*(1+cw2)/8.0-17.0/4.0+5.0/8.0/cw2-rw/8.0
+(9.0/4.0+3.0/4.0/cw2-3/cw2)*log(cw2)+3.0*rw/4.0/(1.0-rw)*log(rw);
	}
//Eq. 71
gsl_complex EWPOSM::sigmaBFMw()
{gsl_complex j;
complex lwh(LL(-(mW()*mW()),(mW()*mW()),(mHl*mHl)));
complex lwz(LL((mW()*mW()),(mW()*mW()),(mZ*mZ)));
double rw=mHl*mHl/mW()/mW();

lwh=lwh*(0.5-rw/6.0+rw*rw/24.0)/(mW()*mW());
lwz=lwz*(-2*cw2-17.0/6.0+2.0/3.0/cw2+1.0/24.0/cw2/cw2)/(mW()*mW());
//total formula
lwh=lwh	+lwz+(-157.0/9.0+23.0/12.0/cw2+1.0/12.0/cw2/cw2-rw/2.0
+rw*rw/12.0+1.0/cw2*(-7.0/2.0+7.0/12.0/cw2+1.0/24.0/cw2/cw2)*log(cw2)
+rw*(-3.0/4.0+rw/4.0-rw*rw/24.0)*log(rw));
GSL_SET_COMPLEX(&j,lwh.real(),lwh.imag());
return j;
}


//Eq. 72
gsl_complex EWPOSM::sigmaBFMz()
{gsl_complex j;
  double rz=mHl*mHl/mZ/mZ;
  double rw=mHl*mHl/mW()/mW();

complex lzh(LL(-(mZ*mZ),(mZ*mZ),(mHl*mHl)));
complex lww(LL(-(mZ*mZ),(mW()*mW()),(mW()*mW())));
lzh=lzh*(0.5-rz/6.0+rz*rz/24.0)/mW()/mW();
lww=lww*(-2.0*cw2*cw2*cw2-17.0/6.0*cw2*cw2+2.0/3.0*cw2+1.0/24.0)/mW()/mW();
//total formula
lww=lww+lzh-8.0*cw2*cw2-34.0*cw2/3.0+35.0/18.0*(1.0+1.0/cw2)-rw/2.0
+rz*rz/12.0/cw2+rw*(-3.0/4.0+rz/4.0-rz*rz/24.0)*log(rz)
+5.0/6.0/cw2*log(cw2);
GSL_SET_COMPLEX(&j,lww.real(),lww.imag());
return j;
}


//Eq. 59
gsl_complex EWPOSM::rhoZf(const std::string& ferm)
{gsl_complex rr;
if((ferm=="u")||(ferm=="d")||ferm=="l")
	{
	complex r(deltarhorem(ferm));
        complex drmtot(deltarremtot());
	r=(1.0+r)/(1.0+r*(1.0-drmtot));
	GSL_SET_COMPLEX(&rr,r.real(),r.imag());
	}

else {GSL_SET_COMPLEX(&rr,0.0,0.0);std::cout
		<<"wrong input for fermion parameter,  it can be equal only to 'u','d', 'l'   rhoZf"<<"\n";}
	return rr;}

//Eq.60
gsl_complex EWPOSM::kZf(const std::string& ferm)
{
gsl_complex kk;
if((ferm=="u")||(ferm=="d")||ferm=="l")
	{
    complex k(deltakrem(ferm));
    complex drhorem(deltarhorem(ferm));
    complex drrtot(deltarremtot());
    k=(1.0+k)*(1.0-cw2/sw2*drhorem*(1-drrtot));
   GSL_SET_COMPLEX(&kk,k.real(),k.imag());}
   else {GSL_SET_COMPLEX(&kk,0.0,0.0);std::cout
    <<"wrong input for fermion parameter "<< ferm<<"  it can be equal only to 'u','d', 'l'"<<"\n";}
	return kk;}


//Eq. 73
gsl_complex EWPOSM::sigmaBFMz1()
{gsl_complex j;
double rz=mHl*mHl/mZ/mZ;
double rw=mHl*mHl/mW()/mW();
complex lzh(LL(-(mZ*mZ),(mZ*mZ),(mHl*mHl)));
complex lww(LL(-(mZ*mZ),(mW()*mW()),(mW()*mW())));
lzh=lzh*(0.5-5.0*rz/24.0+rz*rz/12.0-1.0/2.0/(rz-4.0))/mW()/mW();
lww=lww*(-cw2*cw2*cw2+7.0/6.0*cw2*cw2-17.0/12.0*cw2-1.0/8.0)/mW()/mW();
//total formula
lww=lww+lzh-4*cw2*cw2+17.0/3.0*cw2-23.0/9.0+
5.0/18.0/cw2-rw/2.0+rw*rz/6.0
+rw*(-3.0/4.0+3.0*rz/8.0-rz*rz/12.0)*log(rz)-1.0/12.0/cw2*log(cw2)
+log(rz)/2.0/cw2;
GSL_SET_COMPLEX(&j,lww.real(),lww.imag());
return j;
}



//Eq. 74
// at the end we will have to write it as a sum ,
// but right now I will write a code only for the first generation;
gsl_complex EWPOSM::sigmaFFMw()
{gsl_complex j;
double muq[3]={mu,mc,mt};
double mdq[3]={md,ms,mb};

//not clear but I guess s=(mW()*mW()) in this case
complex sig(0.0);
for(int i=0;i<3;i++){
complex i3(I3(-(mW()*mW()),muq[i]*muq[i],mdq[i]*mdq[i]));
complex i1(I1(-(mW()*mW()),muq[i]*muq[i],mdq[i]*mdq[i]));
sig=sig+3*(-2.0*i3+i1/(mW()*mW())*(muq[i]*muq[i]+mdq[i]*mdq[i]));	};
GSL_SET_COMPLEX(&j,sig.real(),sig.imag());
	return j;}





//Eq. 75
// at the end we will have to write it as a sum ,
// but right now I will write a code only for the first generation;
gsl_complex EWPOSM::sigmaFFMz()
{gsl_complex j;
//not clear but I guess s=(mZ*mZ) in this case

const double muq[3]={mu,mc,mt};
const double mdq[3]={md,ms,mb};
const double ml[3]={me,mmu,mtau};
complex sig(0.0);
for(int i=0;i<3;i++)
{complex i3u(I3(-(mZ*mZ),muq[i]*muq[i],muq[i]*muq[i]));
 complex i0u(I1(-(mZ*mZ),muq[i]*muq[i],muq[i]*muq[i]));
 complex i3d(I3(-(mZ*mZ),mdq[i]*mdq[i],mdq[i]*mdq[i]));
 complex i0d(I1(-(mZ*mZ),mdq[i]*mdq[i],mdq[i]*mdq[i]));
 complex i3l(I3(-(mZ*mZ),ml[i]*ml[i],ml[i]*ml[i]));
 complex i0l(I1(-(mZ*mZ),ml[i]*ml[i],ml[i]*ml[i]));

sig=sig+1/2.0/cw2*(3.0*(-(1.0+vf("u")*vf("u"))*i3u+i0u*muq[i]*muq[i]/(mZ*mZ))
                   +3.0*(-(1.0+vf("d")*vf("d"))*i3d+i0d*mdq[i]*mdq[i]/(mZ*mZ))
                   +(-(1.0+vf("l")*vf("l"))*i3l+i0l*ml[i]*ml[i]/(mZ*mZ)));

}

GSL_SET_COMPLEX(&j,sig.real(),sig.imag());
	return j;}


//Eq. 76
// at the end we will have to write it as a sum ,
// but right now I will write a code only for the first generation;
gsl_complex EWPOSM::sigmaFFMz1()
{gsl_complex j;
//not clear but I guess s=(mW()*mW()) in this case
complex sig(0.0);
double ru;
double rd;
double rl;
const double muq[3]={mu,mc,mt};
const double mdq[3]={md,ms,mb};
const double ml[3]={me,mmu,mtau};
for(int i=0;i<3;i++){

complex fu(fcurve_c(-(mZ*mZ),muq[i]*muq[i],muq[i]*muq[i]));
complex fd(fcurve_c(-(mZ*mZ),mdq[i]*mdq[i],mdq[i]*mdq[i]));
complex fl(fcurve_c(-(mZ*mZ),ml[i]*ml[i],ml[i]*ml[i]));
ru=muq[i]*muq[i]/(mW()*mW());
rd=mdq[i]*mdq[i]/(mW()*mW());
rl=ml[i]*ml[i]/(mW()*mW());

sig= sig
-3.0*(ru/2.0*(1.0-ru*(mW()*mW())*fu)+1.0/6.0/cw2*(1+vf("u")*vf("u"))*
(0.5*log(ru*cw2)+ru*cw2+(-1.0/4.0/cw2+ru/2.0-ru*ru*cw2)*(mW()*mW())*fu))

-3.0*(rd/2.0*(1.0-rd*(mW()*mW())*fd)+1.0/6.0/cw2*(1+vf("d")*vf("d"))*
(0.5*log(rd*cw2)+rd*cw2+(-1.0/4.0/cw2+rd/2.0-rd*rd*cw2)*(mW()*mW())*fd))

-(rl/2.0*(1.0-rl*(mW()*mW())*fl)+1.0/6.0/cw2*(1+vf("l")*vf("l"))*
(0.5*log(rl*cw2)+rl*cw2+(-1.0/4.0/cw2+rl/2.0-rl*rl*cw2)*(mW()*mW())*fl))
;}
GSL_SET_COMPLEX(&j,sig.real(),sig.imag());
	return j;}

////////////////////////////////////////////////
/**
 *
 * QCD strory from B.A.Kniehl Nucl.Phys. B347 ,86 (1990)
 */

double EWPOSM::F1(double x)
{double b=log(1-x);
const double c1=0.2933;
const double c2=-0.0192;
const double c3=0.4417;
const double c4=-0.2837;
const double f1=-1.1881;
const double f2=-2.0979;
const double f3=4.1157;
const double f4=-2.2082;
const double f5=3.6968;
const double f6=-2.1815;
const double zeta2=1.64493;
	return (1-x)*(1-x)*b*(3.0/8.0*b-zeta2-9.0/8.0)+
	f1+f2*x+f3*x*x+f4*x*x*x+f5*x*x*x*x+f6*x*x*x*x*x;}

double EWPOSM::F1inf(double x)
{ 	const double c1=0.2933;
const double c2=-0.0192;
const double c3=0.4417;
const double c4=-0.2837;
const double f1=-1.1881;
const double f2=-2.0979;
const double f3=4.1157;
const double f4=-2.2082;
const double f5=3.6968;
const double f6=-2.1815;
const double zeta2=1.64493;

    return
	f1+(f2+zeta2+9.0/8.0)*x;}
/**
 * 
 * @return deltaremQCD for equation (61) of gfitter 
 */

double EWPOSM::deltarremQCD()
{       const double c1=0.2933;
const double c2=-0.0192;
const double c3=0.4417;
const double c4=-0.2837;
const double f1=-1.1881;
const double f2=-2.0979;
const double f3=4.1157;
const double f4=-2.2082;
const double f5=3.6968;
const double f6=-2.1815;
const double zeta2=1.64493;
	double ut=mt*mt/(mZ*mZ);
	double ub=mb*mb/(mZ*mZ);
return alsMz/M_PI*(-cw2/sw2*(-1.0/4.0)*(zeta2+0.5)*ale/M_PI*mt*mt/sw2/(mW()*mW())
                 +ale/M_PI*(1.0/2.0/sw2*(-1.0/4.0/sw2+1.0/3.0)*log(ut)
                 -1.0/9.0*log(ub)+c1/sw2/sw2+c2/sw2+c3 +1.0/sw2*(1/sw2-2.0)*
                  (F1(cw2/ut)-F1inf(cw2/ut))));
}
///////////////////////////////////////////////////////////////////////////////////////////////
/**
 *
 * @param lept can we equal to "e", "mu", "tau"
 * @return Gammal partial width of the Z decay to the given lepton
 */
double EWPOSM::Gamma_f(const std::string& lept)
// formula (25) of GFitter
{
double g=0;
double mlept;
	if((lept=="e")||(lept=="mu")||lept=="tau")
        {double ml[3]={me,mmu,mtau};
        double G0=GF*mZ*mZ*mZ/24.0/M_PI/sqrt(2.0);
	if(lept=="e") mlept=ml[0];
	if(lept=="mu") mlept=ml[1];
	if(lept=="tau") mlept=ml[2];
	double rlz=gsl_complex_abs(rhoZf("l"));

	double seff2=sw2;//I will have to refer to the sintheata eff form the code
	seff2=seff2+ale*ale*35.0/18.0*(1-8.0/3.0*GSL_REAL(kZf("l"))*sw2);
	g=G0*rlz*sqrt(1-4*mlept*mlept/(mZ*mZ))*((1+2.0*mlept*mlept/(mZ*mZ))*
	(1.0+/* this ratio of couplings*/4.0*(-0.5+2*seff2)*(-0.5+2*seff2))
	-6.0*mlept*mlept/(mZ*mZ))*(1.0+3.0/4.0*ale*ale/M_PI*1.0);
	  	}
	else {std::cout<<"wrong input, lepton can be only 'e','mu','tau'"<<"\n";}
	return g;}