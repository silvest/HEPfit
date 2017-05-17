#ifndef hpl_h
#define hpl_h

#include <limits>
#include <cmath>
#include <complex>
using namespace std;
typedef complex<double> cd;

const numeric_limits<double> nld = *new numeric_limits<double>;

inline cd operator*(int i, cd x) {return (double)i*x;}
inline cd operator*(cd x, int i) {return (double)i*x;}

inline cd hpl_base(int i, cd x){
if(i==0&&norm(x)<nld.min()) return log(nld.min());
if(i==0) return log(x);
if(i==1&&norm(x)<nld.epsilon()) return x;
if(i==1) return -log(1.-x);
else 
 return 0;
}

inline cd hpl_base(int i1, int i2, cd x){
cd u;if(norm(x)<nld.epsilon()) u = -x;
   else u = log(1.-x);

if(i1==0&&i2==1) return 
-1.*u - 0.25*pow(u,2) - 0.027777777777777776*pow(u,3) + 
   0.0002777777777777778*pow(u,5) - 
   4.72411186696901e-6*pow(u,7) + 
   9.185773074661964e-8*pow(u,9) - 
   1.8978869988971e-9*pow(u,11) + 
   4.0647616451442256e-11*pow(u,13) - 
   8.921691020456452e-13*pow(u,15) + 
   1.9939295860721074e-14*pow(u,17) - 
   4.518980029619918e-16*pow(u,19) + 
   1.0356517612181247e-17*pow(u,21);
else 
 return 0;
}

inline cd hpl_base(int i1, int i2, int i3, cd x){
cd u;if(norm(x)<nld.epsilon()) u = -x;
   else u = log(1.-x);

if(i1==0&&i2==0&&i3==1)
return 
-1.*u - 0.375*pow(u,2) - 0.0787037037037037*pow(u,3) - 
   0.008680555555555556*pow(u,4) - 
   0.00012962962962962963*pow(u,5) + 
   0.00008101851851851852*pow(u,6) + 
   3.4193571608537595e-6*pow(u,7) - 
   1.328656462585034e-6*pow(u,8) - 
   8.660871756109851e-8*pow(u,9) + 
   2.52608759553204e-8*pow(u,10) + 
   2.144694468364065e-9*pow(u,11) - 
   5.140110622012979e-10*pow(u,12) - 
   5.24958211460083e-11*pow(u,13) + 
   1.0887754406636318e-11*pow(u,14) + 
   1.2779396094493695e-12*pow(u,15) - 
   2.369824177308745e-13*pow(u,16) - 
   3.104357887965462e-14*pow(u,17) + 
   5.261758629912506e-15*pow(u,18) + 
   7.538479549949265e-16*pow(u,19) - 
   1.1862322577752286e-16*pow(u,20) - 
   1.8316979965491384e-17*pow(u,21);
if(i1==0&&i2==1&&i3==1)
return 
0.25*pow(u,2) + 0.08333333333333333*pow(u,3) + 
   0.010416666666666666*pow(u,4) - 
   0.00011574074074074075*pow(u,6) + 
   2.066798941798942e-6*pow(u,8) - 
   4.1335978835978836e-8*pow(u,10) + 
   8.698648744945042e-10*pow(u,12) - 
   1.887210763816962e-11*pow(u,14) + 
   4.182042665838962e-13*pow(u,16) - 
   9.415778600896063e-15*pow(u,18) + 
   2.146515514069461e-16*pow(u,20);
else 
 return 0;
}

inline cd hpl_base(int i1, int i2, int i3, int i4, cd x){
cd u;if(norm(x)<nld.epsilon()) u = -x;
   else u = log(1.-x);

if(i1==0&&i2==0&&i3==0&&i4==1)
return 
-1.*u - 0.4375*pow(u,2) - 0.11651234567901235*pow(u,3) - 
   0.019820601851851853*pow(u,4) - 
   0.001927932098765432*pow(u,5) - 
   0.000031057098765432096*pow(u,6) + 
   0.000015624009114857836*pow(u,7) + 
   8.485123546773206e-7*pow(u,8) - 
   2.290961660318971e-7*pow(u,9) - 
   2.1832614218526917e-8*pow(u,10) + 
   3.882824879172015e-9*pow(u,11) + 
   5.446292103220332e-10*pow(u,12) - 
   6.960805210682725e-11*pow(u,13) - 
   1.3375737686445216e-11*pow(u,14) + 
   1.2784852685266572e-12*pow(u,15) + 
   3.260562858024892e-13*pow(u,16) - 
   2.364757116861826e-14*pow(u,17) - 
   7.923135122031162e-15*pow(u,18) + 
   4.3452915709984186e-16*pow(u,19) + 
   1.923627006253592e-16*pow(u,20) - 
   7.812414333195955e-18*pow(u,21);
if(i1==0&&i2==1&&i3==0&&i4==1)
return 
0.25*pow(u,2) + 0.1111111111111111*pow(u,3) + 
   0.022569444444444444*pow(u,4) + 
   0.0020833333333333333*pow(u,5) - 
   0.000027006172839506174*pow(u,6) - 
   0.00001984126984126984*pow(u,7) + 
   4.527273872511968e-7*pow(u,8) + 
   3.389987682315725e-7*pow(u,9) - 
   7.939132443100697e-9*pow(u,10) - 
   6.6805622361177916e-9*pow(u,11) + 
   1.4490216610627064e-10*pow(u,12) + 
   1.39908336457158e-10*pow(u,13) - 
   2.7425719106565973e-12*pow(u,14) - 
   3.032441227329819e-12*pow(u,15) + 
   5.358569182999823e-14*pow(u,16) + 
   6.724068599976371e-14*pow(u,17) - 
   1.0756816626218996e-15*pow(u,18) - 
   1.5158529016922455e-15*pow(u,19) + 
   2.208955024323606e-17*pow(u,20) + 
   3.460964863954937e-17*pow(u,21);
if(i1==0&&i2==1&&i3==1&&i4==1)
return 
-0.05555555555555555*pow(u,3) - 
   0.020833333333333332*pow(u,4) - 
   0.002777777777777778*pow(u,5) + 
   0.00003306878306878307*pow(u,7) - 
   6.123848716441309e-7*pow(u,9) + 
   1.252605419272086e-8*pow(u,11) - 
   2.6765073061369356e-10*pow(u,13) + 
   5.871322376319437e-12*pow(u,15) - 
   1.312013385361243e-13*pow(u,17) + 
   2.97340376870402e-15*pow(u,19) - 
   6.814334965299877e-17*pow(u,21);
else 
 return 0;
}

inline cd hpl(int i, cd x){
if(i==0&&norm(x)<nld.min()) return log(nld.min());
if(i==0) return log(x);
if(i==1&&norm(x)<nld.epsilon()) return x;
if(i==1) return -log(1.-x);
else 
 return 0;
}

inline cd hpl(int i1, int i2, cd x){
const double Pi = M_PI;
if(i1==0&&i2==0){
if(norm(x)>1) return
pow(hpl(0,1./x),2)/2.;
if(real(x)>.5) return
pow(hpl(1,1. - x),2)/2.;
return 
pow(hpl_base(0,x),2)/2.;
}
if(i1==0&&i2==1){
if(norm(x)>1) return
Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + pow(Pi,2)/3. - 
   pow(hpl(0,1./x),2)/2.;
if(real(x)>.5) return
hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + pow(Pi,2)/6.;
return 
hpl_base(0,1,x);
}
if(i1==1&&i2==0){
if(norm(x)>1) return
Pi*cd(0,1)*hpl(0,1./x) - 
   hpl(0,1./x)*(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x)) + 
   hpl(0,1,1./x) - pow(Pi,2)/3. + pow(hpl(0,1./x),2)/2.;
if(real(x)>.5) return
hpl(0,1,1. - x) - pow(Pi,2)/6.;
return 
hpl_base(0,x)*hpl_base(1,x) - hpl_base(0,1,x);
}
if(i1==1&&i2==1){
if(norm(x)>1) return
pow(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x),2)/2.;
if(real(x)>.5) return
pow(hpl(0,1. - x),2)/2.;
return 
pow(hpl_base(1,x),2)/2.;
}
else 
 return 0;
}

inline cd hpl(int i1, int i2, int i3, cd x){
const double Pi = M_PI;
if(i1==0&&i2==0&&i3==0){
if(norm(x)>1) return
-pow(hpl(0,1./x),3)/6.;
if(real(x)>.5) return
-pow(hpl(1,1. - x),3)/6.;
return 
pow(hpl_base(0,x),3)/6.;
}
if(i1==0&&i2==0&&i3==1){
if(norm(x)>1) return
hpl(0,0,1,1./x) - (hpl(0,1./x)*pow(Pi,2))/3. + 
   Pi*cd(0,0.5)*pow(hpl(0,1./x),2) + pow(hpl(0,1./x),3)/6.;
if(real(x)>.5) return
1.2020569031595942 + hpl(1,1. - x)*hpl(0,1,1. - x) - 
   hpl(0,1,1,1. - x) - (hpl(1,1. - x)*pow(Pi,2))/6. - 
   (hpl(0,1. - x)*pow(hpl(1,1. - x),2))/2.;
return 
hpl_base(0,0,1,x);
}
if(i1==0&&i2==1&&i3==0){
if(norm(x)>1) return
-(hpl(0,1./x)*(Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
        pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.)) - 
   2*(hpl(0,0,1,1./x) - (hpl(0,1./x)*pow(Pi,2))/3. + 
      Pi*cd(0,0.5)*pow(hpl(0,1./x),2) + 
      pow(hpl(0,1./x),3)/6.);
if(real(x)>.5) return
-(hpl(1,1. - x)*(hpl(0,1. - x)*hpl(1,1. - x) - 
        hpl(0,1,1. - x) + pow(Pi,2)/6.)) - 
   2*(1.2020569031595942 + hpl(1,1. - x)*hpl(0,1,1. - x) - 
      hpl(0,1,1,1. - x) - (hpl(1,1. - x)*pow(Pi,2))/6. - 
      (hpl(0,1. - x)*pow(hpl(1,1. - x),2))/2.);
return 
hpl_base(0,x)*hpl_base(0,1,x) - 2*hpl_base(0,0,1,x);
}
if(i1==0&&i2==1&&i3==1){
if(norm(x)>1) return
1.2020569031595942 + Pi*cd(0,-1)*hpl(0,1,1./x) - 
   hpl(0,1./x)*hpl(0,1,1./x) + hpl(0,0,1,1./x) - 
   hpl(0,1,1,1./x) + (hpl(0,1./x)*pow(Pi,2))/2. + 
   cd(0,0.16666666666666666)*pow(Pi,3) + 
   Pi*cd(0,-0.5)*pow(hpl(0,1./x),2) - pow(hpl(0,1./x),3)/6.;
if(real(x)>.5) return
1.2020569031595942 + hpl(0,1. - x)*hpl(0,1,1. - x) - 
   hpl(0,0,1,1. - x) - (hpl(1,1. - x)*pow(hpl(0,1. - x),2))/
    2.;
return 
hpl_base(0,1,1,x);
}
if(i1==1&&i2==0&&i3==0){
if(norm(x)>1) return
hpl(0,0,1,1./x) - (hpl(0,1./x)*pow(Pi,2))/3. + 
   hpl(0,1./x)*(Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
      pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.) + 
   Pi*cd(0,0.5)*pow(hpl(0,1./x),2) + 
   ((Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
      pow(hpl(0,1./x),2))/2. + pow(hpl(0,1./x),3)/6.;
if(real(x)>.5) return
1.2020569031595942 + hpl(1,1. - x)*hpl(0,1,1. - x) - 
   hpl(0,1,1,1. - x) + hpl(1,1. - x)*
    (hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
      pow(Pi,2)/6.) - (hpl(1,1. - x)*pow(Pi,2))/6. - 
   hpl(0,1. - x)*pow(hpl(1,1. - x),2);
return 
-(hpl_base(0,x)*hpl_base(0,1,x)) + hpl_base(0,0,1,x) + 
   (hpl_base(1,x)*pow(hpl_base(0,x),2))/2.;
}
if(i1==1&&i2==0&&i3==1){
if(norm(x)>1) return
(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
    (Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
      pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.) - 
   2*(1.2020569031595942 + Pi*cd(0,-1)*hpl(0,1,1./x) - 
      hpl(0,1./x)*hpl(0,1,1./x) + hpl(0,0,1,1./x) - 
      hpl(0,1,1,1./x) + (hpl(0,1./x)*pow(Pi,2))/2. + 
      cd(0,0.16666666666666666)*pow(Pi,3) + 
      Pi*cd(0,-0.5)*pow(hpl(0,1./x),2) - 
      pow(hpl(0,1./x),3)/6.);
if(real(x)>.5) return
-(hpl(0,1. - x)*(hpl(0,1. - x)*hpl(1,1. - x) - 
        hpl(0,1,1. - x) + pow(Pi,2)/6.)) - 
   2*(1.2020569031595942 + hpl(0,1. - x)*hpl(0,1,1. - x) - 
      hpl(0,0,1,1. - x) - 
      (hpl(1,1. - x)*pow(hpl(0,1. - x),2))/2.);
return 
hpl_base(1,x)*hpl_base(0,1,x) - 2*hpl_base(0,1,1,x);
}
if(i1==1&&i2==1&&i3==0){
if(norm(x)>1) return
1.2020569031595942 + Pi*cd(0,-1)*hpl(0,1,1./x) - 
   hpl(0,1./x)*hpl(0,1,1./x) + hpl(0,0,1,1./x) - 
   hpl(0,1,1,1./x) + (hpl(0,1./x)*pow(Pi,2))/2. + 
   cd(0,0.16666666666666666)*pow(Pi,3) - 
   (Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
    (Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
      pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.) + 
   Pi*cd(0,-0.5)*pow(hpl(0,1./x),2) - 
   pow(hpl(0,1./x),3)/6. - 
   (hpl(0,1./x)*pow(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x),
       2))/2.;
if(real(x)>.5) return
1.2020569031595942 + hpl(0,1. - x)*hpl(0,1,1. - x) - 
   hpl(0,0,1,1. - x) + hpl(0,1. - x)*
    (hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
      pow(Pi,2)/6.) - hpl(1,1. - x)*pow(hpl(0,1. - x),2);
return 
-(hpl_base(1,x)*hpl_base(0,1,x)) + hpl_base(0,1,1,x) + 
   (hpl_base(0,x)*pow(hpl_base(1,x),2))/2.;
}
if(i1==1&&i2==1&&i3==1){
if(norm(x)>1) return
pow(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x),3)/6.;
if(real(x)>.5) return
-pow(hpl(0,1. - x),3)/6.;
return 
pow(hpl_base(1,x),3)/6.;
}
else 
 return 0;
}

inline cd hpl(int i1, int i2, int i3, int i4, cd x){
const double Pi = M_PI;
if(i1==0&&i2==0&&i3==0&&i4==0){
if(norm(x)>1){
 return
pow(hpl(0,1./x),4)/24.;
}
if(real(x)>.5){
 return
pow(hpl(1,1. - x),4)/24.;
}
return 
pow(hpl_base(0,x),4)/24.;
}
if(i1==0&&i2==0&&i3==0&&i4==1){
if(norm(x)>1){
 return
-hpl(0,0,0,1,1./x) + pow(Pi,4)/45. + 
   (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
   Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
   pow(hpl(0,1./x),4)/24.;
}
if(real(x)>.5){
 return
-1.2020569031595942*hpl(1,1. - x) + 
   hpl(1,1. - x)*hpl(0,1,1,1. - x) - hpl(0,1,1,1,1. - x) + 
   pow(Pi,4)/90. - (hpl(0,1,1. - x)*pow(hpl(1,1. - x),2))/
    2. + (pow(Pi,2)*pow(hpl(1,1. - x),2))/12. + 
   (hpl(0,1. - x)*pow(hpl(1,1. - x),3))/6.;
}
return 
hpl_base(0,0,0,1,x);
}
if(i1==0&&i2==0&&i3==1&&i4==0){
if(norm(x)>1){
 return
-(hpl(0,1./x)*(hpl(0,0,1,1./x) - 
        (hpl(0,1./x)*pow(Pi,2))/3. + 
        Pi*cd(0,0.5)*pow(hpl(0,1./x),2) + 
        pow(hpl(0,1./x),3)/6.)) - 
   3*(-hpl(0,0,0,1,1./x) + pow(Pi,4)/45. + 
      (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
      Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
      pow(hpl(0,1./x),4)/24.);
}
if(real(x)>.5){
 return
-(hpl(1,1. - x)*(1.2020569031595942 + 
        hpl(1,1. - x)*hpl(0,1,1. - x) - hpl(0,1,1,1. - x) - 
        (hpl(1,1. - x)*pow(Pi,2))/6. - 
        (hpl(0,1. - x)*pow(hpl(1,1. - x),2))/2.)) - 
   3*(-1.2020569031595942*hpl(1,1. - x) + 
      hpl(1,1. - x)*hpl(0,1,1,1. - x) - 
      hpl(0,1,1,1,1. - x) + pow(Pi,4)/90. - 
      (hpl(0,1,1. - x)*pow(hpl(1,1. - x),2))/2. + 
      (pow(Pi,2)*pow(hpl(1,1. - x),2))/12. + 
      (hpl(0,1. - x)*pow(hpl(1,1. - x),3))/6.);
}
return 
hpl_base(0,x)*hpl_base(0,0,1,x) - 3*hpl_base(0,0,0,1,x);
}
if(i1==0&&i2==0&&i3==1&&i4==1){
if(norm(x)>1){
 return
(-2*(3.7763731361630786*cd(0,2) + 
        2.4041138063191885*hpl(0,1./x) + 
        Pi*cd(0,1)*hpl(0,1./x)*hpl(0,1,1./x) + 
        Pi*cd(0,-2)*hpl(0,0,1,1./x) - 
        2*hpl(0,1./x)*hpl(0,0,1,1./x) + 
        4*hpl(0,0,0,1,1./x) + hpl(0,1,0,1,1./x) - 
        (hpl(0,1,1./x)*pow(Pi,2))/3. + pow(Pi,4)/90. + 
        (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. - 
        (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
        Pi*cd(0,0.16666666666666666)*pow(hpl(0,1./x),3) + 
        pow(hpl(0,1./x),4)/24.) + 
     pow(Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
       pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.,2))/4.;
}
if(real(x)>.5){
 return
(-2*(2.4041138063191885*hpl(1,1. - x) + 
        hpl(0,1. - x)*hpl(1,1. - x)*hpl(0,1,1. - x) - 
        hpl(0,1,0,1,1. - x) + 
        (hpl(0,1. - x)*hpl(1,1. - x)*pow(Pi,2))/6. - 
        (hpl(0,1,1. - x)*pow(Pi,2))/6. + pow(Pi,4)/120. + 
        pow(hpl(0,1,1. - x),2) - 
        2*(hpl(1,1. - x)*hpl(0,0,1,1. - x) + 
           (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/
            2. + (-2*hpl(0,1,0,1,1. - x) + 
              pow(hpl(0,1,1. - x),2))/2.) - 
        2*(hpl(0,1. - x)*hpl(0,1,1,1. - x) + 
           (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/
            2. + (-2*hpl(0,1,0,1,1. - x) + 
              pow(hpl(0,1,1. - x),2))/2.)) + 
     pow(hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
       pow(Pi,2)/6.,2))/4.;
}
return 
(-2*hpl_base(0,1,0,1,x) + pow(hpl_base(0,1,x),2))/4.;
}
if(i1==0&&i2==1&&i3==0&&i4==0){
if(norm(x)>1){
 return
((Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + pow(Pi,2)/3. - 
        pow(hpl(0,1./x),2)/2.)*pow(hpl(0,1./x),2) + 
     4*hpl(0,1./x)*(hpl(0,0,1,1./x) - 
        (hpl(0,1./x)*pow(Pi,2))/3. + 
        Pi*cd(0,0.5)*pow(hpl(0,1./x),2) + 
        pow(hpl(0,1./x),3)/6.) + 
     6*(-hpl(0,0,0,1,1./x) + pow(Pi,4)/45. + 
        (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
        Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
        pow(hpl(0,1./x),4)/24.))/2.;
}
if(real(x)>.5){
 return
((hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
        pow(Pi,2)/6.)*pow(hpl(1,1. - x),2) + 
     4*hpl(1,1. - x)*(1.2020569031595942 + 
        hpl(1,1. - x)*hpl(0,1,1. - x) - hpl(0,1,1,1. - x) - 
        (hpl(1,1. - x)*pow(Pi,2))/6. - 
        (hpl(0,1. - x)*pow(hpl(1,1. - x),2))/2.) + 
     6*(-1.2020569031595942*hpl(1,1. - x) + 
        hpl(1,1. - x)*hpl(0,1,1,1. - x) - 
        hpl(0,1,1,1,1. - x) + pow(Pi,4)/90. - 
        (hpl(0,1,1. - x)*pow(hpl(1,1. - x),2))/2. + 
        (pow(Pi,2)*pow(hpl(1,1. - x),2))/12. + 
        (hpl(0,1. - x)*pow(hpl(1,1. - x),3))/6.))/2.;
}
return 
(-4*hpl_base(0,x)*hpl_base(0,0,1,x) + 6*hpl_base(0,0,0,1,x) + 
     hpl_base(0,1,x)*pow(hpl_base(0,x),2))/2.;
}
if(i1==0&&i2==1&&i3==0&&i4==1){
if(norm(x)>1){
 return
3.7763731361630786*cd(0,2) + 
   2.4041138063191885*hpl(0,1./x) + 
   Pi*cd(0,1)*hpl(0,1./x)*hpl(0,1,1./x) + 
   Pi*cd(0,-2)*hpl(0,0,1,1./x) - 
   2*hpl(0,1./x)*hpl(0,0,1,1./x) + 4*hpl(0,0,0,1,1./x) + 
   hpl(0,1,0,1,1./x) - (hpl(0,1,1./x)*pow(Pi,2))/3. + 
   pow(Pi,4)/90. + (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. - 
   (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
   Pi*cd(0,0.16666666666666666)*pow(hpl(0,1./x),3) + 
   pow(hpl(0,1./x),4)/24.;
}
if(real(x)>.5){
 return
2.4041138063191885*hpl(1,1. - x) + 
   hpl(0,1. - x)*hpl(1,1. - x)*hpl(0,1,1. - x) - 
   hpl(0,1,0,1,1. - x) + 
   (hpl(0,1. - x)*hpl(1,1. - x)*pow(Pi,2))/6. - 
   (hpl(0,1,1. - x)*pow(Pi,2))/6. + pow(Pi,4)/120. + 
   pow(hpl(0,1,1. - x),2) - 
   2*(hpl(1,1. - x)*hpl(0,0,1,1. - x) + 
      (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/2. + 
      (-2*hpl(0,1,0,1,1. - x) + pow(hpl(0,1,1. - x),2))/2.)\
    - 2*(hpl(0,1. - x)*hpl(0,1,1,1. - x) + 
      (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/2. + 
      (-2*hpl(0,1,0,1,1. - x) + pow(hpl(0,1,1. - x),2))/2.);
}
return 
hpl_base(0,1,0,1,x);
}
if(i1==0&&i2==1&&i3==1&&i4==0){
if(norm(x)>1){
 return
3.7763731361630786*cd(0,-2) - 
   2.4041138063191885*hpl(0,1./x) + 
   Pi*cd(0,-1)*hpl(0,1./x)*hpl(0,1,1./x) + 
   Pi*cd(0,2)*hpl(0,0,1,1./x) + 
   2*hpl(0,1./x)*hpl(0,0,1,1./x) - 4*hpl(0,0,0,1,1./x) - 
   hpl(0,1,0,1,1./x) + (hpl(0,1,1./x)*pow(Pi,2))/3. - 
   pow(Pi,4)/90. - (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. + 
   (pow(Pi,2)*pow(hpl(0,1./x),2))/6. - 
   hpl(0,1./x)*(1.2020569031595942 + 
      Pi*cd(0,-1)*hpl(0,1,1./x) - 
      hpl(0,1./x)*hpl(0,1,1./x) + hpl(0,0,1,1./x) - 
      hpl(0,1,1,1./x) + (hpl(0,1./x)*pow(Pi,2))/2. + 
      cd(0,0.16666666666666666)*pow(Pi,3) + 
      Pi*cd(0,-0.5)*pow(hpl(0,1./x),2) - 
      pow(hpl(0,1./x),3)/6.) + 
   Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
   pow(hpl(0,1./x),4)/24. + 
   (2*(3.7763731361630786*cd(0,2) + 
         2.4041138063191885*hpl(0,1./x) + 
         Pi*cd(0,1)*hpl(0,1./x)*hpl(0,1,1./x) + 
         Pi*cd(0,-2)*hpl(0,0,1,1./x) - 
         2*hpl(0,1./x)*hpl(0,0,1,1./x) + 
         4*hpl(0,0,0,1,1./x) + hpl(0,1,0,1,1./x) - 
         (hpl(0,1,1./x)*pow(Pi,2))/3. + pow(Pi,4)/90. + 
         (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. - 
         (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
         Pi*cd(0,0.16666666666666666)*pow(hpl(0,1./x),3) + 
         pow(hpl(0,1./x),4)/24.) - 
      pow(Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
        pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.,2))/2.;
}
if(real(x)>.5){
 return
-2.4041138063191885*hpl(1,1. - x) - 
   hpl(0,1. - x)*hpl(1,1. - x)*hpl(0,1,1. - x) + 
   hpl(0,1,0,1,1. - x) - 
   (hpl(0,1. - x)*hpl(1,1. - x)*pow(Pi,2))/6. + 
   (hpl(0,1,1. - x)*pow(Pi,2))/6. - pow(Pi,4)/120. - 
   hpl(1,1. - x)*(1.2020569031595942 + 
      hpl(0,1. - x)*hpl(0,1,1. - x) - hpl(0,0,1,1. - x) - 
      (hpl(1,1. - x)*pow(hpl(0,1. - x),2))/2.) - 
   pow(hpl(0,1,1. - x),2) + 
   2*(hpl(1,1. - x)*hpl(0,0,1,1. - x) + 
      (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/2. + 
      (-2*hpl(0,1,0,1,1. - x) + pow(hpl(0,1,1. - x),2))/2.)\
    + 2*(hpl(0,1. - x)*hpl(0,1,1,1. - x) + 
      (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/2. + 
      (-2*hpl(0,1,0,1,1. - x) + pow(hpl(0,1,1. - x),2))/2.)\
    + (2*(2.4041138063191885*hpl(1,1. - x) + 
         hpl(0,1. - x)*hpl(1,1. - x)*hpl(0,1,1. - x) - 
         hpl(0,1,0,1,1. - x) + 
         (hpl(0,1. - x)*hpl(1,1. - x)*pow(Pi,2))/6. - 
         (hpl(0,1,1. - x)*pow(Pi,2))/6. + pow(Pi,4)/120. + 
         pow(hpl(0,1,1. - x),2) - 
         2*(hpl(1,1. - x)*hpl(0,0,1,1. - x) + 
            (2*hpl(0,1,0,1,1. - x) - 
               pow(hpl(0,1,1. - x),2))/2. + 
            (-2*hpl(0,1,0,1,1. - x) + 
               pow(hpl(0,1,1. - x),2))/2.) - 
         2*(hpl(0,1. - x)*hpl(0,1,1,1. - x) + 
            (2*hpl(0,1,0,1,1. - x) - 
               pow(hpl(0,1,1. - x),2))/2. + 
            (-2*hpl(0,1,0,1,1. - x) + 
               pow(hpl(0,1,1. - x),2))/2.)) - 
      pow(hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
        pow(Pi,2)/6.,2))/2.;
}
return 
hpl_base(0,x)*hpl_base(0,1,1,x) - hpl_base(0,1,0,1,x) + 
   (2*hpl_base(0,1,0,1,x) - pow(hpl_base(0,1,x),2))/2.;
}
if(i1==0&&i2==1&&i3==1&&i4==1){
if(norm(x)>1){
 return
Pi*cd(0,-1)*hpl(0,1./x)*hpl(0,1,1./x) + 
   Pi*cd(0,1)*hpl(0,0,1,1./x) + 
   hpl(0,1./x)*hpl(0,0,1,1./x) + 
   Pi*cd(0,-1)*hpl(0,1,1,1./x) - 
   hpl(0,1./x)*hpl(0,1,1,1./x) - hpl(0,0,0,1,1./x) - 
   hpl(0,1,0,1,1./x)/2. - hpl(0,1,1,1,1./x) + 
   (hpl(0,1,1./x)*pow(Pi,2))/2. + 
   cd(0,0.16666666666666666)*hpl(0,1./x)*pow(Pi,3) - 
   (19*pow(Pi,4))/360. - 
   (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. + 
   (pow(Pi,2)*pow(hpl(0,1./x),2))/4. + 
   Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
   pow(hpl(0,1./x),4)/24. + 
   (2*hpl(0,1,0,1,1./x) - pow(hpl(0,1,1./x),2))/2. + 
   pow(hpl(0,1,1./x),2)/4. + 
   (-2*hpl(0,1,0,1,1./x) + pow(hpl(0,1,1./x),2))/2.;
}
if(real(x)>.5){
 return
hpl(0,1. - x)*hpl(0,0,1,1. - x) - hpl(0,0,0,1,1. - x) + 
   pow(Pi,4)/90. - (hpl(0,1,1. - x)*pow(hpl(0,1. - x),2))/
    2. + (hpl(1,1. - x)*pow(hpl(0,1. - x),3))/6.;
}
return 
hpl_base(0,1,1,1,x);
}
if(i1==1&&i2==0&&i3==0&&i4==0){
if(norm(x)>1){
 return
(-3*(Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
        pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.)*
      pow(hpl(0,1./x),2) - 
     6*hpl(0,1./x)*(hpl(0,0,1,1./x) - 
        (hpl(0,1./x)*pow(Pi,2))/3. + 
        Pi*cd(0,0.5)*pow(hpl(0,1./x),2) + 
        pow(hpl(0,1./x),3)/6.) - 
     (Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
      pow(hpl(0,1./x),3) - 
     6*(-hpl(0,0,0,1,1./x) + pow(Pi,4)/45. + 
        (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
        Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
        pow(hpl(0,1./x),4)/24.))/6.;
}
if(real(x)>.5){
 return
(-3*(hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
        pow(Pi,2)/6.)*pow(hpl(1,1. - x),2) - 
     6*hpl(1,1. - x)*(1.2020569031595942 + 
        hpl(1,1. - x)*hpl(0,1,1. - x) - hpl(0,1,1,1. - x) - 
        (hpl(1,1. - x)*pow(Pi,2))/6. - 
        (hpl(0,1. - x)*pow(hpl(1,1. - x),2))/2.) + 
     hpl(0,1. - x)*pow(hpl(1,1. - x),3) - 
     6*(-1.2020569031595942*hpl(1,1. - x) + 
        hpl(1,1. - x)*hpl(0,1,1,1. - x) - 
        hpl(0,1,1,1,1. - x) + pow(Pi,4)/90. - 
        (hpl(0,1,1. - x)*pow(hpl(1,1. - x),2))/2. + 
        (pow(Pi,2)*pow(hpl(1,1. - x),2))/12. + 
        (hpl(0,1. - x)*pow(hpl(1,1. - x),3))/6.))/6.;
}
return 
(6*hpl_base(0,x)*hpl_base(0,0,1,x) - 6*hpl_base(0,0,0,1,x) - 
     3*hpl_base(0,1,x)*pow(hpl_base(0,x),2) + hpl_base(1,x)*pow(hpl_base(0,x),3)
     )/6.;
}
if(i1==1&&i2==0&&i3==0&&i4==1){
if(norm(x)>1){
 return
3.7763731361630786*cd(0,-2) - 
   2.4041138063191885*hpl(0,1./x) + 
   Pi*cd(0,-1)*hpl(0,1./x)*hpl(0,1,1./x) + 
   Pi*cd(0,2)*hpl(0,0,1,1./x) + 
   2*hpl(0,1./x)*hpl(0,0,1,1./x) - 4*hpl(0,0,0,1,1./x) - 
   hpl(0,1,0,1,1./x) + (hpl(0,1,1./x)*pow(Pi,2))/3. - 
   pow(Pi,4)/90. - (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. + 
   (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
   (Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
    (hpl(0,0,1,1./x) - (hpl(0,1./x)*pow(Pi,2))/3. + 
      Pi*cd(0,0.5)*pow(hpl(0,1./x),2) + 
      pow(hpl(0,1./x),3)/6.) + 
   Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
   pow(hpl(0,1./x),4)/24. + 
   (2*(3.7763731361630786*cd(0,2) + 
         2.4041138063191885*hpl(0,1./x) + 
         Pi*cd(0,1)*hpl(0,1./x)*hpl(0,1,1./x) + 
         Pi*cd(0,-2)*hpl(0,0,1,1./x) - 
         2*hpl(0,1./x)*hpl(0,0,1,1./x) + 
         4*hpl(0,0,0,1,1./x) + hpl(0,1,0,1,1./x) - 
         (hpl(0,1,1./x)*pow(Pi,2))/3. + pow(Pi,4)/90. + 
         (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. - 
         (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
         Pi*cd(0,0.16666666666666666)*pow(hpl(0,1./x),3) + 
         pow(hpl(0,1./x),4)/24.) - 
      pow(Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
        pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.,2))/2.;
}
if(real(x)>.5){
 return
-2.4041138063191885*hpl(1,1. - x) - 
   hpl(0,1. - x)*hpl(1,1. - x)*hpl(0,1,1. - x) + 
   hpl(0,1,0,1,1. - x) - 
   (hpl(0,1. - x)*hpl(1,1. - x)*pow(Pi,2))/6. + 
   (hpl(0,1,1. - x)*pow(Pi,2))/6. - pow(Pi,4)/120. - 
   hpl(0,1. - x)*(1.2020569031595942 + 
      hpl(1,1. - x)*hpl(0,1,1. - x) - hpl(0,1,1,1. - x) - 
      (hpl(1,1. - x)*pow(Pi,2))/6. - 
      (hpl(0,1. - x)*pow(hpl(1,1. - x),2))/2.) - 
   pow(hpl(0,1,1. - x),2) + 
   2*(hpl(1,1. - x)*hpl(0,0,1,1. - x) + 
      (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/2. + 
      (-2*hpl(0,1,0,1,1. - x) + pow(hpl(0,1,1. - x),2))/2.)\
    + 2*(hpl(0,1. - x)*hpl(0,1,1,1. - x) + 
      (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/2. + 
      (-2*hpl(0,1,0,1,1. - x) + pow(hpl(0,1,1. - x),2))/2.)\
    + (2*(2.4041138063191885*hpl(1,1. - x) + 
         hpl(0,1. - x)*hpl(1,1. - x)*hpl(0,1,1. - x) - 
         hpl(0,1,0,1,1. - x) + 
         (hpl(0,1. - x)*hpl(1,1. - x)*pow(Pi,2))/6. - 
         (hpl(0,1,1. - x)*pow(Pi,2))/6. + pow(Pi,4)/120. + 
         pow(hpl(0,1,1. - x),2) - 
         2*(hpl(1,1. - x)*hpl(0,0,1,1. - x) + 
            (2*hpl(0,1,0,1,1. - x) - 
               pow(hpl(0,1,1. - x),2))/2. + 
            (-2*hpl(0,1,0,1,1. - x) + 
               pow(hpl(0,1,1. - x),2))/2.) - 
         2*(hpl(0,1. - x)*hpl(0,1,1,1. - x) + 
            (2*hpl(0,1,0,1,1. - x) - 
               pow(hpl(0,1,1. - x),2))/2. + 
            (-2*hpl(0,1,0,1,1. - x) + 
               pow(hpl(0,1,1. - x),2))/2.)) - 
      pow(hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
        pow(Pi,2)/6.,2))/2.;
}
return 
hpl_base(1,x)*hpl_base(0,0,1,x) - hpl_base(0,1,0,1,x) + 
   (2*hpl_base(0,1,0,1,x) - pow(hpl_base(0,1,x),2))/2.;
}
if(i1==1&&i2==0&&i3==1&&i4==0){
if(norm(x)>1){
 return
3.7763731361630786*cd(0,-2) - 
   2.4041138063191885*hpl(0,1./x) + 
   Pi*cd(0,-1)*hpl(0,1./x)*hpl(0,1,1./x) + 
   Pi*cd(0,2)*hpl(0,0,1,1./x) + 
   2*hpl(0,1./x)*hpl(0,0,1,1./x) - 4*hpl(0,0,0,1,1./x) - 
   hpl(0,1,0,1,1./x) + (hpl(0,1,1./x)*pow(Pi,2))/3. - 
   pow(Pi,4)/90. - hpl(0,1./x)*
    (Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
    (Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
      pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.) - 
   (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. + 
   (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
   2*hpl(0,1./x)*(1.2020569031595942 + 
      Pi*cd(0,-1)*hpl(0,1,1./x) - 
      hpl(0,1./x)*hpl(0,1,1./x) + hpl(0,0,1,1./x) - 
      hpl(0,1,1,1./x) + (hpl(0,1./x)*pow(Pi,2))/2. + 
      cd(0,0.16666666666666666)*pow(Pi,3) + 
      Pi*cd(0,-0.5)*pow(hpl(0,1./x),2) - 
      pow(hpl(0,1./x),3)/6.) - 
   2*(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
    (hpl(0,0,1,1./x) - (hpl(0,1./x)*pow(Pi,2))/3. + 
      Pi*cd(0,0.5)*pow(hpl(0,1./x),2) + 
      pow(hpl(0,1./x),3)/6.) + 
   Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
   pow(hpl(0,1./x),4)/24. + 
   pow(Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
     pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.,2);
}
if(real(x)>.5){
 return
-2.4041138063191885*hpl(1,1. - x) - 
   hpl(0,1. - x)*hpl(1,1. - x)*hpl(0,1,1. - x) + 
   hpl(0,1,0,1,1. - x) + 
   hpl(0,1. - x)*hpl(1,1. - x)*
    (hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
      pow(Pi,2)/6.) - (hpl(0,1. - x)*hpl(1,1. - x)*
      pow(Pi,2))/6. + (hpl(0,1,1. - x)*pow(Pi,2))/6. - 
   pow(Pi,4)/120. + 2*hpl(1,1. - x)*
    (1.2020569031595942 + hpl(0,1. - x)*hpl(0,1,1. - x) - 
      hpl(0,0,1,1. - x) - 
      (hpl(1,1. - x)*pow(hpl(0,1. - x),2))/2.) + 
   2*hpl(0,1. - x)*(1.2020569031595942 + 
      hpl(1,1. - x)*hpl(0,1,1. - x) - hpl(0,1,1,1. - x) - 
      (hpl(1,1. - x)*pow(Pi,2))/6. - 
      (hpl(0,1. - x)*pow(hpl(1,1. - x),2))/2.) - 
   pow(hpl(0,1,1. - x),2) + 
   2*(hpl(1,1. - x)*hpl(0,0,1,1. - x) + 
      (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/2. + 
      (-2*hpl(0,1,0,1,1. - x) + pow(hpl(0,1,1. - x),2))/2.)\
    + 2*(hpl(0,1. - x)*hpl(0,1,1,1. - x) + 
      (2*hpl(0,1,0,1,1. - x) - pow(hpl(0,1,1. - x),2))/2. + 
      (-2*hpl(0,1,0,1,1. - x) + pow(hpl(0,1,1. - x),2))/2.)\
    + pow(hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
     pow(Pi,2)/6.,2);
}
return 
hpl_base(0,x)*hpl_base(1,x)*hpl_base(0,1,x) - 2*hpl_base(1,x)*hpl_base(0,0,1,x) - 
   2*hpl_base(0,x)*hpl_base(0,1,1,x) - hpl_base(0,1,0,1,x) + 
   pow(hpl_base(0,1,x),2);
}
if(i1==1&&i2==0&&i3==1&&i4==1){
if(norm(x)>1){
 return
(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
    (1.2020569031595942 + Pi*cd(0,-1)*hpl(0,1,1./x) - 
      hpl(0,1./x)*hpl(0,1,1./x) + hpl(0,0,1,1./x) - 
      hpl(0,1,1,1./x) + (hpl(0,1./x)*pow(Pi,2))/2. + 
      cd(0,0.16666666666666666)*pow(Pi,3) + 
      Pi*cd(0,-0.5)*pow(hpl(0,1./x),2) - 
      pow(hpl(0,1./x),3)/6.) - 
   3*(Pi*cd(0,-1)*hpl(0,1./x)*hpl(0,1,1./x) + 
      Pi*cd(0,1)*hpl(0,0,1,1./x) + 
      hpl(0,1./x)*hpl(0,0,1,1./x) + 
      Pi*cd(0,-1)*hpl(0,1,1,1./x) - 
      hpl(0,1./x)*hpl(0,1,1,1./x) - hpl(0,0,0,1,1./x) - 
      hpl(0,1,0,1,1./x)/2. - hpl(0,1,1,1,1./x) + 
      (hpl(0,1,1./x)*pow(Pi,2))/2. + 
      cd(0,0.16666666666666666)*hpl(0,1./x)*pow(Pi,3) - 
      (19*pow(Pi,4))/360. - 
      (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. + 
      (pow(Pi,2)*pow(hpl(0,1./x),2))/4. + 
      Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
      pow(hpl(0,1./x),4)/24. + 
      (2*hpl(0,1,0,1,1./x) - pow(hpl(0,1,1./x),2))/2. + 
      pow(hpl(0,1,1./x),2)/4. + 
      (-2*hpl(0,1,0,1,1./x) + pow(hpl(0,1,1./x),2))/2.);
}
if(real(x)>.5){
 return
-(hpl(0,1. - x)*(1.2020569031595942 + 
        hpl(0,1. - x)*hpl(0,1,1. - x) - hpl(0,0,1,1. - x) - 
        (hpl(1,1. - x)*pow(hpl(0,1. - x),2))/2.)) - 
   3*(hpl(0,1. - x)*hpl(0,0,1,1. - x) - 
      hpl(0,0,0,1,1. - x) + pow(Pi,4)/90. - 
      (hpl(0,1,1. - x)*pow(hpl(0,1. - x),2))/2. + 
      (hpl(1,1. - x)*pow(hpl(0,1. - x),3))/6.);
}
return 
hpl_base(1,x)*hpl_base(0,1,1,x) - 3*hpl_base(0,1,1,1,x);
}
if(i1==1&&i2==1&&i3==0&&i4==0){
if(norm(x)>1){
 return
hpl(0,1./x)*(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
    (Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
      pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.) - 
   hpl(0,1./x)*(1.2020569031595942 + 
      Pi*cd(0,-1)*hpl(0,1,1./x) - 
      hpl(0,1./x)*hpl(0,1,1./x) + hpl(0,0,1,1./x) - 
      hpl(0,1,1,1./x) + (hpl(0,1./x)*pow(Pi,2))/2. + 
      cd(0,0.16666666666666666)*pow(Pi,3) + 
      Pi*cd(0,-0.5)*pow(hpl(0,1./x),2) - 
      pow(hpl(0,1./x),3)/6.) + 
   (Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
    (hpl(0,0,1,1./x) - (hpl(0,1./x)*pow(Pi,2))/3. + 
      Pi*cd(0,0.5)*pow(hpl(0,1./x),2) + 
      pow(hpl(0,1./x),3)/6.) + 
   (pow(hpl(0,1./x),2)*pow(Pi*cd(0,1) + hpl(0,1./x) + 
        hpl(1,1./x),2))/4. + 
   (2*(3.7763731361630786*cd(0,2) + 
         2.4041138063191885*hpl(0,1./x) + 
         Pi*cd(0,1)*hpl(0,1./x)*hpl(0,1,1./x) + 
         Pi*cd(0,-2)*hpl(0,0,1,1./x) - 
         2*hpl(0,1./x)*hpl(0,0,1,1./x) + 
         4*hpl(0,0,0,1,1./x) + hpl(0,1,0,1,1./x) - 
         (hpl(0,1,1./x)*pow(Pi,2))/3. + pow(Pi,4)/90. + 
         (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. - 
         (pow(Pi,2)*pow(hpl(0,1./x),2))/6. + 
         Pi*cd(0,0.16666666666666666)*pow(hpl(0,1./x),3) + 
         pow(hpl(0,1./x),4)/24.) - 
      pow(Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
        pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.,2))/4.;
}
if(real(x)>.5){
 return
-(hpl(0,1. - x)*hpl(1,1. - x)*
      (hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
        pow(Pi,2)/6.)) - 
   hpl(1,1. - x)*(1.2020569031595942 + 
      hpl(0,1. - x)*hpl(0,1,1. - x) - hpl(0,0,1,1. - x) - 
      (hpl(1,1. - x)*pow(hpl(0,1. - x),2))/2.) + 
   (pow(hpl(0,1. - x),2)*pow(hpl(1,1. - x),2))/4. - 
   hpl(0,1. - x)*(1.2020569031595942 + 
      hpl(1,1. - x)*hpl(0,1,1. - x) - hpl(0,1,1,1. - x) - 
      (hpl(1,1. - x)*pow(Pi,2))/6. - 
      (hpl(0,1. - x)*pow(hpl(1,1. - x),2))/2.) + 
   (2*(2.4041138063191885*hpl(1,1. - x) + 
         hpl(0,1. - x)*hpl(1,1. - x)*hpl(0,1,1. - x) - 
         hpl(0,1,0,1,1. - x) + 
         (hpl(0,1. - x)*hpl(1,1. - x)*pow(Pi,2))/6. - 
         (hpl(0,1,1. - x)*pow(Pi,2))/6. + pow(Pi,4)/120. + 
         pow(hpl(0,1,1. - x),2) - 
         2*(hpl(1,1. - x)*hpl(0,0,1,1. - x) + 
            (2*hpl(0,1,0,1,1. - x) - 
               pow(hpl(0,1,1. - x),2))/2. + 
            (-2*hpl(0,1,0,1,1. - x) + 
               pow(hpl(0,1,1. - x),2))/2.) - 
         2*(hpl(0,1. - x)*hpl(0,1,1,1. - x) + 
            (2*hpl(0,1,0,1,1. - x) - 
               pow(hpl(0,1,1. - x),2))/2. + 
            (-2*hpl(0,1,0,1,1. - x) + 
               pow(hpl(0,1,1. - x),2))/2.)) - 
      pow(hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
        pow(Pi,2)/6.,2))/4.;
}
return 
-(hpl_base(0,x)*hpl_base(1,x)*hpl_base(0,1,x)) + hpl_base(1,x)*hpl_base(0,0,1,x) + 
   hpl_base(0,x)*hpl_base(0,1,1,x) + 
   (pow(hpl_base(0,x),2)*pow(hpl_base(1,x),2))/4. + 
   (2*hpl_base(0,1,0,1,x) - pow(hpl_base(0,1,x),2))/4.;
}
if(i1==1&&i2==1&&i3==0&&i4==1){
if(norm(x)>1){
 return
(-4*(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
      (1.2020569031595942 + Pi*cd(0,-1)*hpl(0,1,1./x) - 
        hpl(0,1./x)*hpl(0,1,1./x) + hpl(0,0,1,1./x) - 
        hpl(0,1,1,1./x) + (hpl(0,1./x)*pow(Pi,2))/2. + 
        cd(0,0.16666666666666666)*pow(Pi,3) + 
        Pi*cd(0,-0.5)*pow(hpl(0,1./x),2) - 
        pow(hpl(0,1./x),3)/6.) + 
     (Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
        pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.)*
      pow(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x),2) + 
     6*(Pi*cd(0,-1)*hpl(0,1./x)*hpl(0,1,1./x) + 
        Pi*cd(0,1)*hpl(0,0,1,1./x) + 
        hpl(0,1./x)*hpl(0,0,1,1./x) + 
        Pi*cd(0,-1)*hpl(0,1,1,1./x) - 
        hpl(0,1./x)*hpl(0,1,1,1./x) - hpl(0,0,0,1,1./x) - 
        hpl(0,1,0,1,1./x)/2. - hpl(0,1,1,1,1./x) + 
        (hpl(0,1,1./x)*pow(Pi,2))/2. + 
        cd(0,0.16666666666666666)*hpl(0,1./x)*pow(Pi,3) - 
        (19*pow(Pi,4))/360. - 
        (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. + 
        (pow(Pi,2)*pow(hpl(0,1./x),2))/4. + 
        Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
        pow(hpl(0,1./x),4)/24. + 
        (2*hpl(0,1,0,1,1./x) - pow(hpl(0,1,1./x),2))/2. + 
        pow(hpl(0,1,1./x),2)/4. + 
        (-2*hpl(0,1,0,1,1./x) + pow(hpl(0,1,1./x),2))/2.))/
   2.;
}
if(real(x)>.5){
 return
((hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
        pow(Pi,2)/6.)*pow(hpl(0,1. - x),2) + 
     4*hpl(0,1. - x)*(1.2020569031595942 + 
        hpl(0,1. - x)*hpl(0,1,1. - x) - hpl(0,0,1,1. - x) - 
        (hpl(1,1. - x)*pow(hpl(0,1. - x),2))/2.) + 
     6*(hpl(0,1. - x)*hpl(0,0,1,1. - x) - 
        hpl(0,0,0,1,1. - x) + pow(Pi,4)/90. - 
        (hpl(0,1,1. - x)*pow(hpl(0,1. - x),2))/2. + 
        (hpl(1,1. - x)*pow(hpl(0,1. - x),3))/6.))/2.;
}
return 
(-4*hpl_base(1,x)*hpl_base(0,1,1,x) + 6*hpl_base(0,1,1,1,x) + 
     hpl_base(0,1,x)*pow(hpl_base(1,x),2))/2.;
}
if(i1==1&&i2==1&&i3==1&&i4==0){
if(norm(x)>1){
 return
(6*(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x))*
      (1.2020569031595942 + Pi*cd(0,-1)*hpl(0,1,1./x) - 
        hpl(0,1./x)*hpl(0,1,1./x) + hpl(0,0,1,1./x) - 
        hpl(0,1,1,1./x) + (hpl(0,1./x)*pow(Pi,2))/2. + 
        cd(0,0.16666666666666666)*pow(Pi,3) + 
        Pi*cd(0,-0.5)*pow(hpl(0,1./x),2) - 
        pow(hpl(0,1./x),3)/6.) - 
     3*(Pi*cd(0,-1)*hpl(0,1./x) - hpl(0,1,1./x) + 
        pow(Pi,2)/3. - pow(hpl(0,1./x),2)/2.)*
      pow(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x),2) - 
     hpl(0,1./x)*pow(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x),
       3) - 6*(Pi*cd(0,-1)*hpl(0,1./x)*hpl(0,1,1./x) + 
        Pi*cd(0,1)*hpl(0,0,1,1./x) + 
        hpl(0,1./x)*hpl(0,0,1,1./x) + 
        Pi*cd(0,-1)*hpl(0,1,1,1./x) - 
        hpl(0,1./x)*hpl(0,1,1,1./x) - hpl(0,0,0,1,1./x) - 
        hpl(0,1,0,1,1./x)/2. - hpl(0,1,1,1,1./x) + 
        (hpl(0,1,1./x)*pow(Pi,2))/2. + 
        cd(0,0.16666666666666666)*hpl(0,1./x)*pow(Pi,3) - 
        (19*pow(Pi,4))/360. - 
        (hpl(0,1,1./x)*pow(hpl(0,1./x),2))/2. + 
        (pow(Pi,2)*pow(hpl(0,1./x),2))/4. + 
        Pi*cd(0,-0.16666666666666666)*pow(hpl(0,1./x),3) - 
        pow(hpl(0,1./x),4)/24. + 
        (2*hpl(0,1,0,1,1./x) - pow(hpl(0,1,1./x),2))/2. + 
        pow(hpl(0,1,1./x),2)/4. + 
        (-2*hpl(0,1,0,1,1./x) + pow(hpl(0,1,1./x),2))/2.))/
   6.;
}
if(real(x)>.5){
 return
(-3*(hpl(0,1. - x)*hpl(1,1. - x) - hpl(0,1,1. - x) + 
        pow(Pi,2)/6.)*pow(hpl(0,1. - x),2) - 
     6*hpl(0,1. - x)*(1.2020569031595942 + 
        hpl(0,1. - x)*hpl(0,1,1. - x) - hpl(0,0,1,1. - x) - 
        (hpl(1,1. - x)*pow(hpl(0,1. - x),2))/2.) + 
     hpl(1,1. - x)*pow(hpl(0,1. - x),3) - 
     6*(hpl(0,1. - x)*hpl(0,0,1,1. - x) - 
        hpl(0,0,0,1,1. - x) + pow(Pi,4)/90. - 
        (hpl(0,1,1. - x)*pow(hpl(0,1. - x),2))/2. + 
        (hpl(1,1. - x)*pow(hpl(0,1. - x),3))/6.))/6.;
}
return 
(6*hpl_base(1,x)*hpl_base(0,1,1,x) - 6*hpl_base(0,1,1,1,x) - 
     3*hpl_base(0,1,x)*pow(hpl_base(1,x),2) + hpl_base(0,x)*pow(hpl_base(1,x),3)
     )/6.;
}
if(i1==1&&i2==1&&i3==1&&i4==1){
if(norm(x)>1){
 return
pow(Pi*cd(0,1) + hpl(0,1./x) + hpl(1,1./x),4)/24.;
}
if(real(x)>.5){
 return
pow(hpl(0,1. - x),4)/24.;
}
return 
pow(hpl_base(1,x),4)/24.;
}
else 
 return 0;
}

inline cd Li2(cd x) {return hpl(0,1,x);}
inline cd Li3(cd x) {return hpl(0,0,1,x);}
inline cd Li4(cd x) {return hpl(0,0,0,1,x);}
inline double Cl2(double x) {return imag(Li2(exp(cd(0,x))));}
inline double Cl3(double x) {return real(Li3(exp(cd(0,x))));}

#endif
