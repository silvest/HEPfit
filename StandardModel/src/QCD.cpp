/* 
 * File:   QCD.cpp
 * Author: marco
 * 
 * Created on February 17, 2011, 2:13 PM
 */

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>
#include "QCD.h"

//void QCD::init(double AlsM_i, double M_i, double mu_i, double md_i, double ms_i,
//        double mc_i, double mb_i, double mt_i, double mu1_i, double mu2_i,
//        double mu3_i) {
//    Nc = 3.;
//    AlsM = AlsM_i;
//    M = M_i;
//    setMass(UP, mu_i);
//    setMass(DOWN, md_i);
//    setMass(CHARM, mc_i);
//    setMass(STRANGE, ms_i);
//    setMass(TOP, mt_i);
//    setMass(BOTTOM, mb_i);
//    mu1 = mu1_i;
//    mu2 = mu2_i;
//    mu3 = mu3_i;
//}

const std::string QCD::QCDvars[NQCDvars] = {"AlsMz","Mz","mup","mdown","mcharm",
"mstrange","mtop","mbottom","mu1_qcd","mu2_qcd","mu3_qcd","MBd","MBs","MBp",
"MK0","MKp","FBd","BBd"};

void QCD::update(const Parameters& Par) {
    double m_i;
    std::vector <double> v;
    if(Par.Find("AlsMz")!=-1) Par.Get("AlsMz",AlsM);
    if(Par.Find("Mz")!=-1) Par.Get("Mz",M);
    if(Par.Find("mup")!=-1) {
        Par.Get("mup",m_i);
        setMass(UP, m_i);
    }
    if(Par.Find("mdown")!=-1) {
        Par.Get("mdown",m_i);
        setMass(DOWN, m_i);
    }
    if(Par.Find("mcharm")!=-1) {
        Par.Get("mcharm",m_i);
        setMass(CHARM, m_i);
    }
    if(Par.Find("mstrange")!=-1) {
        Par.Get("mstrange",m_i);
        setMass(STRANGE, m_i);
    }
    if(Par.Find("mtop")!=-1) {
        Par.Get("mtop",m_i);
        setMass(TOP, m_i);
    }
    if(Par.Find("mbottom")!=-1) {
        Par.Get("mbottom",m_i);
        setMass(BOTTOM, m_i);
    }
    if(Par.Find("mu1_qcd")!=-1) Par.Get("mu1_qcd",mu1);
    if(Par.Find("mu2_qcd")!=-1) Par.Get("mu2_qcd",mu2);
    if(Par.Find("mu3_qcd")!=-1) Par.Get("mu3_qcd",mu3);
    if(Par.Find("MBd")!=-1) {
        Par.Get("MBd",m_i);
        setMass(B_D, m_i);
    }
    if(Par.Find("FBd")!=-1) {
        Par.Get("FBd",m_i);
        particles[B_D].setDecayconst(m_i);
    }
    if(Par.Find("BBd")!=-1) {
        Par.Get("BBd",v);
        particles[B_D].setBpars(v);
    }
    if(Par.Find("MBs")!=-1) {
        Par.Get("MBs",m_i);
        setMass(B_S, m_i);
    }
    if(Par.Find("MBp")!=-1) {
        Par.Get("MBp",m_i);
        setMass(B_P, m_i);
    }
    if(Par.Find("MK0")!=-1) {
        Par.Get("MK0",m_i);
        setMass(K_0, m_i);
    }
    if(Par.Find("MKp")!=-1) {
        Par.Get("MKp",m_i);
        setMass(K_P, m_i);
    }

}

QCD::QCD(const Parameters& Par) : Model() {
    for(int i=0; i<NQCDvars; i++)
        if(Par.Find(QCDvars[i])==-1) {
            std::cout << "missing " << QCDvars[i] << "initialization in QCD" << std::endl;
            exit(EXIT_FAILURE);
        }
    update(Par);
}

//QCD::QCD(double AlsM_i, double M_i, double mu_i, double md_i, double ms_i,
//        double mc_i, double mb_i, double mt_i) : Model() {
//    init(AlsM_i, M_i, mu_i, md_i, ms_i, mc_i, mb_i, mt_i, mt_i, mb_i, mc_i);
//}
//
//QCD::QCD(double AlsM_i, double M_i, double mu_i, double md_i, double ms_i,
//        double mc_i, double mb_i, double mt_i, double mu1_i, double mu2_i,
//        double mu3_i) : Model() {
//    init(AlsM_i, M_i, mu_i, md_i, ms_i, mc_i, mb_i, mt_i, mu1_i, mu2_i, mu3_i);
//}

//QCD::QCD(const QCD& orig) : Model() {
//   init(orig.AlsM, orig.M, orig.getMass(UP),
//           orig.getMass(DOWN), orig.getMass(STRANGE),
//           orig.getMass(CHARM), orig.getMass(BOTTOM),
//           orig.getMass(TOP), orig.mu1, orig.mu2, orig.mu3);
//}

QCD::~QCD() {
}

double QCD::beta0(double nf) const {
    return((11.*Nc-2.*nf)/3.);
}

double QCD::beta1(double nf) const {
    return(34./3.*Nc*Nc-10./3.*Nc*nf-(Nc-1./Nc)*nf);
}

double QCD::thresholds(int i) const
{
    if(!(mu1 > mu2 && mu2 > mu3)) {
        std::cout << "inverted thresholds!" << std::endl;
        exit(EXIT_FAILURE);
    }
    switch(i) {
        case 0: return(1.0E10);
        case 1: return(mu1);
        case 2: return(mu2);
        case 3: return(mu3);
        default: return(0.);
    }
}

double QCD::Nf(double mu) const {
    int i;

    for(i=1; i<5; i++)
        if(mu >= thresholds(i))
            return(7.-i);

    exit(EXIT_FAILURE);
}

double QCD::aboveth(double mu) const {
    int i;

    for(i=4; i>=0; i--)
        if(mu < thresholds(i)) return(thresholds(i));

    exit(EXIT_FAILURE);
}

double QCD::belowth(double mu) const {
  int i;

  for(i=0; i<5; i++)
    if(mu >= thresholds(i)) return(thresholds(i));

  exit(EXIT_FAILURE);
}

double QCD::als(double mu, double lam, double nf, int le) const {
  double ll = 2.*log(mu/lam);
  switch(le) {
      case 0:
          return(4.*M_PI/beta0(nf)/ll);
      case 1:
          return(4.*M_PI/beta0(nf)/ll*(1.-beta1(nf)*log(ll)/pow(beta0(nf),2.)/ll));
      case -1:
          return(4.*M_PI/beta0(nf)/ll*(-beta1(nf)*log(ll)/pow(beta0(nf),2.)/ll));
      default:
          std::cout << "als: order not defined\n" << std::endl;
          exit(EXIT_FAILURE);
  }
}

double QCD::als(double mu, double nf, double alsi, double mi, int le) const {
    double v = 1.-beta0(nf)*alsi/2./M_PI*log(mi/mu);
    switch(le) {
        case 0:
            return(alsi/v);
        case 1:
            return(alsi/v*(1.-beta1(nf)/beta0(nf)*alsi/4./M_PI*log(v)/v));
        case -1:
            return(alsi/v*(-beta1(nf)/beta0(nf)*alsi/4./M_PI*log(v)/v));
        default:
            std::cout << "alsf: order not defined\n" << std::endl;
            exit(EXIT_FAILURE);
    }
}

double QCD::als(double mu, double nfmu, int le) const {
    double nfz = Nf(M), m, nfs;

    if(nfmu == nfz) return(als(mu, nfmu, AlsM, M, le));
    if(nfmu > nfz) {
        m = belowth(mu);
        nfs = nfmu-1;
    }
    else {
        m = aboveth(mu);
        nfs = nfmu+1;
    };

    return(als(mu, nfmu, als(m, nfs, le), m, le));
}

void QCD::CacheShift(double cache[][5], int n) const {
    int i,j;
    for(i=4;i>0;i--)
        for(j=0;j<n;j++)
            cache[j][i] = cache[j][i-1];
}

double QCD::als(double mu, int le) const {
    int i;
    for(i=0;i<5;i++)
        if((mu ==  als_cache[0][i]) && (le == als_cache[1][i]) &&
                (AlsM == als_cache[2][i]) && (M == als_cache[3][i]))
            return als_cache[4][i];

    CacheShift(als_cache,5);
    als_cache[0][0] = mu;
    als_cache[1][0] = le;
    als_cache[2][0] = AlsM;
    als_cache[3][0] = M;
    als_cache[4][0] = als(mu, Nf(mu), le);

    return(als_cache[4][0]);
}

double QCD::als(double mu) const {
    return als(mu,1);
}

double QCD::zero(double *x, double *y) const {
    return(als(mu2, *x, 4., (int) *y) - als(mu2, (int) *y));
}

double QCD::lambda4(int le) const {

    if(le!=0 && le!=1){
        std::cout<< "error in order selection in lambda4" << std::endl;
        exit(EXIT_FAILURE);
    }

    int i;
    double alsmu2 = als(mu2,le);
    for(i=0;i<5;i++)
        if(alsmu2 ==  lambda4_cache[0][i])
            return lambda4_cache[1][i];

    CacheShift(lambda4_cache,2);
    lambda4_cache[0][0] = alsmu2;

    double xmin = 0.01, xmax = 0.8;
    TF1 f = TF1("f",this,&QCD::zero,xmin,xmax,1,"QCD","zero");

    ROOT::Math::WrappedTF1 wf1(f);
    double ledouble = le;
    wf1.SetParameters(&ledouble);

    ROOT::Math::BrentRootFinder brf;
    brf.SetFunction( wf1, xmin, xmax );

    if(brf.Solve()) lambda4_cache[1][0]=brf.Root();
    else {
        std::cout << "error in QCD::findroot" << std::endl;
        exit(EXIT_FAILURE);
    }
    return(lambda4_cache[1][0]);
}

// running da m(m) a m(mu)
double QCD::mrun(double mu, double m, double nf, int le) const
{
    double j;

    j = -(4./3.*(4.+97.-10./3.*nf))/2./beta0(nf)+4.*beta1(nf)/pow(beta0(nf),2.);

    switch(le) {
        case 0:
            return(m*pow(als(mu,le)/als(m,le),4./beta0(nf)));
        case 1:
            return(m*pow(als(mu,le)/als(m,le),4./beta0(nf))*(1.+(als(m,le)-
         als(mu,le))/4./M_PI*j));
        case -1:
            return(m*pow(als(mu,le)/als(m,le),4./beta0(nf))*(als(m,le)-
         als(mu,le))/4./M_PI*j);
        default:
            std::cout << "alsf: order not defined\n" << std::endl;
            exit(EXIT_FAILURE);
    }
}

double QCD::mrun(double mu, double m, double nf) const {
    return mrun(mu, m, nf, 1);
}

double QCD::mp2mbara(double * mu, double * mp) const
{
  return(*mp-mbar2mp(*mu));
}

double QCD::mp2mbar(double mp) const {

    int i;
    double ms = getMass(STRANGE), mc = getMass(CHARM);
    double alsmp = als(mp);
    for(i=0;i<5;i++)
        if(alsmp == mp2mbar_cache[0][i] || ms == mp2mbar_cache[1][i] ||
               mc == mp2mbar_cache[2][i] )
            return mp2mbar_cache[3][i];

    CacheShift(mp2mbar_cache,4);
    mp2mbar_cache[0][0] = alsmp;
    mp2mbar_cache[1][0] = ms;
    mp2mbar_cache[2][0] = mc;

    TF1 f("f",this,&QCD::mp2mbara,mp/2.,2.*mp,1,"QCD","mp2mbara");

    ROOT::Math::WrappedTF1 wf1(f);
    wf1.SetParameters(&mp);

    ROOT::Math::BrentRootFinder brf;

    brf.SetFunction( wf1, .7*mp, 1.3*mp );
    if(brf.Solve()) mp2mbar_cache[3][0] = brf.Root();
    else
    {
        std::cout << "error in QCD::mp2mbar" << std::endl;
        exit(EXIT_FAILURE);
    }
    return(mp2mbar_cache[3][0]);
}

double QCD::mbar2mp(double mbar) const {
    if(mbar > 3.)
    {
        double a,D=5.;
        a=als(mbar)/M_PI;
        if(mbar < 50.)
            D=4.-4./3.*(getMass(STRANGE)+getMass(CHARM))/mbar; //only for the b quark

        return(mbar*(1.+4./3.*a+a*a*(13.44434-1.0414*D)));
    }
    else
    {
        std::cout << "can convert only top and bottom masses" << std::endl;
        exit(EXIT_FAILURE);
    }
}

