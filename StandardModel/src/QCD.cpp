/* 
 * File:   QCD.cpp
 * Author: marco
 * 
 * Created on February 17, 2011, 2:13 PM
 */

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <map>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>
#include "QCD.h"

const std::string QCD::QCDvars[NQCDvars] = {
    "AlsMz","Mz","mup","mdown","mcharm","mstrange",
    "mtop","mbottom","mut","mub","muc","MBd",
    "MBs","MBp","MK0","MKp","MD","FBs","FBsoFBd", "FD",
    "BBsoBBd","BBs1","BBs2","BBs3","BBs4","BBs5", "BBsscale", "BBsscheme",
    "BD1","BD2","BD3","BD4","BD5", "BDscale", "BDscheme",
    "BK1","BK2","BK3","BK4","BK5", "BKscale", "BKscheme"
};

void QCD::Update(const std::map<std::string, double>& DPars) {
    computeYu = false;
    computeYd = false;
    computeBd = false; 
    computeFBd = false;
    computemt = false;
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++) 
        SetParameter(it->first,it->second);
    if(computeFBd)
        mesons[B_D].setDecayconst(mesons[B_S].getDecayconst()/FBsoFBd);
    if(computeBd)
        BBd.setBpars(0, BBs.getBpars()(0)/BBsoBBd);
    if(computemt)
        quarks[TOP].setMass(Mp2Mbar(mtpole));
}

void QCD::SetParameter(const std::string name, const double& value) {
    if(name.compare("AlsMz")==0){
        AlsMz = value;
        computemt = true;
    }
    else if(name.compare("Mz")==0)
        Mz = value;
    else if(name.compare("mup")==0){
        quarks[UP].setMass(value);
        computeYu = true;
    }
    else if(name.compare("mdown")==0){
        quarks[DOWN].setMass(value);
        computeYd = true;
    }
    else if(name.compare("mcharm")==0){
        quarks[CHARM].setMass(value);
        computeYu = true;
    }
    else if(name.compare("mstrange")==0){
        quarks[STRANGE].setMass(value);
        computeYd = true;
    }
    else if(name.compare("mtop")==0){
        mtpole = value;
        computeYu = true;
        computemt = true;
    }
    else if(name.compare("mbottom")==0){
        quarks[BOTTOM].setMass(value);
        computeYd = true;
    }
    else if(name.compare("mut")==0)
        mut = value;
    else if(name.compare("mub")==0)
        mub = value;
    else if(name.compare("muc")==0)
        muc = value;
    else if(name.compare("MBd")==0)
        mesons[B_D].setMass(value);
    else if(name.compare("MBs")==0)
        mesons[B_S].setMass(value);
    else if(name.compare("MBp")==0)
        mesons[B_P].setMass(value);
    else if(name.compare("MK0")==0)
        mesons[K_0].setMass(value);
    else if(name.compare("MKp")==0)
        mesons[K_P].setMass(value);
    else if(name.compare("MD")==0)
        mesons[D_0].setMass(value);
    else if(name.compare("FBs")==0) {
        mesons[B_S].setDecayconst(value);
        computeFBd = true;
    }
    else if(name.compare("FBsoFBd")==0) {
        FBsoFBd = value;
        computeFBd = true;
    }
    else if(name.compare("FD")==0) {
        mesons[D_0].setDecayconst(value);
    }
    else if(name.compare("BBsoBBd")==0) {
        BBsoBBd = value;
        computeBd = true;
    }
    else if(name.compare("BBs1")==0) {
        BBs.setBpars(0,value);
        computeBd = true;
    }
   else if(name.compare("BBs2")==0) {
        BBd.setBpars(1,value);
        BBs.setBpars(1,value);
    }
    else if(name.compare("BBs3")==0) {
        BBd.setBpars(2,value);
        BBs.setBpars(2,value);
    }
    else if(name.compare("BBs4")==0){
        BBd.setBpars(3,value);
        BBs.setBpars(3,value);
    }
    else if(name.compare("BBs5")==0){
        BBd.setBpars(4,value);
        BBs.setBpars(4,value);
    }    
    else if(name.compare("BBsscale")==0){
        BBd.setMu(value);
        BBs.setMu(value);
    }    
    else if(name.compare("BBsscheme")==0){
        BBd.setScheme((schemes) value);
        BBs.setScheme((schemes) value);
    }
    else if(name.compare("BD1")==0) {
        BD.setBpars(0,value);
    }
    else if(name.compare("BD2")==0) {
        BD.setBpars(1,value);
    }
    else if(name.compare("BD3")==0) {
        BD.setBpars(2,value);
    }
    else if(name.compare("BD4")==0){
        BD.setBpars(3,value);}
    else if(name.compare("BD5")==0){
        BD.setBpars(4,value);
    }    
    else if(name.compare("BDscale")==0){
        BD.setMu(value);
    }    
    else if(name.compare("BDscheme")==0){
        BD.setScheme((schemes) value);
    }
    else if(name.compare("BK1")==0) {
        BK.setBpars(0,value);
    }
   else if(name.compare("BK2")==0) {
        BK.setBpars(1,value);
   }
    else if(name.compare("BK3")==0) {
        BK.setBpars(2,value);
    }
    else if(name.compare("BK4")==0){
        BK.setBpars(3,value);}
    else if(name.compare("BK5")==0){
        BK.setBpars(4,value);
    }    
    else if(name.compare("BKscale")==0){
        BK.setMu(value);
    }    
    else if(name.compare("BKscheme")==0){
        BK.setScheme((schemes) value);
    }
//    else {
//        std::cout << "cannot set parameter " << name << " in SetQCDParameter" << std::endl;
//        exit(EXIT_FAILURE);
//    }

}

bool QCD::CheckParameters(const std::map<std::string, double>& DPars) {
    for(int i=0;i<NQCDvars;i++)
        if(DPars.find(QCDvars[i])==DPars.end()) {
            std::cout << "missing mandatory QCD parameter " << QCDvars[i] << std::endl;
            return false;
        }
    return true;
}

bool QCD::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return(CheckParameters(DPars));
}

double QCD::Beta0(double nf) const {
    return((11.*Nc-2.*nf)/3.);
}

double QCD::Beta1(double nf) const {
    return(34./3.*Nc*Nc-10./3.*Nc*nf-(Nc-1./Nc)*nf);
}

double QCD::Beta2(double nf) const{
    return (2857./2.-5033./18.*nf+325./54.*nf*nf);
}

double QCD::Thresholds(int i) const {
    if(!(mut > mub && mub > muc)) {
        std::cout << "inverted thresholds!" << std::endl;
        exit(EXIT_FAILURE);
    }
    switch(i) {
        case 0: return(1.0E10);
        case 1: return(mut);
        case 2: return(mub);
        case 3: return(muc);
        default: return(0.);
    }
}

double QCD::Nf(double mu) const {
    int i;

    for(i=1; i<5; i++)
        if(mu >= Thresholds(i))
            return(7.-i);
    std::cout << "error in QCD::Nf" << std::endl;
    exit(EXIT_FAILURE);
}

double QCD::AboveTh(double mu) const {
    int i;

    for(i=4; i>=0; i--)
        if(mu < Thresholds(i)) return(Thresholds(i));

    exit(EXIT_FAILURE);
}

double QCD::BelowTh(double mu) const {
  int i;

  for(i=0; i<5; i++)
    if(mu >= Thresholds(i)) return(Thresholds(i));

  exit(EXIT_FAILURE);
}

double QCD::Als(double mu, double lam, double nf, orders order) const {
  double ll = 2.*log(mu/lam);
  double log_ll;
  if (order!=LO)
      log_ll = log(ll);
  
  switch(order) {
      case LO:
          return(4.*M_PI/Beta0(nf)/ll);
      case FULLNLO:
          return(4.*M_PI/Beta0(nf)/ll*(1.-Beta1(nf)*log_ll/pow(Beta0(nf),2.)/ll));
      case NLO:
          return(4.*M_PI/Beta0(nf)/ll*(-Beta1(nf)*log_ll/pow(Beta0(nf),2.)/ll));
      case FULLNNLO:
          return(4.*M_PI*(1./Beta0(nf)/ll-Beta1(nf)*log_ll/pow(Beta0(nf)*ll,2.)/Beta0(nf)
                          + 1./pow(Beta0(nf)*ll, 3.)
                            *( Beta1(nf)*Beta1(nf)/Beta0(nf)/Beta0(nf)
                               *(log_ll*log_ll-log_ll-1.) + Beta2(nf)/Beta0(nf) )));
      case NNLO:
          return(4.*M_PI*(1./pow(Beta0(nf)*ll, 3.)
                         *( Beta1(nf)*Beta1(nf)/Beta0(nf)/Beta0(nf)
                            *(log_ll*log_ll-log_ll-1.) + Beta2(nf)/Beta0(nf) )));
      default:
          std::cerr << "als: order " << order <<" not defined\n" << std::endl;
          exit(EXIT_FAILURE);
  }
}

double QCD::Als(double mu, double nf, double alsi, double mi, orders order) const {
    double v = 1.-Beta0(nf)*alsi/2./M_PI*log(mi/mu);
    switch(order) {
        case LO:
            return(alsi/v);
        case FULLNLO:
            return(alsi/v*(1.-Beta1(nf)/Beta0(nf)*alsi/4./M_PI*log(v)/v));
        case NLO:
            return(alsi/v*(-Beta1(nf)/Beta0(nf)*alsi/4./M_PI*log(v)/v));
        default:
           std::cerr << "als: order " << order <<" not defined\n" << std::endl;
           exit(EXIT_FAILURE);
    }
}

double QCD::Als(double mu, double nfmu, orders order) const {
    double nfz = Nf(Mz), m, nfs;

    // Note: Threshold correction should be taken into account 
    //       in the cases of FULLNLO/NLO.
    
    switch (order) {
        case LO:
        case FULLNLO:
        case NLO:
            if(nfmu == nfz) 
                return(Als(mu, nfmu, AlsMz, Mz, order));
            else if(nfmu > nfz) {
                m = BelowTh(mu);
                nfs = nfmu - 1.;
            } else {
                m = AboveTh(mu);
                nfs = nfmu + 1.;
            }
            return(Als(mu, nfmu, Als(m, nfs, order), m, order));  
        case FULLNNLO:
        case NNLO:
            return(Als(mu, Lambda(mu, FULLNNLO), nfmu, order));
        default:
            throw "Error in QCD::Als(mu,nfmu,order)";
    }
}

void QCD::CacheShift(double cache[][5], int n) const {
    int i,j;
    for(i=4;i>0;i--)
        for(j=0;j<n;j++)
            cache[j][i] = cache[j][i-1];
}

double QCD::Als(double mu, orders order) const {
    return Als(mu, Nf(mu), order);

// Note: Threshold scales should be checked in addition to mu, order, AlsMz and Mz.
//    int i;
//    for(i=0;i<5;i++)
//        if((mu ==  als_cache[0][i]) && (order == als_cache[1][i]) &&
//                (AlsMz == als_cache[2][i]) && (Mz == als_cache[3][i]))
//            return als_cache[4][i];
//
//    CacheShift(als_cache,5);
//    als_cache[0][0] = mu;
//    als_cache[1][0] = order;
//    als_cache[2][0] = AlsMz;
//    als_cache[3][0] = Mz;
//    als_cache[4][0] = Als(mu, Nf(mu), order);
//
//    return(als_cache[4][0]);
}

double QCD::ZeroNf5(double *x, double *y)const{
    return ( Als(Mz, *x, 5., (orders) *y) - AlsMz );
}

double QCD::Lambda5(orders order) const {
    if(order!=LO && order!=FULLNLO && order!=FULLNNLO)
        throw "Error in order selection in QCD::Lambda5()";

    for(int i=0; i<5; i++)
        if( (AlsMz ==  lambda5_cache[0][i]) && (Mz == lambda5_cache[1][i]) 
             && ((double)order == lambda5_cache[2][i]) )
            return lambda5_cache[3][i];

    CacheShift(lambda5_cache,4);
    lambda5_cache[0][0] = AlsMz;
    lambda5_cache[1][0] = Mz;
    lambda5_cache[2][0] = (double)order;
    
    //------ TEST for NNLO -------------------
    //double beta05 = Beta0(5.);
    //double beta15 = Beta1(5.);
    //double beta25 = Beta2(5.);
    //double LogLambda5 
    //         = log(Mz)  
    //           - 1./2./beta05*( 4.*M_PI/AlsMz + beta15/beta05*log(AlsMz/M_PI) 
    //                           + (beta25/beta05 - pow(beta15/beta05,2.))*AlsMz/4./M_PI )
    //           - beta15/2./beta05/beta05*log(beta05/4.);
    //lambda5_cache[3][0] = exp(LogLambda5);
    //std::cout << "Als(Mz)= " 
    //          << Als(Mz, lambda5_cache[3][0], 5., FULLNNLO) << std::endl;// TEST
    //return ( lambda5_cache[3][0] );
    //----------------------------------------
    
    double xmin = 0.01, xmax = 0.8;
    TF1 f = TF1("f",this,&QCD::ZeroNf5,xmin,xmax,1,"QCD","zeroNf5");

    ROOT::Math::WrappedTF1 wf1(f);
    double ledouble = (double)order;
    wf1.SetParameters(&ledouble);

    ROOT::Math::BrentRootFinder brf;
    brf.SetFunction( wf1, xmin, xmax );

    if(brf.Solve()) lambda5_cache[3][0]=brf.Root();
    else
        throw "Error in QCD::Lambda5()";

    return ( lambda5_cache[3][0] );
}

double QCD::Lambda(double mu, orders order) const {
    if(order!=LO && order!=FULLNLO && order!=FULLNNLO)
        throw "Error in QCD::Lambda()";
    
    double muMatching, mf, LambdaORG, LambdaTMP;
    if (Nf(mu)==5.) 
        LambdaTMP = Lambda5(order);
    else if (Nf(mu)==6.) {
        muMatching = Thresholds(1);
        mf = getQuarks(TOP).getMass(); // NNLO codes should be implemented!!!
        LambdaORG = Lambda5(order);
        LambdaTMP = Lambda(muMatching, mf, 6., 5., LambdaORG, order);
    } else if (Nf(mu)==4. || Nf(mu)==3.) { 
        muMatching = Thresholds(2);
        mf = getQuarks(BOTTOM).getMass();
        LambdaORG = Lambda5(order);
        LambdaTMP = Lambda(muMatching, mf, 4., 5., LambdaORG, order);
        if (Nf(mu)==3.) { 
            muMatching = Thresholds(3);
            mf = getQuarks(CHARM).getMass();
            LambdaORG = LambdaTMP;
            LambdaTMP = Lambda(muMatching, mf, 3., 4., LambdaORG, order);
        }
    } else
        throw "Error in QCD::Lambda()";

    return LambdaTMP;
}

double QCD::Lambda(double muMatching, double mf, double nfNEW, double nfORG, 
                   double LambdaORG, orders order) const {
    if (fabs(nfNEW-nfORG)!=1.)
        throw "Error in QCD::Lambda()";
    
    double rNEW = Beta1(nfNEW)/Beta0(nfNEW);
    double rORG = Beta1(nfORG)/Beta0(nfORG);    
    double Lf = log(muMatching*muMatching/LambdaORG/LambdaORG);
    double log_Lf, log_mu2_mf2;
    double C1, C2; // threshold corrections        
    double LogLambda;    
    switch (order) {
        case LO:
            LogLambda = 1./2./Beta0(nfNEW)
                         *(Beta0(nfNEW) - Beta0(nfORG))*Lf + log(LambdaORG);
            break;
        case FULLNLO:
            log_mu2_mf2 = log(muMatching*muMatching/mf/mf);
            C1 = 2./3.*log_mu2_mf2;
            if (nfNEW > nfORG)
                C1 *= - 1.; 
            
            LogLambda = 1./2./Beta0(nfNEW)
                        *( (Beta0(nfNEW) - Beta0(nfORG))*Lf 
                            + (rNEW - rORG)*log(Lf)
                            - rNEW*log(Beta0(nfNEW)/Beta0(nfORG)) - C1 ) 
                        + log(LambdaORG);
            break;
        case FULLNNLO:
            log_mu2_mf2 = log(muMatching*muMatching/mf/mf);
            log_Lf = log(Lf);
            C1 = 2./3.*log_mu2_mf2;
            C2 = - 16.*(log_mu2_mf2*log_mu2_mf2/36. - 19./24.*log_mu2_mf2 + 11./72.);                
            if (nfNEW > nfORG) {
                C1 *= - 1.; 
                C2 *= -16.*(log_mu2_mf2*log_mu2_mf2/36. + 19./24.*log_mu2_mf2 - 11./72.);
            }
            
            LogLambda = 1./2./Beta0(nfNEW)
                        *( (Beta0(nfNEW) - Beta0(nfORG))*Lf 
                            + (rNEW - rORG)*log_Lf 
                            - rNEW*log(Beta0(nfNEW)/Beta0(nfORG)) - C1
                            + 1./Beta0(nfORG)/Lf
                              *( rORG*(rNEW - rORG)*log_Lf
                                 + rNEW*rNEW - rORG*rORG 
                                 - Beta2(nfNEW)/Beta0(nfNEW) 
                                 + Beta2(nfORG)/Beta0(nfORG) 
                               + rNEW*C1 - C1*C1 - C2 ) ) 
                        + log(LambdaORG);
            break;
        default:
            throw "Error in QCD::Lambda()";
    }
    
    return ( exp( LogLambda ) );
}

// running da m(m) a m(mu)
double QCD::Mrun(double mu, double m, orders order) const {
    return(Mrun(mu,m,m,order));
}

// running da m(mu_i) a m(mu_f)
double QCD::Mrun(double mu_f, double mu_i, double m, orders order) const {
    int i;
    for(i=0;i<5;i++)
        if((mu_f ==  mrun_cache[0][i]) && (mu_i == mrun_cache[1][i]) &&
                (m == mrun_cache[2][i]) && (order == mrun_cache[3][i])&&
                (AlsMz == mrun_cache[4][i]))
            return mrun_cache[5][i];

    double mx, nfx;
    double nfi = Nf(mu_i), nff = Nf(mu_f);

    if(nff == nfi || (nff == nfi - 1 && mu_i == AboveTh(mu_f))) 
            return(Mrun(mu_f, mu_i, m, nff, order));
    if(nff > nfi) {
        mx = AboveTh(mu_i);
        nfx = nfi;
    }
    else {
        mx = BelowTh(mu_i);
        nfx = nfi-1;
    }

    double mrun = Mrun(mu_f, mx, Mrun(mx, mu_i, m, nfx, order), order);

    CacheShift(mrun_cache,6);
    mrun_cache[0][0] = mu_f;
    mrun_cache[1][0] = mu_i;
    mrun_cache[2][0] = m;
    mrun_cache[3][0] = order;
    mrun_cache[4][0] = AlsMz;
    mrun_cache[5][0] = mrun;
    
    return(mrun);
}

// running da m(mu_i) a m(mu_f) a nf fissato
double QCD::Mrun(double mu_f, double mu_i, double m, double nf, orders order) const {
    double j;

    j = -(4./3.*(4.+97.-10./3.*nf))/2./Beta0(nf)+4.*Beta1(nf)/(Beta0(nf)*Beta0(nf));

    switch(order) {
        case LO:
            return(m*pow(Als(mu_f,order)/Als(mu_i,order),4./Beta0(nf)));
        case FULLNLO:
            return(m*pow(Als(mu_f,order)/Als(mu_i,order),4./Beta0(nf))*(1.+(Als(mu_i,order)-
                    Als(mu_f,order))/4./M_PI*j));
        case NLO:
            return(m*pow(Als(mu_f,order)/Als(mu_i,order),4./Beta0(nf))*(Als(mu_i,order)-
                    Als(mu_f,order))/4./M_PI*j);
        default:
            std::cerr << "mrun: order not defined\n" << std::endl;
            exit(EXIT_FAILURE);
    }
}

double QCD::Mp2Mbara(double * mu, double * mp) const {
    return(*mp-Mbar2Mp(*mu));
}

double QCD::Mp2Mbar(double mp) const {

    int i;
    double ms = quarks[STRANGE].getMass(), mc = quarks[CHARM].getMass();
    double alsmp = Als(mp);
    for(i=0;i<5;i++)
        if(alsmp == mp2mbar_cache[0][i] && ms == mp2mbar_cache[1][i] &&
               mc == mp2mbar_cache[2][i] )
            return mp2mbar_cache[3][i];

    CacheShift(mp2mbar_cache,4);
    mp2mbar_cache[0][0] = alsmp;
    mp2mbar_cache[1][0] = ms;
    mp2mbar_cache[2][0] = mc;

    TF1 f("f",this,&QCD::Mp2Mbara,mp/2.,2.*mp,1,"QCD","mp2mbara");

    ROOT::Math::WrappedTF1 wf1(f);
    wf1.SetParameters(&mp);

    ROOT::Math::BrentRootFinder brf;

    brf.SetFunction( wf1, .7*mp, 1.3*mp );
    if(brf.Solve()) mp2mbar_cache[3][0] = brf.Root();
    else
    {
        std::cerr << "error in QCD::mp2mbar" << std::endl;
        exit(EXIT_FAILURE);
    }
    return(mp2mbar_cache[3][0]);
}

double QCD::Mbar2Mp(double mbar) const {
    if(mbar > 3.)
    {
        double a,D=5.;
        a=Als(mbar)/M_PI;
        if(mbar < 50.)
            D=4.-4./3.*(quarks[STRANGE].getMass()+quarks[CHARM].getMass())/mbar; //only for the b quark

        return(mbar*(1.+4./3.*a+a*a*(13.44434-1.0414*D)));
    }
    else
    {
        std::cerr << "can convert only top and bottom masses" << std::endl;
        exit(EXIT_FAILURE);
    }
}

double QCD::MS2DRqmass(const double& MSbar) const {
    return(MSbar/(1.+Als(MSbar)/4./M_PI*CF));
}