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

bool QCD::SetFlag(const std::string name , const bool& value){  
    return (false);
}

bool QCD::PreUpdate() {
    computeYu = false;
    computeYd = false;
    computeBd = false;
    computeFBd = false;
    computemt = false;
    
    return (true);
}

bool QCD::PostUpdate(){
    
    if (computeFBd)
        mesons[B_D].setDecayconst(mesons[B_S].getDecayconst() / FBsoFBd);
    if (computeBd)
        BBd.setBpars(0, BBs.getBpars()(0) / BBsoBBd);
    if (computemt)
        quarks[TOP].setMass(Mp2Mbar(mtpole));
    
    return (true);
}

bool QCD::Update(const std::map<std::string, double>& DPars) {
       
   if (!PreUpdate()) return (false);

   UpdateError = false; 
   
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);

    if (UpdateError) return (false);
    
    if (!PostUpdate()) return (false);

    return (true);
}

void QCD::SetParameter(const std::string name, const double& value) {
    if(name.compare("AlsMz")==0){
        AlsMz = value;
        computemt = true;
    }
    else if(name.compare("Mz")==0)
        Mz = value;
    else if(name.compare("mup")==0){
        if(value < MEPS) UpdateError = true; 
        quarks[UP].setMass(value);
        computeYu = true;
    }
    else if(name.compare("mdown")==0){
        if(value < MEPS) UpdateError = true;
        quarks[DOWN].setMass(value);
        computeYd = true;
    }
    else if(name.compare("mcharm")==0){
        quarks[CHARM].setMass(value);
        computeYu = true;
    }
    else if(name.compare("mstrange")==0){
        if(value < MEPS) UpdateError = true;
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

////////////////////////////////////////////////////////////////////////

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

double QCD::Nf(double mu) const {
    int i;
    for(i=1; i<5; i++)
        if(mu >= Thresholds(i))
            return (7. - (double)i);
    std::cout << "error in QCD::Nf" << std::endl;
    exit(EXIT_FAILURE);
}

void QCD::CacheShift(double cache[][5], int n) const {
    int i,j;
    for(i=4;i>0;i--)
        for(j=0;j<n;j++)
            cache[j][i] = cache[j][i-1];
}

////////////////////////////////////////////////////////////////////////

double QCD::Beta0(double nf) const {
    return ((11.*Nc-2.*nf)/3.);
}

double QCD::Beta1(double nf) const {
    return (34./3.*Nc*Nc-10./3.*Nc*nf-(Nc-1./Nc)*nf);
}

double QCD::Beta2(double nf) const{
    return (2857./2.-5033./18.*nf+325./54.*nf*nf);
}

double QCD::AlsWithLambda(double mu, double logLambda, double nf, orders order) const {
    double L = 2.*(log(mu) - logLambda);
    double log_L, b0=Beta0(nf), b1, b2;
    double b0L = b0*L; 
    if (order!=LO) { 
        log_L = log(L);
        b1 = Beta1(nf);
    }
    if (order!=LO && order!=FULLNLO && order != NLO) b2 = Beta2(nf);
    
    switch (order) {
        case LO:
            return ( 4.*M_PI/b0L );
        case FULLNLO:
            return ( 4.*M_PI/b0L*(1. - b1*log_L/b0/b0L) );
        case NLO:
            return ( 4.*M_PI/b0L*(   - b1*log_L/b0/b0L) );
        case FULLNNLO:
            return ( 4.*M_PI/b0L
                     *(1. - b1*log_L/b0/b0L
                         + 1./b0L/b0L
                           *(b1*b1/b0/b0*(log_L*log_L-log_L-1.)+b2/b0)) );
        case NNLO:
            return ( 4.*M_PI/b0L
                     *(     1./b0L/b0L
                            *(b1*b1/b0/b0*(log_L*log_L-log_L-1.)+b2/b0)) );
        default:
           std::cerr << "QCD::AlsFromLambda(): order " << order 
                     << " not defined" << std::endl;
           exit(EXIT_FAILURE);
    }
}

double QCD::Als(double mu, double nf, double alsi, double mi, orders order) const {
    double v = 1.-Beta0(nf)*alsi/2./M_PI*log(mi/mu);
    switch (order) {
        case LO:
            return (alsi/v);
        case FULLNLO:
            return (alsi/v*(1.-Beta1(nf)/Beta0(nf)*alsi/4./M_PI*log(v)/v));
        case NLO:
            return (alsi/v*(-Beta1(nf)/Beta0(nf)*alsi/4./M_PI*log(v)/v));
        default:
           std::cerr << "QCD::Als(): order " << order 
                     << " not defined" << std::endl;
           exit(EXIT_FAILURE);
    }
}

double QCD::Als(double mu, double nfmu, orders order) const {
    double nfz = Nf(Mz), m, nfs;

    // Note: Threshold correction should be taken into account 
    //       in the cases of FULLNLO/NLO, when using the function 
    //       Als(mu, nfmu, Als(m, nfs, order), m, order). 
    
    switch (order) {
        case LO:
        case FULLNLO:
        case NLO:
            if(nfmu == nfz) 
                return Als(mu, nfmu, AlsMz, Mz, order);
            else if(nfmu > nfz) {
                m = BelowTh(mu);
                nfs = nfmu - 1.;
            } else {
                m = AboveTh(mu);
                nfs = nfmu + 1.;
            }
            return Als(mu, nfmu, Als(m, nfs, order), m, order);
        case FULLNNLO:
        case NNLO:
            return AlsWithLambda(mu, logLambda(mu, FULLNNLO), nfmu, order);
        default:
            throw "Error in QCD::Als(mu,nfmu,order)";
    }
}

double QCD::Als(double mu, orders order) const {
    return Als(mu, Nf(mu), order);

    //[Caching]
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
    return ( AlsWithLambda(Mz, *x, 5., (orders) *y) - AlsMz );
}

double QCD::logLambda5(orders order) const {
    if (order==NLO) order = FULLNLO;
    if (order==NNLO) order = FULLNNLO;
    
    for (int i=0; i<5; ++i)
        if ( (AlsMz == logLambda5_cache[0][i]) 
              && (Mz == logLambda5_cache[1][i]) 
              && ((double)order == logLambda5_cache[2][i]) )
            return logLambda5_cache[3][i];

    CacheShift(logLambda5_cache,4);
    logLambda5_cache[0][0] = AlsMz;
    logLambda5_cache[1][0] = Mz;
    logLambda5_cache[2][0] = (double)order;
    
    if (order==LO)
        logLambda5_cache[3][0] = log(Mz) - 2.*M_PI/Beta0(5.)/AlsMz;
    else {
        double xmin = -4., xmax = -0.2;
        TF1 f = TF1("f",this,&QCD::ZeroNf5,xmin,xmax,1,"QCD","zeroNf5");
        
        ROOT::Math::WrappedTF1 wf1(f);
        double ledouble = (double)order;
        wf1.SetParameters(&ledouble);
        
        ROOT::Math::BrentRootFinder brf;
        brf.SetFunction( wf1, xmin, xmax );
        
        if (brf.Solve()) logLambda5_cache[3][0] = brf.Root();
        else
            throw "Error in QCD::logLambda5()";
    }
    return ( logLambda5_cache[3][0] );
}

double QCD::logLambda(double muMatching, double mf, double nfNEW, double nfORG, 
                      double logLambdaORG, orders order) const {
    if (fabs(nfNEW-nfORG)!=1.)
        throw "Error in QCD::logLambda()";
    if (order==NLO) order = FULLNLO;
    if (order==NNLO) order = FULLNNLO;
    
    double rNEW = Beta1(nfNEW)/Beta0(nfNEW);
    double rORG = Beta1(nfORG)/Beta0(nfORG);
    double logMuMatching = log(muMatching);    
    double L = 2.*(logMuMatching - logLambdaORG);
    double log_L, log_mu2_mf2;
    double C1, C2; // threshold corrections        
    double logLambdaNEW; 
    
    // LO contribution
    logLambdaNEW = 1./2./Beta0(nfNEW)
                   *(Beta0(nfNEW) - Beta0(nfORG))*L + logLambdaORG;

    // NLO contribution
    if (order==FULLNLO || order==FULLNNLO) {
        log_mu2_mf2 = 2.*(logMuMatching - log(mf));
        log_L = log(L);
        if (nfNEW < nfORG)
            C1 = 2./3.*log_mu2_mf2;
        else
            C1 = -2./3.*log_mu2_mf2;
        logLambdaNEW += 1./2./Beta0(nfNEW)
                        *( (rNEW - rORG)*log_L
                           - rNEW*log(Beta0(nfNEW)/Beta0(nfORG)) - C1 );
    }

    // NNLO contribution
    if (order==FULLNNLO) {
        if (nfNEW < nfORG) {
            C2 = -16.*(log_mu2_mf2*log_mu2_mf2/36. - 19./24.*log_mu2_mf2 + 11./72.);                
        } else {
            C2 = -16.*(log_mu2_mf2*log_mu2_mf2/36. + 19./24.*log_mu2_mf2 - 11./72.);
        }
        logLambdaNEW += 1./2./Beta0(nfNEW)/Beta0(nfORG)/L
                        *( rORG*(rNEW - rORG)*log_L + rNEW*rNEW - rORG*rORG 
                          - Beta2(nfNEW)/Beta0(nfNEW) + Beta2(nfORG)/Beta0(nfORG) 
                          + rNEW*C1 - C1*C1 - C2 );
    }
    
    return logLambdaNEW;
}

double QCD::logLambda(double mu, orders order) const {
    if (order==NLO) order = FULLNLO;
    if (order==NNLO) order = FULLNNLO;
    
    double muMatching, mf, logLambdaORG, logLambdaNEW;
    if (Nf(mu)==5.) 
        return logLambda5(order);
    else if (Nf(mu)==6.) {
        muMatching = Thresholds(1); // mut

        //!!!! NNLO codes should be implemented !!!!
        mf = getQuarks(TOP).getMass(); // m_t(m_t)

        return logLambda(muMatching, mf, 6., 5., logLambda5(order), order);
    } else if (Nf(mu)==4. || Nf(mu)==3.) { 
        muMatching = Thresholds(2); // mub
        mf = getQuarks(BOTTOM).getMass(); // m_b(m_b)
        logLambdaORG = logLambda5(order);
        logLambdaNEW = logLambda(muMatching, mf, 4., 5., logLambdaORG, order);
        if (Nf(mu)==3.) { 
            muMatching = Thresholds(3); // muc
            mf = getQuarks(CHARM).getMass(); // m_c(m_c)
            logLambdaORG = logLambdaNEW;
            logLambdaNEW = logLambda(muMatching, mf, 3., 4., logLambdaORG, order);
        }
        return logLambdaNEW;
    } else
        throw "Error in QCD::logLambda()";
}

////////////////////////////////////////////////////////////////////////

// running da m(m) a m(mu)
double QCD::Mrun(double mu, double m, orders order) const {
    return(Mrun(mu,m,m,order));
}

// running da m(mu_i) a m(mu_f)
double QCD::Mrun(double mu_f, double mu_i, double m, orders order) const {
    int i;
//    if(fabs(mu_i - m) < MEPS){
//        m = mu_i;
//    }
    for (i = 0; i < 5; i++) {
        if ((mu_f == mrun_cache[0][i]) && (mu_i == mrun_cache[1][i]) &&
                (m == mrun_cache[2][i]) && (order == mrun_cache[3][i]) &&
                (AlsMz == mrun_cache[4][i]))
            return mrun_cache[5][i];
    }
    

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


double QCD::MS2DRqmass(const double & MSscale, const double& MSbar) const {
    return(MSbar/(1.+Als(MSscale)/4./M_PI*CF));
}
