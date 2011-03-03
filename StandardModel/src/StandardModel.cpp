/* 
 * File:   StandardModel.cpp
 * Author: silvest
 * 
 * Created on November 30, 2010, 1:27 PM
 */


//#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include <iostream>
#include <math.h>
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>
#include "StandardModel.h"
#include <gsl/gsl_sf_dilog.h>

#define LEPS 1.e-5 // tolerance in the limit of S(x,y) to S(x)

//const std::vector<std::string> pino = boost::assign::list_of("A")("BC");
//
//const std::map<std::string,std::vector<std::string> > StandardModel::Deps =
//   boost::assign::map_list_of
////    ("mt",pino)
//    ("mZ",(std::vector<std::string>)boost::assign::list_of("mW")("sin2thw"))
//    ("dAle5Mz",boost::assign::list_of("mW")("sin2thw"))
//    ("alsMz",boost::assign::list_of("mW")("sin2thw"))
//    ("mHl",boost::assign::list_of("mW")("sin2thw"))
//    ("GF",boost::assign::list_of("v")("pino"));

//void StandardModel::init(
//        double mu_i, double md_i, double ms_i, double mc_i, double mb_i,
//        double mt_i, double me_i,
//        double mmu_i, double mtau_i, double mnu1_i, double mnu2_i,
//        double mnu3_i, double GF_i, double alsMz_i, double ale_i, double mZ_i,
//        double dAle5Mz_i, double mHl_i, double mu1_i, double mu2_i,
//        double mu3_i) {
//    QCD::init(alsMz_i, mZ_i, mu_i, md_i, ms_i, mc_i, mb_i, mt_i, mu1_i, mu2_i,
//            mu3_i);
//    setMass(NEUTRINO_1, mnu1_i);
//    setMass(NEUTRINO_2, mnu2_i);
//    setMass(NEUTRINO_3, mnu3_i);
//    setMass(ELECTRON, me_i);
//    setMass(MU, mmu_i);
//    setMass(TAU, mtau_i);
//    GF = GF_i;
//    alsMz = alsMz_i;
//    ale = ale_i;
//    mZ = mZ_i;
//    dAle5Mz = dAle5Mz_i;
//    mHl = mHl_i;
//    Yd.assign(0,0,md_i/v()*sqrt(2.)); // cambiare le scale!!!
//    Yd.assign(1,1,ms_i/v()*sqrt(2.));
//    Yd.assign(2,2,mb_i/v()*sqrt(2.));
//    Yu.assign(0,0,mu_i/v()*sqrt(2.));
//    Yu.assign(1,1,mc_i/v()*sqrt(2.));
//    Yu.assign(2,2,mt_i/v()*sqrt(2.));
//    Yu = Yu*VCKM;
//    Ye.assign(0,0,me_i/v()*sqrt(2.));
//    Ye.assign(1,1,mmu_i/v()*sqrt(2.));
//    Ye.assign(2,2,mtau_i/v()*sqrt(2.));
//    Yn.assign(0,0,mnu1_i/v()*sqrt(2.));
//    Yn.assign(1,1,mnu2_i/v()*sqrt(2.));
//    Yn.assign(2,2,mnu3_i/v()*sqrt(2.));
//    Yn = Yn*UPMNS.hconjugate();
//}

const std::string StandardModel::SMvars[NSMvars] = {"GF","mneutrino_1","mneutrino_2",
"mneutrino_3","melectron","mmu","mtau","VCKM","UPMNS","ale","dAle5Mz","mHl","muw"};

void StandardModel::update(Parameters& Par) {
    QCD::update(Par);
    // within the SM and in the MonteCarlo we force AlsM=alsMz and M=Mz
    alsMz = AlsM;
    mZ = M;
    if(Par.Find("GF")!=-1) Par.Get("GF",GF);
    if(Par.Find("muw")!=-1) Par.Get("muw",muw);
    double m_i;
    bool compute_Yn=false;

    if(Par.Find("mneutrino_1")!=-1) {
        Par.Get("mneutrino_1",m_i);
        setMass(NEUTRINO_1, m_i);
        Yn.assign(0,0,m_i/v()*sqrt(2.));
        compute_Yn=true;
    }
    if(Par.Find("mneutrino_2")!=-1) {
        Par.Get("mneutrino_2",m_i);
        setMass(NEUTRINO_2, m_i);
        Yn.assign(1,1,m_i/v()*sqrt(2.));
        compute_Yn=true;
    }
    if(Par.Find("mneutrino_3")!=-1) {
        Par.Get("mneutrino_3",m_i);
        setMass(NEUTRINO_3, m_i);
        Yn.assign(2,2,m_i/v()*sqrt(2.));
        compute_Yn=true;
    }
    if(Par.Find("melectron")!=-1) {
        Par.Get("melectron",m_i);
        setMass(ELECTRON, m_i);
        Ye.assign(0,0,m_i/v()*sqrt(2.));
    }
    if(Par.Find("mmu")!=-1) {
        Par.Get("mmu",m_i);
        setMass(MU, m_i);
        Ye.assign(1,1,m_i/v()*sqrt(2.));
    }
    if(Par.Find("mtau")!=-1) {
        Par.Get("mtau",m_i);
        setMass(TAU, m_i);
        Ye.assign(2,2,m_i/v()*sqrt(2.));
    }
    if(Par.Find("ale")!=-1) Par.Get("ale",ale);
    if(Par.Find("dAle5Mz")!=-1) Par.Get("dAle5Mz",dAle5Mz);
    if(Par.Find("mHl")!=-1) Par.Get("mHl",mHl);
    for(int i=0; i<3; i++)
    {
        Yu.assign(i,i,this->QCD::particles[quark(UP)+2*i].getMass()/v()*sqrt(2.));
        Yd.assign(i,i,this->QCD::particles[quark(UP)+1+2*i].getMass()/v()*sqrt(2.));
    }
    if(Par.Find("VCKM")!=-1) Par.Get("VCKM",VCKM);
    Yu = Yu*VCKM;
    if(Par.Find("UPMNS")!=-1) {
        Par.Get("UPMNS",UPMNS);
        compute_Yn=true;
    }
    if(compute_Yn) Yn = Yn*UPMNS.hconjugate();
}

//StandardModel::StandardModel(const matrix<complex>& VCKM_i,
//        double mu_i, double md_i, double ms_i, double mc_i, double mb_i,
//        double mt_i, const matrix<complex>& UPMNS_i, double me_i,
//        double mmu_i, double mtau_i, double mnu1_i, double mnu2_i,
//        double mnu3_i, double GF_i, double alsMz_i, double ale_i, double mZ_i,
//        double dAle5Mz_i, double mHl_i, double mu1_i, double mu2_i,
//        double mu3_i) : UPMNS(UPMNS_i), VCKM(VCKM_i), Yd(3,3,0.),
//        Yu(3,3,0.), Ye(3,3,0.), Yn(3,3,0.) {
//    init(mu_i, md_i, ms_i, mc_i, mb_i, mt_i, me_i, mmu_i,
//            mtau_i, mnu1_i, mnu2_i, mnu3_i, GF_i, alsMz_i, ale_i, mZ_i,
//            dAle5Mz_i, mHl_i, mu1_i, mu2_i, mu3_i);
//}
//
//StandardModel::StandardModel(Parameters& Par): UPMNS(3,3,0.), VCKM(3,3,0.),
//        Yd(3,3,0.), Yu(3,3,0.), Ye(3,3,0.), Yn(3,3,0.) {
//    double mu, md, mc, ms, mt, mb, me, mmu, mtau, mnu1, mnu2, mnu3;
//    Par.Get("mu",mu);
//    Par.Get("md",md);
//    Par.Get("mc",mc);
//    Par.Get("ms",ms);
//    Par.Get("mt",mt);
//    Par.Get("mb",mb);
//    Par.Get("me",me);
//    Par.Get("mmu",mmu);
//    Par.Get("mtau",mtau);
//    Par.Get("mnu1",mnu1);
//    Par.Get("mnu2",mnu2);
//    Par.Get("mnu3",mnu3);
//    Par.Get("VCKM",VCKM);
//    Par.Get("UPMNS",UPMNS);
//    Par.Get("GF", GF);
//    Par.Get("alsMz", alsMz);
//    Par.Get("ale", ale);
//    Par.Get("mZ", mZ);
//    Par.Get("dAle5Mz", dAle5Mz);
//    Par.Get("mHl", mHl);
//    if(Par.Find("mu1") == Parameters::DOUBLE) Par.Get("mu1", mu1);
//    else mu1 = mt;
//    if(Par.Find("mu2") == Parameters::DOUBLE) Par.Get("mu2", mu2);
//    else mu2 = mb;
//    if(Par.Find("mu3") == Parameters::DOUBLE) Par.Get("mu3", mu3);
//    else mu3 = mc;
//    init(mu, md, ms, mc, mb, mt, me, mmu, mtau, mnu1,
//            mnu2, mnu3, GF, alsMz, ale, mZ, dAle5Mz, mHl, mu1, mu2, mu3);
//}


//StandardModel::StandardModel(const StandardModel& orig) :
//        UPMNS(orig.getUPMNS()), VCKM(orig.getVCKM()), Yd(3,3,0.), Yu(3,3,0.),
//        Ye(3,3,0.), Yn(3,3,0.) {
//    Init(orig.getMass(UP), orig.getMass(DOWN), orig.getMass(STRANGE),
//            orig.getMass(CHARM), orig.getMass(BOTTOM), orig.getMass(TOP),
//            orig.getMass(ELECTRON), orig.getMass(MU), orig.getMass(TAU),
//            orig.getMass(NEUTRINO_1), orig.getMass(NEUTRINO_2),
//            orig.getMass(NEUTRINO_3), orig.getGF(), orig.getAlsMz(),
//            orig.getAle(), orig.getMZ(), orig.getDAle5Mz(), orig.getMHl(),
//            orig.getMu1(), orig.getMu2(), orig.getMu3());
//}

StandardModel::StandardModel(Parameters& Par) : QCD(Par), UPMNS(3,3,0.),
        VCKM(3,3,0.), Yd(3,3,0.), Yu(3,3,0.), Ye(3,3,0.), Yn(3,3,0.){
    for(int i=0; i<NSMvars; i++)
        if(Par.Find(SMvars[i])==-1) {
            std::cout << "missing " << SMvars[i] << "initialization in SM" << std::endl;
            exit(EXIT_FAILURE);
        }
    update(Par);
}

StandardModel::~StandardModel() {
}

double StandardModel::v() const {
    return 1./sqrt(sqrt(2.)*GF);
}

double StandardModel::mW() const {
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
    double dt = pow((getMass(TOP)/174.3), 2.0) - 1.0;
    double dZ = mZ/91.1875 - 1.0;
    double dalphae = dAle5Mz/0.05907 - 1.0;
    double dalphas = alsMz/0.119 - 1.0;

    return (Mw0 - c1*dH - c2*dH*dH + c3*pow(dH, 4.0)
                + c4*(dh - 1.0) - c5*dalphae + c6*dt - c7*dt*dt
                - c8*dH*dt + c9*dh*dt - c10*dalphas + c11*dZ);
}

double StandardModel::sin2thw() const {
    // Effective leptonic weak mixing angle
    // hep-ph/0407317, hep-ph/0608099
    // applicable for 10 GeV <= mHl <= 1 TeV
    if(mHl<10.||mHl>1000.)
    {
        std::cout << "Higgs mass out of range in sin2thw()" << std::endl;
        exit(EXIT_FAILURE);
    }

    const double s0 = 0.2312527;
    const double d1 = 4.729*0.0001;
    const double d2 = 2.07*0.00001;
    const double d3 = 3.85*0.000001;
    const double d4 = -1.85*0.000001;
    const double d5 = 0.0207;
    const double d6 = -0.002851;
    const double d7 = 1.82*0.0001;
    const double d8 = -9.74*0.000001;
    const double d9 = 3.98*0.0001;
    const double d10 = -0.655;

    // mt, mZ, dAle5Mz and alsMz have to be varied within their combined
    // 2 sigma region around their central values (year 2003) adopted below.
    double L_H = log(mHl/100.0);
    double Delta_H = mHl/100.0;
    double Delta_alphae = dAle5Mz/0.05907 - 1.0;
    double Delta_t = pow((getMass(TOP)/178.0), 2.0) - 1.0;
    double Delta_alphas = alsMz/0.117 - 1.0;
    double Delta_Z = mZ/91.1876 - 1.0;

    return(s0 + d1*L_H + d2*L_H*L_H + d3*L_H*L_H*L_H*L_H
                 + d4*(Delta_H*Delta_H - 1.0) + d5*Delta_alphae + d6*Delta_t
                 + d7*Delta_t*Delta_t + d8*Delta_t*(Delta_H - 1.0)
                 + d9*Delta_alphas + d10*Delta_Z);
}

double StandardModel::sin2thwb() const{
    //Effective mixing angle for b quarks
    //http://arXiv.org/abs/0811.1364v2
    // applicable for 10 GeV <= mHl <= 1 TeV
       const double  s0=0.2327580;
       const double  d1=4.749*0.0001;
       const double  d2=2.03*0.00001;
       const double  d3=3.94*0.000001;
       const double  d4=-1.84*0.000001;
       const double  d5=2.08*0.01;
       const double  d6=-9.93*0.0001;
       const double  d7=7.08*0.00001;
       const double  d8=-7.61*0.000001;
       const double  d9=4.03*0.0001;
       const double  d10=0.661;


    double L_H = log(mHl/100.0);
    double Delta_H = mHl/100.0;
    double Delta_alphae = dAle5Mz/0.05907 - 1.0;
    double Delta_t = pow((getMass(TOP)/178.0), 2.0) - 1.0;
    double Delta_alphas = alsMz/0.117 - 1.0;
    double Delta_Z = mZ/91.1876 - 1.0;

    double sin2b = s0 + d1*L_H + d2*L_H*L_H + d3*L_H*L_H*L_H*L_H
                 + d4*(Delta_H*Delta_H - 1.0) + d5*Delta_alphae + d6*Delta_t
                 + d7*Delta_t*Delta_t + d8*Delta_t*(Delta_H - 1.0)
                 + d9*Delta_alphas + d10*Delta_Z;
    return sin2b;
}


double StandardModel::sin2thwall(const std::string& ferm) const {
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
    double Delta_t = pow((getMass(TOP)/178.0), 2.0) - 1.0;
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



double StandardModel::GammaW() const {
    return 0.0;
}

double StandardModel::GammaZ() const {
    return 0.0;
}

double StandardModel::sigma_had() const {
    return 0.0;
}

double StandardModel::R_l() const {
    return 0.0;
}

double StandardModel::R_c() const {
    return 0.0;
}

double StandardModel::R_b() const {
    return 0.0;
}

double StandardModel::AFB_l() const {
    return 0.0;
}

double StandardModel::AFB_c() const {
    return 0.0;
}

double StandardModel::AFB_b() const {
    return 0.0;
}

double StandardModel::A_l() const {
    return 0.0;
}

double StandardModel::A_c() const {
    return 0.0;
}

double StandardModel::A_b() const {
    return 0.0;
}

double StandardModel::S() const {
    return 0.0;
}

double StandardModel::T() const {
    return 0.0;
}

double StandardModel::U() const {
    return 0.0;
}

// Angles

double StandardModel::getBeta() const {
  return (-VCKM(1,0)*VCKM(1,2).conjugate()/(VCKM(2,0)*VCKM(2,2).conjugate())).arg();
}

double StandardModel::getGamma() const {
  return (-VCKM(0,0)*VCKM(0,2).conjugate()/(VCKM(1,0)*VCKM(1,2).conjugate())).arg();
}

double StandardModel::getAlpha() const {
  return (-VCKM(2,0)*VCKM(2,2).conjugate()/(VCKM(0,0)*VCKM(0,2).conjugate())).arg();
}

double StandardModel::getBetas() const {
  return (-VCKM(2,1)*VCKM(2,2).conjugate()/(VCKM(1,1)*VCKM(1,2).conjugate())).arg();
}

// Lambda_q

gslpp::complex StandardModel::getlamt() const {
  return VCKM(2,0)*VCKM(2,1).conjugate();
}

gslpp::complex StandardModel::getlamc() const {
  return VCKM(1,0)*VCKM(1,1).conjugate();
}

gslpp::complex StandardModel::getlamu() const {
  return VCKM(0,0)*VCKM(0,1).conjugate();
}


gslpp::complex StandardModel::getlamt_d() const {
  return VCKM(2,0)*VCKM(2,2).conjugate();
}

gslpp::complex StandardModel::getlamc_d() const {
  return VCKM(1,0)*VCKM(1,2).conjugate();
}

gslpp::complex StandardModel::getlamu_d() const {
  return VCKM(0,0)*VCKM(0,2).conjugate();
}


gslpp::complex StandardModel::getlamt_s() const {
  return VCKM(2,1)*VCKM(2,2).conjugate();
}

gslpp::complex StandardModel::getlamc_s() const {
  return VCKM(1,1)*VCKM(1,2).conjugate();
}

gslpp::complex StandardModel::getlamu_s() const {
  return VCKM(0,1)*VCKM(0,2).conjugate();
}

double StandardModel::S(double x, double y) const { // Buras 2000 Appendix
  if(fabs(1.-y/x)<LEPS)
    return((x*(-4. + 15.*x - 12.*x*x + pow(x,3.) +
       6.*x*x*log(x)))/(4.*pow(-1. + x,3.)));
  else
    return(x*y*((1./4. + 3./2./(1. - x) - 3./4./pow(1. - x,2.))*
		log(x)/(x - y) +
		(1./4. + 3./2./(1. - y) - 3./4./pow(1. - y,2.))*
		log(y)/(y - x) -
		3./4./(1. - x)/(1.-y)));
}

double StandardModel::kt_sing_a(const double x) const
{
  return(x*(-4.+18.*x+3.*pow(x,2.)+pow(x,3.))/4./pow(x-1.,3.)
         -9.*pow(x,3.)/2./pow(x-1.,4.)*log(x));
}

double StandardModel::kt_sing(const double x) const
{
  return(x*(4.-39.*x+168.*pow(x,2.)+11.*pow(x,3.))/4./pow(x-1.,3.)
         +3.*pow(x,3.)*gsl_sf_dilog(1.-x)*(5.+x)/pow(x-1.,3.)
         +3.*x*log(x)*(-4.+24.*x-36.*pow(x,2.)-7.*pow(x,3.)-pow(x,4.))/2.
         /pow(x-1.,4.)+3.*pow(x,3.)*pow(log(x),2.)*(13.+4.*x+pow(x,2.))/2.
         /pow(x-1.,4.));
}

double StandardModel::kt_oct(const double x) const
{
  return((-64.+68.*x+17.*pow(x,2.)-11.*pow(x,3.))/4./pow(x-1.,2.)
         +pow(M_PI,2.)*8./3./x+2.*gsl_sf_dilog(1.-x)*(8.-24.*x+20.*pow(x,2.)
           -pow(x,3.)+7.*pow(x,4.)-pow(x,5.))/x/pow(x-1.,3.)
         +log(x)*(-32.+68.*x-32.*pow(x,2.)+28.*pow(x,3.)-3.*pow(x,4.))
         /2./pow(x-1.,3.)+pow(x,2.)*pow(log(x),2.)*(4.-7.*x+7.*pow(x,2.)
         -2.*pow(x,3.))/2./pow(x-1.,4.));
}

double StandardModel::kt(const double x, const double mu) const
{
//  printf("%e\n",4./3.*kt_sing(x)+1./3.*kt_oct(x));
  return(4./3.*(kt_sing(x)+12.*log(mu/mW())*kt_sing_a(x))
         +1./3.*(kt_oct(x)+12.*log(mu/mW())*S(x,x))
         +(4.+5./3.)*S(x,x));
}

double StandardModel::eta2bbar(const int LE) const
{
  double eta2b,xt;
  double g0p,g1p5,jp5;

  g0p = 4.;
  g1p5 = -7.+20./9.;
  jp5 = g0p*beta1(5.)/2./pow(beta0(5.),2.)-g1p5/2./beta0(5.);

  xt = pow(mrun(muw,particles[TOP].getMass(),5.,LE)/mW(),2.);

  eta2b = pow(als(muw,LE)/als(particles[BOTTOM].getMass(),LE),
          g0p/2./beta0(5.))*(1.+(LE == 1 ? 1. : 0.)/4./M_PI*
          (als(muw,LE)*(-jp5+kt(xt,muw)/S(xt,xt))+
          als(particles[BOTTOM].getMass(),LE)*jp5));

//  printf("%e %e\n",eta2b*pow(als(Mb_v,le),-g0p/2./beta0(5.))*(1.+le/4./PI*als(Mb_v,le)*jp5),als(Mb_v,le));
  return(eta2b);
}

gslpp::complex StandardModel::getDBD2Amplitude(const int LE) const {
    double xtW = mrun(muw,particles[TOP].getMass(),5.,1)/mW();
    xtW*=xtW;
    return((GF/4./M_PI*mW())*(GF/4./M_PI*mW())*eta2bbar(LE)*getlamt_d()*
          getlamt_d()*S(xtW,xtW)*
          particles[B_D].getDecayconst()*
          particles[B_D].getDecayconst()*particles[B_D].getBpars()[0]*4./3.*
          particles[B_D].getMass());
}
