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
#include <gsl/gsl_sf_zeta.h>

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
            std::cout << "missing " << SMvars[i] << " initialization in SM" << std::endl;
            exit(EXIT_FAILURE);
        }
    update(Par);
}

StandardModel::~StandardModel() {
}


///////////////////////////////////////////////////////////////////////////

double StandardModel::v() const {
    return 1./sqrt(sqrt(2.)*GF);
}


////////////////////////////////////////////////////////////////////////

double StandardModel::dAleLepMz() const {
    /*
     * Eqs.(14-15) in J.H.Kuhn, M.Steinhauser, PLB437,425(1998) [hep-ph/9802241]
     * Eqs.(5-10) in M.Steinhauser, PLB429,158(1998) [hep-ph/9803313]
     *
     * oneLoop=314.19007, twoLoop=0.77617, threeLoop=0.01063
     * sum=314.97686
     *
     * Notes: oneLoop and twoLoop are OK for me=0.00051099907, mmu=0.105658389,
     *        mtau=1.777, ale=1.0/137.0359895 and mZ=91.187, but only threeLoop
     *        differs from the above value by 5% (Why?).
     */
    double xl[3] = {mZ*mZ/getMass(ELECTRON)/getMass(ELECTRON),
                    mZ*mZ/getMass(MU)/getMass(MU),
                    mZ*mZ/getMass(TAU)/getMass(TAU)};
    double log_l[3] = {log(xl[0]), log(xl[1]), log(xl[2])};
    double log2 = log(2.0);

    /* TESTS */
    //std::cout << "Me= " << particles[ELECTRON].getMass() << std::endl; // WRONG!
    //std::cout << "Me= " << getMass(ELECTRON) << std::endl;             // OK!

    /* zeta functions */
    double zeta2 = gsl_sf_zeta_int(2);
    double zeta3 = gsl_sf_zeta_int(3);
    double zeta5 = gsl_sf_zeta_int(5);

    double oneLoop[3], twoLoop[3], threeLoop[3];
    for (int i=0; i<3; i++) {
        oneLoop[i] = - 5.0/9.0 + log_l[i]/3.0 - 2.0/xl[i];
        twoLoop[i] = - 5.0/24.0 + zeta3 + log_l[i]/4.0 + 3.0/xl[i]*log_l[i];
        threeLoop[i] = - 121.0/48.0 + (-5.0 + 8.0*log2)*zeta2 - 99.0/16.0*zeta3
                       + 10.0*zeta5 + log_l[i]/8.0;
        for (int j=0; j<3; j++) {
            if (i>j) { /* Pi^{(2)}_l */
                threeLoop[i] += - 116.0/27.0 + 4.0/3.0*zeta2 + 38.0/9.0*zeta3
                                + 14.0/9.0*log_l[i]
                                + (5.0/18.0 - 4.0/3.0*zeta3)*log_l[j]
                                + log_l[i]*log_l[i]/6.0
                                - log_l[i]*log_l[j]/3.0;
            } else if (i==j) { /* Pi^{(2)}_F */
                threeLoop[i] += - 307.0/216.0 - 8.0/3.0*zeta2 + 545.0/144.0*zeta3
                                + (11.0/6.0 - 4.0/3.0*zeta3)*log_l[i]
                                - log_l[i]*log_l[i]/6.0;
            } else { /* Pi^{(2)}_h */
                threeLoop[i] += - 37.0/6.0 + 38.0/9.0*zeta3
                                + (11.0/6.0 - 4.0/3.0*zeta3)*log_l[j]
                                - log_l[j]*log_l[j]/6.0;
            }
        }
        threeLoop[i] /= -4.0;
    }

    /* TESTS */
    //for (int i=0; i<3; i++) {
    //    std::cout << ale/M_PI*oneLoop[i] << "  "
    //              << ale/M_PI*ale/M_PI*twoLoop[i] << "  "
    //              << ale/M_PI*ale/M_PI*ale/M_PI*threeLoop[i] << std::endl;
    //}

    return ( ale/M_PI*(oneLoop[0] + oneLoop[1] + oneLoop[2])
             + ale/M_PI*ale/M_PI*(twoLoop[0] + twoLoop[1] + twoLoop[2])
             + ale/M_PI*ale/M_PI*ale/M_PI
               *(threeLoop[0] + threeLoop[1] + threeLoop[2]) );
}

double StandardModel::dAleTopMz() const {
    /*
     * Eq.(12) in J.H.Kuhn, M.Steinhauser, PLB437,425(1998) [hep-ph/9802241]
     *
     * (-0.70+-0.05)*10^{-4} for mt=175.6+-5.5 and alpha_s(mZ)=0.118
     */
    double xt = mZ*mZ/getMass(TOP)/getMass(TOP);
    double log_t = log(xt);

    return ( -4.0/45.0*ale/M_PI*xt
             * (1.0 + 5.062*alsMz/M_PI
                + (28.220 + 9.702*log_t)*alsMz/M_PI*alsMz/M_PI
                + xt * (0.1071 + 0.8315*alsMz/M_PI
                        + (6.924 + 1.594*log_t)*alsMz/M_PI*alsMz/M_PI)) );
}

double StandardModel::dAleTotalMz() const {
    return ( dAle5Mz + dAleLepMz() + dAleTopMz() );
}

double StandardModel::aleMz() const {
    return ( ale/(1.0 - dAleTotalMz()) );
}

double StandardModel::mcMz() const {
    double mc_at_mb = mrun(getMass(BOTTOM), getMass(CHARM), 4.0);
    double mc_at_mZ =  mrun(mZ, mc_at_mb, 5.0);

    /* TEST */
    //std::cout << "mc(mc)= " << getMass(CHARM) << std::endl;
    //std::cout << "mc(mb)_LO+NLO= " << mc_at_mb << std::endl;
    //std::cout << "mc(mZ)_LO+NLO= " << mc_at_mZ << std::endl;

    return ( mc_at_mZ );
   //return ( 0.563817 );// <--- used in ZFITTER with the effective mass mc=1.5
}

double StandardModel::mbMz() const {
    /* TEST */
    //std::cout << "mb(mb)= " << getMass(BOTTOM) << std::endl;
    //std::cout << "mb(mZ)_LO= " << mrun(mZ, getMass(BOTTOM), 5.0, 0) << std::endl;
    //std::cout << "mb(mZ)_LO+NLO= " << mrun(mZ, getMass(BOTTOM), 5.0, 1) << std::endl;

    return ( mrun(mZ, getMass(BOTTOM), 5.0) );
    //return (2.819440);// <--- used in ZFITTER with the effective mass mb=4.7
}

double StandardModel::mW() const {

    std::cout << "Write codes for StandardModel::mW() " << std::endl;
    return (80.3613);
}

gslpp::complex StandardModel::gZf(const int INDF) const {

    std::cout << "Write codes for StandardModel::gZf() " << std::endl;
    gslpp::complex tmp(0.0738065, -0.0120949, false);
    return (tmp);
}

gslpp::complex StandardModel::rhoZf(const int INDF) const {

    std::cout << "Write codes for StandardModel::rhoZf() " << std::endl;
    gslpp::complex tmp(1.00516, -0.00473674, false);
    return (tmp);
}

double StandardModel::Delta_r() const {

    std::cout << "Write codes for StandardModel::Delta_r() " << std::endl;
    return (0.0378211);
}


////////////////////////////////////////////////////////////////////////

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

  xt = pow(mrun(muw,getMass(TOP),5.,LE)/mW(),2.);

  eta2b = pow(als(muw,LE)/als(getMass(BOTTOM),LE),
          g0p/2./beta0(5.))*(1.+(LE == 1 ? 1. : 0.)/4./M_PI*
          (als(muw,LE)*(-jp5+kt(xt,muw)/S(xt,xt))+
          als(getMass(BOTTOM),LE)*jp5));

//  printf("%e %e\n",eta2b*pow(als(Mb_v,le),-g0p/2./beta0(5.))*(1.+le/4./PI*als(Mb_v,le)*jp5),als(Mb_v,le));
  return(eta2b);
}

gslpp::complex StandardModel::getDBD2Amplitude(const int LE) const {
    double xtW = mrun(muw,getMass(TOP),5.,1)/mW();
    xtW*=xtW;
    return((GF/4./M_PI*mW())*(GF/4./M_PI*mW())*eta2bbar(LE)*getlamt_d()*
          getlamt_d()*S(xtW,xtW)*
          particles[B_D].getDecayconst()*
          particles[B_D].getDecayconst()*particles[B_D].getBpars()[0]*4./3.*
          particles[B_D].getMass());
}
