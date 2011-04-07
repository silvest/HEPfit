/* 
 * File:   EWPOSM.h
 * Author: aleksandr
 *
 * Created on February 23, 2011, 5:48 PM
 */

#ifndef EWPOSM_H
#define	EWPOSM_H
#include <gslpp_complex.h>
#include <gslpp_vector_double.h>
#include <gslpp_vector_complex.h>
#include <gslpp_matrix_double.h>
#include <gslpp_matrix_complex.h>
#include <gsl/gsl_sf.h>
#include "QCD.h"
#include "StandardModel.h"
#include <math.h>
#define LEPTON   1.0

//we still jhave to decide where to put them.... constants
const double sw2=0.23;
const double cw2=1.0-sw2;
const double Qf[3]={2.0/3.0,-1.0/3.0,-1.0};
const double GammaW0=2.085;
const double GammaZ0=2.4952;
//////////////////////////

class EWPOSM {
public:
    /**
     * @brief copy constructor
     * @param orig reference to a StandardModel object
     */
    EWPOSM(StandardModel& sm);
   
    /**
     * @brief EWPOSM destructor
     */
   virtual ~EWPOSM();



  double mW() const;

         /**
     *
     * @param ferm defines the type of the fermions: leptons -"l", charm -"c", or b quark "b"
     * @return sin theta effective for the given fermion ferm @f$\sin^2\theta_{\mathrm{eff}}^\ell$@f
     */
    double sin2thwall(const std::string& ferm) const;

/**
 * 
 * @return different functions and integrals from Bardin et al  Z.Phys.C 44,493-502 (1989)
 */
double lambda(double Q2, double M12, double M22 );
double lambdas(double s, double Mv2);
double fcurve(double Q2, double M12, double M22);
gslpp::complex fcurve_c(double Q2, double M12, double M22);
gslpp::complex fcurve2(double s,double Mv2);
gslpp::complex fcurve3(double s,double Mv2);
gslpp::complex I0(double Q2, double M12, double M22);
gslpp::complex I1(double Q2, double M12, double M22);
gslpp::complex I3(double Q2, double M12, double M22);
gslpp::complex LL(double Q2,double M12, double M22);
gslpp::complex DI0(double Q2,double AMZ2, double M12, double M22);
gslpp::complex DI3(double Q2, double AMZ2,double M12, double M22);
gslpp::complex DL(double Q2, double AMZ2,double M12, double M22);


gslpp::complex deltarrem ();
gslpp::complex deltarhowF();
gslpp::complex v1W(double s);
gslpp::complex v1Z(double s);
gslpp::complex v2W(double s);
gslpp::complex v2Wa();
gslpp::complex piZg();
double sigmaBF0();
gslpp::complex sigmaBFMw();
gslpp::complex sigmaBFMz();
gslpp::complex sigmaBFMz1();
gslpp::complex sigmaFww(double s);
gslpp::complex sigmaFzz();
gslpp::complex sigmaFzz1();
gslpp::complex sigmaFzz1a();
gslpp::complex uf(const std::string& ferm);
double vf(const std::string& ferm);
gslpp::complex deltakrem(const std::string& ferm);
gslpp::complex deltarhorem(const std::string& ferm);
gslpp::complex deltarhozF();

/**
 *
 * @param ferm can be equal to "u", "d", "l"
 * @return fromfactor rho from eq(59) of GFitter
 */
gslpp::complex rhoZf(const std::string& ferm);

/**
 *
 * @param ferm can be equal to "u", "u", "l"
 * @return fromfactor k from eq(59) of GFitter
 */

gslpp::complex kZf(const std::string& ferm);

double F1(double x);
double F1inf(double x);
double deltarremQCD();
double deltarhoQCD();
gslpp::complex deltarremtot();
/**
 *
 * @param lept can be equal to "e", "mu", "tau"
 * @return partial decay width
 */
double Gamma_f(const std::string& lept);

double Gamma_Z() const;
double Gamma_W() const;
double sw2() const;
double cw2() const;
double rekZf(const std::string& ferm )const;

gslpp::complex deltarho();
gslpp::complex Deltarho();

/**
 *
 * @return delta r from formula (16) of ZFitter
 */
double deltar() const;

gslpp::complex gzf(const std::string& ferm);
gslpp::complex gzf1(const std::string& ferm);
gslpp::complex spence(double);
gslpp::complex fspence(double);
gslpp::complex  XAMM1();
double SCALE();
double ROQCD();
double AKQCD();
double CLQQCD();


protected:
    double me, mmu, mtau, mnu1, mnu2, mnu3,mu,md,ms,mc,mt,mb;
    double mHl, alsMz, ale, mZ, GF, dAletotmz;

};



#endif	/* EWPOSM_H */

