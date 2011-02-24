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

//we still jhave to decide where to put them.... constants
const double sw2=0.23;
const double cw2=1.0-sw2;
const double Qf[3]={2.0/3.0,-1.0/3.0,-1.0};
//////////////////////////

class EWPOSM : public StandardModel {
    /**
     * @brief copy constructor
     * @param orig reference to a StandardModel object
     */
    EWPOSM(const Parameters& Par);

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
gsl_complex fcurve_c(double Q2, double M12, double M22);
gsl_complex fcurve2(double s,double Mv2);
gsl_complex fcurve3(double s,double Mv2);
gsl_complex I0(double Q2, double M12, double M22);
gsl_complex I1(double Q2, double M12, double M22);
gsl_complex I3(double Q2, double M12, double M22);
gsl_complex LL(double Q2,double M12, double M22);


double deltarrem (double Mw2);
gsl_complex deltarhowF();
gsl_complex v1W(double s);
gsl_complex v1Z(double s);
gsl_complex v2W(double s);
gsl_complex piZg();
double sigmaBF0();
gsl_complex sigmaBFMw();
gsl_complex sigmaBFMz();
gsl_complex sigmaBFMz1();
gsl_complex sigmaFFMw();
gsl_complex sigmaFFMz();
gsl_complex sigmaFFMz1();
gsl_complex uf(const std::string& ferm);
double vf(const std::string& ferm);
gsl_complex deltarrem();
gsl_complex deltakrem(const std::string& ferm);
gsl_complex deltarhorem(const std::string& ferm);
gsl_complex deltarhozF();

/**
 *
 * @param ferm can be equal to "u", "d", "l"
 * @return fromfactor rho from eq(59) of GFitter
 */
gsl_complex rhoZf(const std::string& ferm);

/**
 *
 * @param ferm can be equal to "u", "u", "l"
 * @return fromfactor k from eq(59) of GFitter
 */

gsl_complex kZf(const std::string& ferm);

double F1(double x);
double F1inf(double x);
double deltarremQCD();
gsl_complex deltarremtot();
/**
 *
 * @param lept can be equal to "e", "mu", "tau"
 * @return partial decay width
 */
double Gamma_f(const std::string& lept);

};



#endif	/* EWPOSM_H */

