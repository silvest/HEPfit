/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MVLL_H
#define	MVLL_H

#include <math.h>
#include <StandardModel.h>
#include <ThObservable.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <assert.h>


#define CUTOFF 10    //cutoff between LCSR and lattice values for Form Factors, in GeV^2

/*******************************************************************************
 * GSL Function Conversion BEGIN                                                  *
 * ****************************************************************************/

template<class F>
static double gslFunctionAdapter( double x, void* p)
{
    // Here I do recover the "right" pointer, safer to use static_cast
    // than reinterpret_cast.
    F* function = static_cast<F*>( p );
    return (*function)( x );
}

template<class F>
gsl_function convertToGslFunction( const F& f )
{
    gsl_function gslFunction;
    
    const void* p = &f;
    assert (p != 0);
    
    gslFunction.function = &gslFunctionAdapter<F>;
    // Just to eliminate the const.
    gslFunction.params = const_cast<void*>( p );
    
    return gslFunction;
}

/*******************************************************************************
 * GSL Function conversion END                                                     *
 * ****************************************************************************/

/**
 * @class MVll
 * @ingroup Flavour
 * @brief A class for the @f$M \to V l^+ l^-@f$ decay.  
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to 
 * compute the observables relative to the @f$M \to V l^+ l^-@f$ decay. After the
 * parameters are updated in updateParameters() and the cache is checked in 
 * checkCache(), the form factor are build in the transverse basis in the functions
 * V(), A_0(), A_1(), A_2(), T_1(), T_2(), T_3tilde() and T_3() using either the LCSR functions
 * LCSR_fit1(), LCSR_fit2() and LCSR_fit3() or the lattice function lat_fit().
 * The form factor are consequentely translated in the helicity basis through the
 * functions V_L(), V_R(), T_L(), T_R(), S_L() and S_R(). Form factors and parameters
 * are combined together in the functions H_V(), H_A(), H_S() and H_P() in order
 * to build the helicity aplitudes, which are consequentely combined to create
 * the angular coefficients in the function I(). Those coefficients are used to
 * create the CP averaged coefficients in the function Sigma() ad the CP asymmetric
 * coefficients in the function Delta(). Form factors, CP averaged and asymmetric
 * coefficients and hadronic contributions are integrated in the functions 
 * integrateSigma(), integrateDelta(), integrategtilde() and integrateh() in order
 * to be further used to build the observables.
 */
class MVll {
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     * @param[in] lep_i final leptons of the decay
     */
    MVll(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~MVll();
    
    /**
     * @brief The update parameter method for MVll.
     */
    void updateParameters();
    
    /**
     * @brief The caching method for MVll.
     */
    void checkCache();
    
    double GF;            /**<Fermi constant */
    double ale;           /**<alpha electromagnetic */
    double Mlep;          /**<muon mass */
    double MM;            /**<initial meson mass */
    double MV;            /**<final vector meson mass */
    double Mb;            /**<b quark mass */
    double mu_b;          /**<b mass scale */
    double Mc;            /**<c quark mass */
    double Ms;            /**<s quark mass */
    double width;         /**<initial meson width */
    double MW;            /**<W boson mass */
    gslpp::complex lambda_t;     /**<Vckm factor */
    double b;             /**<BF of the decay V -> final states */
    gslpp::complex h_0[3];         /**<parameter that contains the contribution from the hadronic hamiltonian */
    gslpp::complex h_1[3];         /**<parameter that contains the contribution from the hadronic hamiltonian */
    gslpp::complex h_2[3];         /**<parameter that contains the contribution from the hadronic hamiltonian */
    double q2;            /**<\f$q^2\f$ of the decay */
    
    /*lattice fit parameters*/
    double a_0V;/**<lattice fit parameter */
    double a_1V;/**<lattice fit parameter */
    double a_2V;/**<lattice fit parameter */
    double dmV;/**<lattice fit parameter */
    double a_0A0;/**<lattice fit parameter */
    double a_1A0;/**<lattice fit parameter */
    double a_2A0;/**<lattice fit parameter */
    double dmA0;/**<lattice fit parameter */
    double a_0A1;/**<lattice fit parameter */
    double a_1A1;/**<lattice fit parameter */
    double a_2A1;/**<lattice fit parameter */
    double dmA1;/**<lattice fit parameter */
    double a_0A12;/**<lattice fit parameter */
    double a_1A12;/**<lattice fit parameter */
    double a_2A12;/**<lattice fit parameter */
    double dmA12;/**<lattice fit parameter */
    double a_0T1;/**<lattice fit parameter */
    double a_1T1;/**<lattice fit parameter */
    double a_2T1;/**<lattice fit parameter */
    double dmT1;/**<lattice fit parameter */
    double a_0T2;/**<lattice fit parameter */
    double a_1T2;/**<lattice fit parameter */
    double a_2T2;/**<lattice fit parameter */
    double dmT2;/**<lattice fit parameter */
    double a_0T23;/**<lattice fit parameter */
    double a_1T23;/**<lattice fit parameter */
    double a_2T23;/**<lattice fit parameter */
    double dmT23;/**<lattice fit parameter */

    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients */
    gslpp::vector<gslpp::complex> ** allcoeffprime;/**<vector that contains the primed Wilson coeffients */
    
    gslpp::complex C_1;/**<Wilson coeffients @f$C_1@f$*/
    gslpp::complex C_2;/**<Wilson coeffients @f$C_2@f$*/
    gslpp::complex C_3;/**<Wilson coeffients @f$C_3@f$*/
    gslpp::complex C_4;/**<Wilson coeffients @f$C_4@f$*/
    gslpp::complex C_5;/**<Wilson coeffients @f$C_5@f$*/
    gslpp::complex C_6;/**<Wilson coeffients @f$C_6@f$*/
    gslpp::complex C_7;/**<Wilson coeffients @f$C_7@f$*/
    gslpp::complex C_9;/**<Wilson coeffients @f$C_9@f$*/
    gslpp::complex C_10;/**<Wilson coeffients @f$C_{10}@f$*/
    gslpp::complex C_S;/**<Wilson coeffients @f$C_S@f$*/
    gslpp::complex C_P;/**<Wilson coeffients @f$C_P@f$*/
    
    gslpp::complex C_7p;/**<Wilson coeffients @f$C_7'@f$*/
    gslpp::complex C_9p;/**<Wilson coeffients @f$C_9'@f$*/
    gslpp::complex C_10p;/**<Wilson coeffients @f$C_{10}'@f$*/
    gslpp::complex C_Sp;/**<Wilson coeffients @f$C_S'@f$*/
    gslpp::complex C_Pp;/**<Wilson coeffients @f$C_P'@f$*/
    
    
    /**
    * @brief The fit function from arXiv:1503.05534v1, \f$ f^{LCSR} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] a_0 fit parameter
    * @param[in] a_1 fit parameter
    * @param[in] a_2 fit parameter
    * @param[in] dm shift in the initial meson mass
    * @return \f$ f^{lat} \f$
    */
    double LCSR_fit(double q2, double a_0, double a_1, double a_2, double dm);
    
    
    /**
    * @brief The lattice parameter \f$ z \f$ from arXiv:1310.3722v3.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ z \f$
    */
    double z(double q2);
    
    
    /**
    * @brief The fit function from arXiv:1310.3722v3, \f$ f^{lat} \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] a_0 fit parameter
    * @param[in] a_1 fit parameter
    * @param[in] c_01 fit parameter
    * @param[in] c_01s fit parameter
    * @param[in] dm shift in the initial meson mass
    * @return \f$ f^{lat} \f$
    */
    double lat_fit(double q2, double a_0, double a_1, double dm);
    
    
    /**
    * @brief The transverse form factor \f$ V \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V \f$
    */
    double V(double q2);

    
    /**
    * @brief The transverse form factor \f$ A_0 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_0 \f$
    */
    double A_0(double q2);

    
    /**
    * @brief The transverse form factor \f$ A_1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_1 \f$
    */
    double A_1(double q2);

    
    /**
    * @brief The transverse form factor \f$ A_2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ A_2 \f$
    */
    double A_2(double q2);

    
    /**
    * @brief The transverse form factor \f$ T_1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_1 \f$
    */
    double T_1(double q2);

    
    /**
    * @brief The transverse form factor \f$ T_2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_2 \f$
    */
    double T_2(double q2);
    
    /**
    * @brief The helicity form factor \f$ V_L^{\lambda} \f$.
    * @param[in] i polarization: 0 for 0, 1 for +, 2 for -
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_L^{\lambda} \f$
    */
    double V_L(int i, double q2);

    
    /**
    * @brief The helicity form factor \f$ V_R^{\lambda} \f$.
    * @param[in] i polarization: 0 for 0, 1 for +, 2 for -
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_R^{\lambda} \f$
    */
    double V_R(int i, double q2);


    /**
    * @brief The helicity form factor \f$ T_L^{\lambda} \f$.
    * @param[in] i polarization: 0 for 0, 1 for +, 2 for -
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_L^{\lambda} \f$
    */
    double T_L(int i, double q2);


    /**
    * @brief The helicity form factor \f$ T_R^{\lambda} \f$.
    * @param[in] i polarization: 0 for 0, 1 for +, 2 for -
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_R^{\lambda} \f$
    */
    double T_R(int i, double q2);


    /**
    * @brief The helicity form factor \f$ S_L \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ S_L \f$
    */
    double S_L(double q2);


    /**
    * @brief The helicity form factor \f$ S_R \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ S_R \f$
    */
    double S_R(double q2);
    
    
    /**
    * @brief The helicity amplitudes normalization factor \f$ N \f$ .
    * @return \f$ N \f$
    */
    gslpp::complex N();
    
    
    /**
    * @brief The \f$ h(q^2,m) \f$ function involved into \f$ C_9^{eff}\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] m mass
    * @return \f$ h(q^2,m) \f$
    */
    gslpp::complex H(double q2, double m);
    
    
    /**
    * @brief The \f$ Y(q^2) \f$ function involved into \f$ C_9^{eff}\f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ Y(q^2) \f$
    */
    gslpp::complex Y(double q2);
    
    
    /**
    * @brief The helicity amplitude \f$ H_V^{\lambda} \f$ .
    * @param[in] i polarization: 0 for 0, 1 for +, 2 for -
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ H_V^{\lambda} \f$
    */
    gslpp::complex H_V(int i, double q2, int bar);


    /**
    * @brief The helicity amplitude \f$ H_A^{\lambda} \f$ .
    * @param[in] i polarization: 0 for 0, 1 for +, 2 for -
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ H_A^{\lambda} \f$
    */
    gslpp::complex H_A(int i, double q2, int bar);


    /**
    * @brief The helicity amplitude \f$ H_S^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ H_S^{\lambda} \f$
    */
    gslpp::complex H_S(double q2, int bar);


    /**
    * @brief The helicity amplitude \f$ H_P^{\lambda} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ H_P^{\lambda} \f$
    */
    gslpp::complex H_P(double q2, int bar);
    
    
    /**
    * @brief The square of the 3-momentum of the recoiling meson in the M rest frame, \f$ k^2 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ k^2 \f$ 
    */
    double k2 (double q2);
    
    
    /**
    * @brief The factor \f$ \beta \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \beta \f$
    */
    double beta (double q2);
    
    
    /**
    * @brief The factor \f$ \lambda \f$ used in the angular coefficients \f$I_i\f$. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \lambda \f$
    */
    double lambda(double q2);

    
    /**
    * @brief The factor \f$ F \f$ used in the angular coefficients I_i. 
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] b_i BF of the decay \f$ V \to M_1 M_2\f$ 
    * @return \f$ F \f$
    */
    double F(double q2, double b_i);
    
    
    /**
    * @brief The angular coefficient \f$ I_{i} \f$ .
    * @param[in] i index of the angular coefficient: 0 for 1c, 1 for 1s, 2 for 2c,
    *  3 for 2s, 4 for 3, 5 for 4, 6 for 5, 7 for 6s, 8 for 6c, 9 for 7, 10 for 8, 11 for 9
    * @param[in] q2 \f$q^2\f$ of the decay
    * @param[in] bar index to choose betwen regular coefficient (bar=0) and conjugated coefficient (bar=1)
    * @return \f$ I_{i} \f$
    */
    double  I(int i, double q2, int bar);
    
    
    /**
    * @brief The CP average \f$ \Sigma_{i} \f$ .
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{i} \f$
    */
    double Sigma(int i, double q2);
    
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{i} \f$ .
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{i} \f$
    */
    double Delta(int i, double q2);
    
    /**
    * @brief The integral of \f$ \Sigma_{i} \f$ from \f$q_{min}\f$ to \f$q_{max}\f$
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return \f$ <\Sigma_{i}> \f$ 
    */
    double integrateSigma(int i, double q_min, double q_max);
    
    /**
    * @brief The integral of \f$ \Delta_{i} \f$ from \f$q_{min}\f$ to \f$q_{max}\f$
    * @param[in] i index of the angular coefficient \f$ I_{i} \f$
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return \f$ <\Delta_{i}> \f$ 
    */
    double integrateDelta(int i, double q_min, double q_max);
    
    /**
    * @brief The integral of the hadronic correction \f$ \tilde{g}^i \f$ from \f$q_{min}\f$ to \f$q_{max}\f$
    * @param[in] i index of the hadronic correction \f$ \tilde{g}^i \f$
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return \f$ <\tilde{g}^i> \f$ 
    */
    gslpp::complex integrategtilde(int i, double q_min, double q_max);
    
    /**
    * @brief The integral of the hadronic correction \f$ h_{\lambda} \f$ from \f$q_{min}\f$ to \f$q_{max}\f$
    * @param[in] i index of the hadronic correction \f$ h_{\lambda} \f$
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return \f$ <h_{\lambda}> \f$ 
    */
    gslpp::complex integrateh(int i, double q_min, double q_max);
    
    /**
    * @brief The integral of a form factor from \f$q_{min}\f$ to \f$q_{max}\f$.
    * @param[in] i index of the form factor: 0 for \f$ V_0 \f$, 1 for \f$ V_+ \f$,
     *  2 for \f$ V_- \f$, 3 for \f$ T_0 \f$, 4 for \f$ T_+ \f$, 5 for \f$ T_- \f$,
     * 6 for \f$ S \f$
    * @param[in] q_min minimum q^2 of the integral
    * @param[in] q_max maximum q^2 of the integral
    * @return \f$ <FF> \f$ 
    */
    double integrateFF(int i, double q_min, double q_max);
    
    /**
    * @brief The CP average \f$ \Sigma_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{1s} \f$
    */
    double getSigma0(double q2)
    {
        return Sigma(0, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{1c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{1c} \f$
    */
    double getSigma1(double q2)
    {
        return Sigma(1, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{2s} \f$
    */
    double getSigma2(double q2)
    {
        return Sigma(2, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{2c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{2c} \f$
    */
    double getSigma3(double q2)
    {
        return Sigma(3, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{3} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{3} \f$
    */ 
    double getSigma4(double q2)
    {
        return Sigma(4, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{4} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{4} \f$
    */ 
    double getSigma5(double q2)
    {
        return Sigma(5, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{5} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{5} \f$
    */ 
    double getSigma6(double q2)
    {
        return Sigma(6, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{6s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{6s} \f$
    */ 
    double getSigma7(double q2)
    {
        return Sigma(7, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{7} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{7} \f$
    */ 
    double getSigma9(double q2)
    {
        return Sigma(9, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{8} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{8} \f$
    */ 
    double getSigma10(double q2)
    {
        return Sigma(10, q2);
    };
    
    /**
    * @brief The CP average \f$ \Sigma_{9} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Sigma_{9} \f$
    */ 
    double getSigma11(double q2)
    {
        return Sigma(11, q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{1s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{1s} \f$
    */
    double getDelta0(double q2)
    {
        return Delta(0, q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{1c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{1c} \f$
    */
    double getDelta1(double q2)
    {
        return Delta(1, q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{2s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{2s} \f$
    */
    double getDelta2(double q2)
    {
        return Delta(2, q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{2c} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{2c} \f$
    */
    double getDelta3(double q2)
    {
        return Delta(3, q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{6s} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{6s} \f$
    */
    double getDelta7(double q2)
    {
        return Delta(7, q2);
    };
    
    /**
    * @brief The CP asymmetry \f$ \Delta_{9} \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \Delta_{9} \f$
    */
    double getDelta11(double q2)
    {
        return Delta(11, q2);
    };
    
    /**
    * @brief The square of the absolute value of the helicity amplitude \f$ H_V^0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ |H_V^0|^2 \f$
    */
    double getHV0_abs2(double q2)
    {
        return H_V(0, q2, 0).abs2();
    };
    
    /**
    * @brief The square of the absolute value of the helicity amplitude \f$ H_V^+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ |H_V^+|^2 \f$
    */
    double getHV1_abs2(double q2)
    {
        return H_V(1, q2, 0).abs2();
    };
    
    /**
    * @brief The square of the absolute value of the helicity amplitude \f$ H_V^- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ |H_V^-|^2 \f$
    */
    double getHV2_abs2(double q2)
    {
        return H_V(2, q2, 0).abs2();
    };
    
    /**
    * @brief The square of the absolute value of the helicity amplitude \f$ H_A^0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ |H_A^0|^2 \f$
    */
    double getHA0_abs2(double q2)
    {
        return H_A(0, q2, 0).abs2();
    };
    
    /**
    * @brief The square of the absolute value of the helicity amplitude \f$ H_A^+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ |H_A^+|^2 \f$
    */
    double getHA1_abs2(double q2)
    {
        return H_A(1, q2, 0).abs2();
    };
    
    /**
    * @brief The square of the absolute value of the helicity amplitude \f$ H_A^- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ |H_A^-|^2 \f$
    */
    double getHA2_abs2(double q2)
    {
        return H_A(2, q2, 0).abs2();
    };
    
    /**
    * @brief The square of the absolute value of the helicity amplitude \f$ H_S \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ |H_S|^2 \f$
    */
    double getHS_abs2(double q2)
    {
        return H_S(q2, 0).abs2();
    }
    
    /**
    * @brief The square of the absolute value of the helicity amplitude \f$ H_P \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ |H_P|^2 \f$
    */
    double getHP_abs2(double q2)
    {
        return H_P(q2, 0).abs2();
    }
    
    /**
    * @brief The form factor \f$ V_0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_0 \f$
    */
    double getV0(double q2)
    {
        return (2. * MM * sqrt(q2))/sqrt(lambda(q2)) * V_L(0,q2);
    };
    
    /**
    * @brief The form factor \f$ V_+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_+ \f$
    */
    double getVp(double q2)
    {
        return V_L(1,q2);
    };
    
    /**
    * @brief The form factor \f$ V_- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ V_- \f$
    */
    double getVm(double q2)
    {
        return V_L(2,q2);
    };
    
    /**
    * @brief The form factor \f$ T_0 \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_0 \f$
    */
    double getT0(double q2)
    {
        return (2. * pow(MM, 3.))/sqrt(q2 * lambda(q2)) * T_L(0,q2);
    };
    
    /**
    * @brief The form factor \f$ T_+ \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_+ \f$
    */
    double getTp(double q2)
    {
        return T_L(1,q2);
    };
    
    /**
    * @brief The form factor \f$ T_- \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ T_- \f$
    */
    double getTm(double q2)
    {
        return T_L(2,q2);
    };
    
    /**
    * @brief The form factor \f$ S \f$ .
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ S \f$
    */
    double getS(double q2)
    {
        return (-2. * MM * (Mb + Ms))/sqrt(lambda(q2)) * S_L(q2);
    };
    
    /**
    * @brief The real part of \f$ \tilde{g}^1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(\tilde{g}^1) \f$ 
    */
    double getgtilde_1_re(double q2)
    {
        updateParameters();
        return 1./(2. * C_2.real()) * (-16.*pow(MM,3.)*(MM + MV)*pow(M_PI,2.)/(sqrt(lambda(q2)) * V(q2)) * (h_0[2]/q2 + h_1[2] + h_2[2] * q2 - (h_0[1]/q2 + h_1[1] + h_2[1] * q2))).real();
    }
    
    /**
    * @brief The immaginary part of \f$ \tilde{g}^1 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(\tilde{g}^1) \f$ 
    */
    double getgtilde_1_im(double q2)
    {
        updateParameters();
        return 1./(2. * C_2.real()) * (-16.*pow(MM,3.)*(MM + MV)*pow(M_PI,2.)/(sqrt(lambda(q2)) * V(q2)) * (h_0[2]/q2 + h_1[2] + h_2[2] * q2 - (h_0[1]/q2 + h_1[1] + h_2[1] * q2))).imag();
    }
    
    /**
    * @brief The real part of \f$ \tilde{g}^2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(\tilde{g}^2) \f$ 
    */
    double getgtilde_2_re(double q2)
    {
        updateParameters();
        return 1./(2. * C_2.real()) * (-16.*pow(MM,3.)*pow(M_PI,2.)/((MM + MV) * A_1(q2)) * (h_0[2]/q2 + h_1[2] + h_2[2] * q2 + h_0[1]/q2 + h_1[1] + h_2[1] * q2)).real();
    }
    
    /**
    * @brief The immaginary part of \f$ \tilde{g}^2 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(\tilde{g}^2) \f$ 
    */
    double getgtilde_2_im(double q2)
    {
        updateParameters();
        return 1./(2. * C_2.real()) * (-16.*pow(MM,3.)*pow(M_PI,2.)/((MM + MV) * A_1(q2)) * (h_0[2]/q2 + h_1[2] + h_2[2] * q2 + h_0[1]/q2 + h_1[1] + h_2[1] * q2)).imag();
    }
    
    /**
    * @brief The real part of \f$ \tilde{g}^3 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(\tilde{g}^3) \f$ 
    */
    double getgtilde_3_re(double q2)
    {
        updateParameters();
        return 1./(2. * C_2.real()) * (64.*pow(MM,3.)*pow(M_PI,2.)*MV*(MM + MV)/(lambda(q2) * A_2(q2)) * (sqrt(q2)*(h_0[0]/q2 + h_1[0] + h_2[0] * q2)-(MM*MM - q2 - MV*MV)/(4.*MV) * (h_0[2]/q2 + h_1[2] + h_2[2] * q2 + h_0[1]/q2 + h_1[1] + h_2[1] * q2))).real();
    }

    /**
    * @brief The immaginary part of \f$ \tilde{g}^3 \f$.
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(\tilde{g}^3) \f$ 
    */
    double getgtilde_3_im(double q2)
    {
        updateParameters();
        return 1./(2. * C_2.real()) * (64.*pow(MM,3.)*pow(M_PI,2.)*MV*(MM + MV)/(lambda(q2) * A_2(q2)) * (sqrt(q2)*(h_0[0]/q2 + h_1[0] + h_2[0] * q2)-(MM*MM - q2 - MV*MV)/(4.*MV) * (h_0[2]/q2 + h_1[2] + h_2[2] * q2 + h_0[1]/q2 + h_1[1] + h_2[1] * q2))).imag();
    }
    
    /**
    * @brief The real part of \f$ h_0 \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(h_0) \f$
    */
    double geth_0_re(double q2)
    {
        return (16 * M_PI * M_PI * MM * MM * (h_0[0]/q2 + h_1[0] + h_2[0] * q2)).real();
        //return (-8.*M_PI*M_PI*pow(MM, 3.)/(sqrt(lambda(q2))*T_1(q2)*(Mb+Ms))*(h_0[2]-h_0[1])).real();
    }
    
    /**
    * @brief The immaginary part of \f$ h_0 \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(h_0) \f$
    */
    double geth_0_im(double q2)
    {
        return (16 * M_PI * M_PI * MM * MM * (h_0[0]/q2 + h_1[0] + h_2[0] * q2)).imag();
        //return (-8.*M_PI*M_PI*pow(MM, 3.)/(sqrt(lambda(q2))*T_1(q2)*(Mb+Ms))*(h_0[2]-h_0[1])).imag();
    }
    
    /**
    * @brief The real part of \f$ h_+ \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(h_+) \f$
    */
    double geth_p_re(double q2)
    {
        return (16 * M_PI * M_PI * MM * MM * (h_0[1]/q2 + h_1[1] + h_2[1] * q2)).real();
        //return (-8.*M_PI*M_PI*pow(MM, 3.)/(T_2(q2)*(Mb+Ms)*(MM*MM-MV*MV))*(h_0[2]+h_0[1])).real();
    }
    
    /**
    * @brief The immaginary part of \f$ h_+ \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(h_+) \f$
    */
    double geth_p_im(double q2)
    {
        return (16 * M_PI * M_PI * MM * MM * (h_0[1]/q2 + h_1[1] + h_2[1] * q2)).imag();
        //return (-8.*M_PI*M_PI*pow(MM, 3.)/(T_2(q2)*(Mb+Ms)*(MM*MM-MV*MV))*(h_0[2]+h_0[1])).imag();
    }
    
    /**
    * @brief The real part of \f$ h_- \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{RE}(h_-) \f$
    */
    double geth_m_re(double q2)
    {
        return (16 * M_PI * M_PI * MM * MM * (h_0[2]/q2 + h_1[2] + h_2[2] * q2)).real();
        //return (32.*M_PI*M_PI*pow(MM, 3.)*MV/(lambda(q2)*(Mb-Ms)*(T_2(q2)+q2/(MM*MM-MV*MV)*T_3(q2)))*(h_0[0]*sqrt(q2)-(MM*MM-q2-MV*MV)/(4.*MV)*(h_0[2]+h_0[1]))).real();
    }

    /**
    * @brief The immaginary part of \f$ h_- \f$.  
    * @param[in] q2 \f$q^2\f$ of the decay
    * @return \f$ \mathrm{IM}(h_-) \f$
    */
    double geth_m_im(double q2)
    {
        return (16 * M_PI * M_PI * MM * MM * (h_0[2]/q2 + h_1[2] + h_2[2] * q2)).imag();
        //return (32.*M_PI*M_PI*pow(MM, 3.)*MV/(lambda(q2)*(Mb-Ms)*(T_2(q2)+q2/(MM*MM-MV*MV)*T_3(q2)))*(h_0[0]*sqrt(q2)-(MM*MM-q2-MV*MV)/(4.*MV)*(h_0[2]+h_0[1]))).imag();
    }

private:
    const StandardModel& mySM;/**< Model type */
    StandardModel::lepton lep;/**< Final leptons type */
    StandardModel::meson meson;/**< Initial meson type */
    StandardModel::meson vectorM;/**< Final vector meson type */
    
    std::map<std::pair<double, double>, double > cacheSigma0;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma1;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma2;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma3;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma4;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma5;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma6;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma7;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma9;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma10;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheSigma11;/**< Cache variable */
    
    std::map<std::pair<double, double>, double > cacheDelta0;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta1;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta2;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta3;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta7;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheDelta11;/**< Cache variable */
    
    std::map<std::pair<double, double>, double > cacheVp;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheVm;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheTp;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheTm;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheT0;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheV0;/**< Cache variable */
    std::map<std::pair<double, double>, double > cacheS;/**< Cache variable */
    
    std::map<std::pair<double, double>, gslpp::complex > cachegtilde_1;/**< Cache variable */
    std::map<std::pair<double, double>, gslpp::complex > cachegtilde_2;/**< Cache variable */
    std::map<std::pair<double, double>, gslpp::complex > cachegtilde_3;/**< Cache variable */
    
    std::map<std::pair<double, double>, gslpp::complex > cacheh_0;/**< Cache variable */
    std::map<std::pair<double, double>, gslpp::complex > cacheh_p;/**< Cache variable */
    std::map<std::pair<double, double>, gslpp::complex > cacheh_m;/**< Cache variable */
    
    double avaSigma0;/**< Gsl integral variable */
    double avaSigma1;/**< Gsl integral variable */
    double avaSigma2;/**< Gsl integral variable */
    double avaSigma3;/**< Gsl integral variable */
    double avaSigma4;/**< Gsl integral variable */
    double avaSigma5;/**< Gsl integral variable */
    double avaSigma6;/**< Gsl integral variable */
    double avaSigma7;/**< Gsl integral variable */
    double avaSigma9;/**< Gsl integral variable */
    double avaSigma10;/**< Gsl integral variable */
    double avaSigma11;/**< Gsl integral variable */
    
    double errSigma0;/**< Gsl integral variable */
    double errSigma1;/**< Gsl integral variable */
    double errSigma2;/**< Gsl integral variable */
    double errSigma3;/**< Gsl integral variable */
    double errSigma4;/**< Gsl integral variable */
    double errSigma5;/**< Gsl integral variable */
    double errSigma6;/**< Gsl integral variable */
    double errSigma7;/**< Gsl integral variable */
    double errSigma9;/**< Gsl integral variable */
    double errSigma10;/**< Gsl integral variable */
    double errSigma11;/**< Gsl integral variable */
    
    double avaDelta0;/**< Gsl integral variable */
    double avaDelta1;/**< Gsl integral variable */
    double avaDelta2;/**< Gsl integral variable */
    double avaDelta3;/**< Gsl integral variable */
    double avaDelta7;/**< Gsl integral variable */
    double avaDelta11;/**< Gsl integral variable */
    
    double errDelta0;/**< Gsl integral variable */
    double errDelta1;/**< Gsl integral variable */
    double errDelta2;/**< Gsl integral variable */
    double errDelta3;/**< Gsl integral variable */
    double errDelta7;/**< Gsl integral variable */
    double errDelta11;/**< Gsl integral variable */
    
    double avaVp;/**< Gsl integral variable */
    double avaTp;/**< Gsl integral variable */
    double avaVm;/**< Gsl integral variable */
    double avaTm;/**< Gsl integral variable */
    double avaT0;/**< Gsl integral variable */
    double avaV0;/**< Gsl integral variable */
    double avaS;/**< Gsl integral variable */
    
    double errVp;/**< Gsl integral variable */
    double errTp;/**< Gsl integral variable */
    double errVm;/**< Gsl integral variable */
    double errTm;/**< Gsl integral variable */
    double errT0;/**< Gsl integral variable */
    double errV0;/**< Gsl integral variable */
    double errS;/**< Gsl integral variable */
    
    double avagtilde_1_re;/**< Gsl integral variable */
    double avagtilde_1_im;/**< Gsl integral variable */
    double avagtilde_2_re;/**< Gsl integral variable */
    double avagtilde_2_im;/**< Gsl integral variable */
    double avagtilde_3_re;/**< Gsl integral variable */
    double avagtilde_3_im;/**< Gsl integral variable */
    
    double errgtilde_1_re;/**< Gsl integral variable */
    double errgtilde_1_im;/**< Gsl integral variable */
    double errgtilde_2_re;/**< Gsl integral variable */
    double errgtilde_2_im;/**< Gsl integral variable */
    double errgtilde_3_re;/**< Gsl integral variable */
    double errgtilde_3_im;/**< Gsl integral variable */
    
    double avah_0_re;/**< Gsl integral variable */
    double avah_0_im;/**< Gsl integral variable */
    double avah_p_re;/**< Gsl integral variable */
    double avah_p_im;/**< Gsl integral variable */
    double avah_m_re;/**< Gsl integral variable */
    double avah_m_im;/**< Gsl integral variable */
    
    double errh_0_re;/**< Gsl integral variable */
    double errh_0_im;/**< Gsl integral variable */
    double errh_p_re;/**< Gsl integral variable */
    double errh_p_im;/**< Gsl integral variable */
    double errh_m_re;/**< Gsl integral variable */
    double errh_m_im;/**< Gsl integral variable */
    
    gsl_function FS0;/**< Gsl integral variable */
    gsl_function FS1;/**< Gsl integral variable */
    gsl_function FS2;/**< Gsl integral variable */
    gsl_function FS3;/**< Gsl integral variable */
    gsl_function FS4;/**< Gsl integral variable */
    gsl_function FS5;/**< Gsl integral variable */
    gsl_function FS6;/**< Gsl integral variable */
    gsl_function FS7;/**< Gsl integral variable */
    gsl_function FS9;/**< Gsl integral variable */
    gsl_function FS10;/**< Gsl integral variable */
    gsl_function FS11;/**< Gsl integral variable */
    
    gsl_function FD0;/**< Gsl integral variable */
    gsl_function FD1;/**< Gsl integral variable */
    gsl_function FD2;/**< Gsl integral variable */
    gsl_function FD3;/**< Gsl integral variable */
    gsl_function FD7;/**< Gsl integral variable */
    gsl_function FD11;/**< Gsl integral variable */
    
    gsl_function FVp;/**< Gsl integral variable */
    gsl_function FVm;/**< Gsl integral variable */
    gsl_function FTp;/**< Gsl integral variable */
    gsl_function FTm;/**< Gsl integral variable */
    gsl_function FV0;/**< Gsl integral variable */
    gsl_function FT0;/**< Gsl integral variable */
    gsl_function FS;/**< Gsl integral variable */
    
    gsl_function Fgtilde_1_re;/**< Gsl integral variable */
    gsl_function Fgtilde_1_im;/**< Gsl integral variable */
    gsl_function Fgtilde_2_re;/**< Gsl integral variable */
    gsl_function Fgtilde_2_im;/**< Gsl integral variable */
    gsl_function Fgtilde_3_re;/**< Gsl integral variable */
    gsl_function Fgtilde_3_im;/**< Gsl integral variable */
    
    gsl_function Fh_0_re;/**< Gsl integral variable */
    gsl_function Fh_0_im;/**< Gsl integral variable */
    gsl_function Fh_p_re;/**< Gsl integral variable */
    gsl_function Fh_p_im;/**< Gsl integral variable */
    gsl_function Fh_m_re;/**< Gsl integral variable */
    gsl_function Fh_m_im;/**< Gsl integral variable */
    
    gsl_integration_workspace * w_sigma0;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma1;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma2;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma3;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma4;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma5;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma6;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma7;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma9;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma10;/**< Gsl integral variable */
    gsl_integration_workspace * w_sigma11;/**< Gsl integral variable */
    
    gsl_integration_workspace * w_delta0;/**< Gsl integral variable */
    gsl_integration_workspace * w_delta1;/**< Gsl integral variable */
    gsl_integration_workspace * w_delta2;/**< Gsl integral variable */
    gsl_integration_workspace * w_delta3;/**< Gsl integral variable */
    gsl_integration_workspace * w_delta7;/**< Gsl integral variable */
    gsl_integration_workspace * w_delta11;/**< Gsl integral variable */
    
    gsl_integration_workspace * w_Vp;/**< Gsl integral variable */
    gsl_integration_workspace * w_Vm;/**< Gsl integral variable */
    gsl_integration_workspace * w_Tp;/**< Gsl integral variable */
    gsl_integration_workspace * w_Tm;/**< Gsl integral variable */
    gsl_integration_workspace * w_V0;/**< Gsl integral variable */
    gsl_integration_workspace * w_T0;/**< Gsl integral variable */
    gsl_integration_workspace * w_S;/**< Gsl integral variable */
    
    gsl_integration_workspace * w_gtilde_1_re;/**< Gsl integral variable */
    gsl_integration_workspace * w_gtilde_1_im;/**< Gsl integral variable */
    gsl_integration_workspace * w_gtilde_2_re;/**< Gsl integral variable */
    gsl_integration_workspace * w_gtilde_2_im;/**< Gsl integral variable */
    gsl_integration_workspace * w_gtilde_3_re;/**< Gsl integral variable */
    gsl_integration_workspace * w_gtilde_3_im;/**< Gsl integral variable */
    
    gsl_integration_workspace * w_h_0_re;/**< Gsl integral variable */
    gsl_integration_workspace * w_h_0_im;/**< Gsl integral variable */
    gsl_integration_workspace * w_h_p_re;/**< Gsl integral variable */
    gsl_integration_workspace * w_h_p_im;/**< Gsl integral variable */
    gsl_integration_workspace * w_h_m_re;/**< Gsl integral variable */
    gsl_integration_workspace * w_h_m_im;/**< Gsl integral variable */
    
    unsigned int N_updated;/**< Cache variable */
    gslpp::vector<double> N_cache;/**< Cache variable */
    gslpp::complex Nc_cache;/**< Cache variable */
    
    unsigned int V_updated;/**< Cache variable */
    gslpp::vector<double> V_cache;/**< Cache variable */
    
    unsigned int A0_updated;/**< Cache variable */
    gslpp::vector<double> A0_cache;/**< Cache variable */
    
    unsigned int A1_updated;/**< Cache variable */
    gslpp::vector<double> A1_cache;/**< Cache variable */
    
    unsigned int T1_updated;/**< Cache variable */
    gslpp::vector<double> T1_cache;/**< Cache variable */
    
    unsigned int T2_updated;/**< Cache variable */
    gslpp::vector<double> T2_cache;/**< Cache variable */
    
    unsigned int k2_updated;/**< Cache variable */
    gslpp::vector<double> k2_cache;/**< Cache variable */
    
    unsigned int z_updated;/**< Cache variable */
    
    unsigned int lambda_updated;/**< Cache variable */
    
    unsigned int beta_updated;/**< Cache variable */
    double beta_cache;/**< Cache variable */
    
    unsigned int F_updated;/**< Cache variable */
    
    unsigned int VL1_updated;/**< Cache variable */
    unsigned int VL2_updated;/**< Cache variable */
    
    unsigned int TL1_updated;/**< Cache variable */
    unsigned int TL2_updated;/**< Cache variable */

    unsigned int VR1_updated;/**< Cache variable */
    unsigned int VR2_updated;/**< Cache variable */
    
    unsigned int TR1_updated;/**< Cache variable */
    unsigned int TR2_updated;/**< Cache variable */
    
    unsigned int VL0_updated;/**< Cache variable */
    gslpp::vector<double> VL0_cache;/**< Cache variable */
    
    unsigned int TL0_updated;/**< Cache variable */
    gslpp::vector<double> TL0_cache;/**< Cache variable */
    
    unsigned int VR0_updated;/**< Cache variable */
    
    unsigned int TR0_updated;/**< Cache variable */
    
    unsigned int Mb_Ms_updated;/**< Cache variable */
    
    unsigned int SL_updated;/**< Cache variable */
    gslpp::vector<double> SL_cache;/**< Cache variable */
    
    unsigned int SR_updated;/**< Cache variable */
    
    unsigned int C_1_updated;/**< Cache variable */
    gslpp::complex C_1_cache;/**< Cache variable */

    unsigned int C_2_updated;/**< Cache variable */
    gslpp::complex C_2_cache;/**< Cache variable */
    
    unsigned int C_3_updated;/**< Cache variable */
    gslpp::complex C_3_cache;/**< Cache variable */
    
    unsigned int C_4_updated;/**< Cache variable */
    gslpp::complex C_4_cache;/**< Cache variable */
    
    unsigned int C_5_updated;/**< Cache variable */
    gslpp::complex C_5_cache;/**< Cache variable */
    
    unsigned int C_6_updated;/**< Cache variable */
    gslpp::complex C_6_cache;/**< Cache variable */
    
    unsigned int C_7_updated;/**< Cache variable */
    gslpp::complex C_7_cache;/**< Cache variable */

    unsigned int C_9_updated;/**< Cache variable */
    gslpp::complex C_9_cache;/**< Cache variable */
    
    unsigned int C_10_updated;/**< Cache variable */
    gslpp::complex C_10_cache;/**< Cache variable */
    
    unsigned int C_7p_updated;/**< Cache variable */
    gslpp::complex C_7p_cache;/**< Cache variable */
    
    unsigned int C_9p_updated;/**< Cache variable */
    gslpp::complex C_9p_cache;/**< Cache variable */
    
    unsigned int C_10p_updated;/**< Cache variable */
    gslpp::complex C_10p_cache;/**< Cache variable */
    
    unsigned int C_S_updated;/**< Cache variable */
    gslpp::complex C_S_cache;/**< Cache variable */
    
    unsigned int C_P_updated;/**< Cache variable */
    gslpp::complex C_P_cache;/**< Cache variable */
    
    unsigned int C_Sp_updated;/**< Cache variable */
    gslpp::complex C_Sp_cache;/**< Cache variable */
    
    unsigned int C_Pp_updated;/**< Cache variable */
    gslpp::complex C_Pp_cache;/**< Cache variable */
    
    unsigned int Yupdated;/**< Cache variable */
    gslpp::vector<double> Ycache;/**< Cache variable */
    
    gslpp::complex h0Ccache[3];/**< Cache variable */
    gslpp::complex h1Ccache[3];/**< Cache variable */
    gslpp::complex h2Ccache[3];/**< Cache variable */
    
    unsigned int h0_updated;/**< Cache variable */
    unsigned int h1_updated;/**< Cache variable */
    unsigned int h2_updated;/**< Cache variable */
    
    unsigned int H_V0updated;/**< Cache variable */
    gslpp::vector<double> H_V0cache;/**< Cache variable */
    
    unsigned int H_V1updated;/**< Cache variable */
    gslpp::vector<double> H_V1cache;/**< Cache variable */
    
    unsigned int H_V2updated;/**< Cache variable */
    gslpp::vector<double> H_V2cache;/**< Cache variable */
    
    unsigned int H_A0updated;/**< Cache variable */
    unsigned int H_A1updated;/**< Cache variable */
    unsigned int H_A2updated;/**< Cache variable */
    
    unsigned int H_Supdated;/**< Cache variable */
    gslpp::vector<double> H_Scache;/**< Cache variable */
    
    unsigned int H_Pupdated;/**< Cache variable */
    gslpp::vector<double> H_Pcache;/**< Cache variable */
    
    unsigned int I0_updated;/**< Cache variable */
    unsigned int I1_updated;/**< Cache variable */
    unsigned int I2_updated;/**< Cache variable */
    unsigned int I3_updated;/**< Cache variable */
    unsigned int I4_updated;/**< Cache variable */
    unsigned int I5_updated;/**< Cache variable */
    unsigned int I6_updated;/**< Cache variable */
    unsigned int I7_updated;/**< Cache variable */
    unsigned int I8_updated;/**< Cache variable */
    unsigned int I9_updated;/**< Cache variable */
    unsigned int I10_updated;/**< Cache variable */
    unsigned int I11_updated;/**< Cache variable */
    
    unsigned int gtilde_1updated;/**< Cache variable */
    unsigned int gtilde_2updated;/**< Cache variable */
    unsigned int gtilde_3updated;/**< Cache variable */
    
    unsigned int h_0updated;/**< Cache variable */
    unsigned int h_pupdated;/**< Cache variable */
    unsigned int h_mupdated;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > sigma0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma1Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma2Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma3Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma4Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma5Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma6Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma7Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma9Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma10Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > sigma11Cached;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > delta0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta1Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta2Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta3Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta7Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > delta11Cached;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > VpCached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > TpCached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > VmCached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > TmCached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > V0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > T0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > SCached;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > gtilde_1Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > gtilde_2Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > gtilde_3Cached;/**< Cache variable */
    
    std::map<std::pair<double, double>, unsigned int > h_0Cached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > h_pCached;/**< Cache variable */
    std::map<std::pair<double, double>, unsigned int > h_mCached;/**< Cache variable */
    
};

#endif	/* MVLL_H */

