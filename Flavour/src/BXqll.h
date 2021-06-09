/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BXQLL_H
#define	BXQLL_H

class StandardModel;

#include "QCD.h"
#include "ThObservable.h"
//#include "Particle.h"
#include "gslpp.h"
#include "Expanded.h"
#include <gsl/gsl_integration.h>
#include "HeffDF1.h"
#include "F_1.h"
#include "F_2.h"

/**
 * @class BXqll
 * @ingroup Flavour
 * @brief A class for the @f$B \to X_q l^+ l^-@f$ decay.  
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to 
 * build the observables relative to the @f$B \to X_q l^+ l^-@f$ decays, where
 * @f$X_q@f$ is an inclusive hadronic state containing a @f$q@f$ quark.
 * Formulae for the effective Hamiltonian are taken from Bobeth et al., hep-ph/9910220.
 * Formulae for the matrix elements are taken from Greub et al., arXiv:0810.4077
 */
class BXqll {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] quark_i quark defining the inclusive final hadronic state @f$X_q@f$ of the decay
     * @param[in] lep_i final leptons of the decay
     */
    BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~BXqll();  

    /**
    * @brief A method for initializing the parameters necessary for BXqll.
    * @return the vector of BXqll specific parameters
    */
    std::vector<std::string> initializeBXqllParameters();
    
    /**
     * @brief dGamma/ds for @f$B \to X_q l^+ l^-@f$ in the low dilepton invariant mass region.
     * @param[in] sh normalized dilepton invariant mass @f$q^2/m_b^2@f$
      */    
    double getR_LOWQ2(double sh);
    
    /**
    * @brief Method to obtain each observable \f$H_I\f$ as defined in @cite Huber:2015sra
    * @param[in] obs the observable in question
    * @param[in] sh normalized dilepton invariant mass @f$q^2/m_b^2@f$
    */
    double getH(std::string obs, double sh);
    
    /**
    * @brief The integral of each observable \f$H_I\f$ as defined in @cite Huber:2015sra
    * @param[in] obs the observable in question
    * @param[in] q_min minimum \f$q^2\f$ of the integral
    * @param[in] q_max maximum \f$q^2\f$ of the integral
    */
    double integrateH(std::string obs, double q_min, double q_max);

private:
    const StandardModel& mySM;/**< Model type */
    F_1 myF_1;
    F_2 myF_2;
    HeffDF1 myHeff;
    QCD::lepton lep;/**< Final leptons type */
    QCD::quark quark;/**< Initial meson type */
    double CF, GF, ale, alsmu, alsmuc, alstilde, aletilde, kappa;
    double Mlep, mu_b, mu_c, Mb, Mc, Mtau, Mb_pole, Mc_pole, Ms, MW;
    double abslambdat_over_Vcb, Vts_over_Vcb, z, muh, lambda_1, lambda_2, Lbl;
    double phi1, phi2, phi00, phi00_2, phi20, phi01;
    double phinv00, phinv10, phinv20, phinv01, phinv11, phinv21, phinv_00_01;
    double BR_BXcenu, C_ratio, pre;
    unsigned int QCD_max, QED_max;
   
    std::vector< gslpp::vector<gslpp::complex> > M_7;
    std::vector< gslpp::vector<gslpp::complex> > M_9;
    std::vector< gslpp::vector<gslpp::complex> > M_10;
    
    std::vector< gslpp::matrix<gslpp::complex> > Hij_T;
    std::vector< gslpp::matrix<gslpp::complex> > Hij_L;
    std::vector< gslpp::matrix<gslpp::complex> > Hij_A;

    std::vector<std::string> BXqllParameters;/**< The string of mandatory BXqll parameters */

//    gslpp::vector<gslpp::complex> ** allcoeff_smm;/**<Vector that contains the Wilson coeffients */
    Expanded<gslpp::vector<gslpp::complex> > allcoeff;/**<Vector that contains the Wilson coeffients */
    
    gslpp::matrix<gslpp::complex> WC;/**<Matrix that contains the Wilson coeffients for each order */
    
    double aveH;/**< Gsl integral variable */
    double errH;/**< Gsl integral variable */
    gsl_function FH;/**< Gsl integral variable */
    gsl_integration_cquad_workspace * w_H;/**< Gsl integral variable */
    gsl_error_handler_t * old_handler; /**< GSL error handler store */

    /**
     * @brief The update parameter method for BXqll.
     */
    void updateParameters();

    double F_17re(double muh, double z, double sh, int maxpow=20);
    double F_17im(double muh, double z, double sh, int maxpow=20);
    double F_19re(double muh, double z, double sh, int maxpow=20);
    double F_19im(double muh, double z, double sh, int maxpow=20);
    double F_27re(double muh, double z, double sh, int maxpow=20);
    double F_27im(double muh, double z, double sh, int maxpow=20);
    double F_29re(double muh, double z, double sh, int maxpow=20);
    double F_29im(double muh, double z, double sh, int maxpow=20);
    
    double DeltaF_19re(double muh, double z, double sh, int maxpow=20);
    double DeltaF_19im(double muh, double z, double sh, int maxpow=20);
    double DeltaF_29re(double muh, double z, double sh, int maxpow=20);
    double DeltaF_29im(double muh, double z, double sh, int maxpow=20);
    /**
    * @brief The correction \f$ F_{17} \f$ from @cite Greub:2008cy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$ F_{17} \f$
    */
    gslpp::complex F17(double sh);

    /**
    * @brief The correction \f$ F_{27} \f$ from @cite Greub:2008cy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$ F_{27} \f$
    */
    gslpp::complex F27(double sh);

    /**
    * @brief The correction \f$ F_{19} \f$ from @cite Greub:2008cy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$ F_{19} \f$
    */
    gslpp::complex F19(double sh);

    /**
    * @brief The correction \f$ F_{29} \f$ from @cite Greub:2008cy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$ F_{29} \f$
    */
    gslpp::complex F29(double sh);

    /**
    * @brief The correction \f$ F_{87} \f$ from @cite Greub:2008cy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$ F_{87} \f$
    */
    gslpp::complex F87(double sh);

    /**
    * @brief The correction \f$ F_{89} \f$ from @cite Greub:2008cy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$ F_{89} \f$
    */
    double F89(double sh);
    
    /**
    * @brief Auxiliary functions \f$S_{NM}^T\f$ from @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    * @param[in] order LO or NLO
    */
    double S77_T(double sh, orders order);
    double S79_T(double sh, orders order);
    double S99_T(double sh, orders order);
    double S1010_T(double sh, orders order);
    
    /**
    * @brief Auxiliary functions \f$S_{NM}^L\f$ from @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    * @param[in] order LO or NLO
    */
    double S77_L(double sh, orders order);
    double S79_L(double sh, orders order);
    double S99_L(double sh, orders order);
    double S1010_L(double sh, orders order);
    
    /**
    * @brief Auxiliary functions \f$S_{NM}^A\f$ from @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    * @param[in] order LO or NLO
    */
    double S710_A(double sh, orders order);
    double S910_A(double sh, orders order);
    
    /**
    * @brief \f$\mathcal{O}(\Lambda_{QCD}^2/m_c^2)\f$ contributions \f$c_{ij}^I\f$ as defined in @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    * @param[in] i,j indices in eq. (4.10) in @cite Huber:2015sra
    * @param[in] ord possible QCD,QED orders: 11, 22, 32
    */
    gslpp::complex cij_T(unsigned int i, unsigned int j, double sh, unsigned int ord);
    gslpp::complex cij_L(unsigned int i, unsigned int j, double sh, unsigned int ord);
    gslpp::complex cij_A(unsigned int i, unsigned int j, double sh);
    
    /**
    * @brief Log-enhanced electromagnetic corrections \f$e_{ij}^I\f$ as defined in @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    * @param[in] i,j indices in eq. (4.10) in @cite Huber:2015sra
    */
    gslpp::complex eij_T(unsigned int i, unsigned int j, double sh);
    gslpp::complex eij_L(unsigned int i, unsigned int j, double sh);
    gslpp::complex eij_A(unsigned int i, unsigned int j, double sh);
    
    /**
    * @brief Auxiliary functions \f$omega_{NM}^T\f$ from @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double omega77_T(double sh);
    double omega79_T(double sh);
    double omega99_T(double sh);
    
    /**
    * @brief Auxiliary functions \f$omega_{NM}^L\f$ from @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double omega77_L(double sh);
    double omega79_L(double sh);
    double omega99_L(double sh);
    
    /**
    * @brief Auxiliary functions \f$omega_{NM}^A\f$ from @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double omega710_A(double sh);
    double omega910_A(double sh);
    
    /**
    * @brief Auxiliary functions \f$omega_{NM,T}^{(em)}\f$ from @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double omega77em_T(double sh);
    double omega79em_T(double sh);
    double omega99em_T(double sh);
    double omega22em_T(double sh);
    gslpp::complex omega27em_T(double sh);
    gslpp::complex omega29em_T(double sh);
    
    /**
    * @brief Auxiliary functions \f$omega_{NM,L}^{(em)}\f$ from @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double omega77em_L(double sh);
    double omega79em_L(double sh);
    double omega99em_L(double sh);
    double omega22em_L(double sh);
    gslpp::complex omega27em_L(double sh);
    gslpp::complex omega29em_L(double sh);
    
    /**
    * @brief Auxiliary functions \f$omega_{NM}^A\f$ from @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double omega710em_A(double sh);
    double omega910em_A(double sh);
    gslpp::complex omega210em_A(double sh);
    
    /**
    * @brief Auxiliary function \f$f_{i}\f$ from @cite Huber:2005ig
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    * @param[in] gamma_9 anomalous dimension matrix \f$gamma_{i9^{(01)}}\f$
    * @param[in] rho_c,b,0,num numbers from Table 7 of @cite Huber:2005ig
    */
    gslpp::complex f_Huber(double sh, double gamma_9, double rho_c, double rho_b, double rho_0, double rho_num);
    
    /**
    * @brief Auxiliary function \f$f_{9}^{pen}\f$ from @cite Huber:2005ig
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    gslpp::complex f9pen_Huber(double sh);
    
    /**
    * @brief Auxiliary function \f$g(y)\f$ from @cite Huber:2005ig
    * @param[in] y fraction of z over sh
    */
    gslpp::complex g_Huber(double y);
    
    /**
    * @brief Kruger-Sehgal factorizable non-perturbative charm contributions following @cite Huber:2007vv 
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    gslpp::complex KS_cc(double sh);
    
    /**
    * @brief Auxiliary function for the Kruger-Sehgal charm contributions
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    * @param[in] m charmonium mass
    * @param[in] Gamma charmonium total decay width
    * @param[in] Br_ll branching fraction of the decay mode \f$V \to l^+ l^-\f$
    * @param[in] Br_had branching fraction of the decay mode \f$V \to {\rm hadrons}\f$
    */
    gslpp::complex KS_aux(double sh, double m, double Gamma, double Br_ll, double Br_had);
    
    /**
    * @brief Auxiliary function \f$F(r)\f$ from @cite Buchalla:1997ky
    * @param[in] r normalized dilepton invariant mass \f$q^2/{4 m_c^2}\f$
    */
    gslpp::complex F_BIR(double r);
    
    /**
    * @brief Vectors of auxiliary functions \f$M_i^{7,9,10}sh)\f$ from Table 6 of @cite Huber:2005ig
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    * @param[in] order LO or NLO
    */
    void computeMi(double sh);
    
    /**
    * @brief Angular observable \f$H_T\f$ as defined in @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double H_T (double sh);
    
    /**
    * @brief Angular observable \f$H_L\f$ as defined in @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double H_L (double sh);
    
    /**
    * @brief Angular observable \f$H_A\f$ as defined in @cite Huber:2015sra
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double H_A (double sh);

    /**
    * @brief Normalization function for \f$B\to X_s\ell\ell\f$ from eq. (4.8) of 1503.04849
    * @param[in] ord/ord_qed order to be returned
    */
    double Phi_u(orders ord);
    double Phi_u(orders_qed ord_qed);
    
    /**
    * @brief Inverse of the normalization function for \f$B\to X_s\ell\ell\f$ from eq. (4.8) of 1503.04849
    * @param[in] (ord_qcd, ord_qed) order to be returned
    */
    double Phi_u_inv(unsigned int ord_qcd, unsigned int ord_qed);
    
    /**
    * @brief The finite bremsstrahlung corrections to dGamma/ds for @f$B \to X_q l^+ l^-@f$ from @cite Asatryan:2002iy
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double PhiTL_brems(double sh);
    
    /**
    * @brief The finite bremsstrahlung corrections to dAFB/ds for @f$B \to X_q l^+ l^-@f$ from @cite Asatryan:2003yk
    * @param[in] sh normalized dilepton invariant mass \f$q^2/m_b^2\f$
    */
    double PhiA_brems(double sh);
    
    /**
    * @brief The modified coefficient \f${\tilde C}_9^{0,\mathrm{mod}}\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f${\tilde C}_9^{0,\mathrm{mod}}\f$
    */
    gslpp::complex C9mod(double sh);
    
    /**
    * @brief Auxiliary function \f$h(z,sh)\f$ from @cite Asatrian:2001zw.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @param[in] z \f$m_c^2/m_b^2\f$
    * @return \f$h(z,sh)\f$
    */
    gslpp::complex h_z(double zed, double sh);
    
    /**
    * @brief The finite bremsstrahlung correction \f$tau_{78}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{78}(\hat s)\f$
    */
    double tau78(double sh);
    
    /**
    * @brief The finite bremsstrahlung correction \f$tau_{89}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{89}(\hat s)\f$
    */
    double tau89(double sh);
    
    /**
    * @brief The finite bremsstrahlung correction \f$tau_{88}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{88}(\hat s)\f$
    */
    double tau88(double sh);
    
    /**
    * @brief The fit of the finite bremsstrahlung correction \f$tau_{22}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{22}(\hat s)\f$
    */
    double tau22fit(double sh);
    
    /**
    * @brief The fit of the real part of finite bremsstrahlung correction \f$tau_{27}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{27}(\hat s)\f$
    */
    double tau27fit_Re(double sh);
    
    /**
    * @brief The fit of the imaginary part of finite bremsstrahlung correction \f$tau_{27}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{27}(\hat s)\f$
    */
    double tau27fit_Im(double sh);
    
    /**
    * @brief The fit of the real part of finite bremsstrahlung correction \f$tau_{28}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{28}(\hat s)\f$
    */
    double tau28fit_Re(double sh);
    
    /**
    * @brief The fit of the imaginary part of finite bremsstrahlung correction \f$tau_{28}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{28}(\hat s)\f$
    */
    double tau28fit_Im(double sh);
    
    /**
    * @brief The fit of the real part of finite bremsstrahlung correction \f$tau_{29}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{29}(\hat s)\f$
    */
    double tau29fit_Re(double sh);
    
    /**
    * @brief The fit of the imaginary part of finite bremsstrahlung correction \f$tau_{29}(\hat s)\f$ from @cite Asatryan:2002iy.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{29}(\hat s)\f$
    */
    double tau29fit_Im(double sh);
    
    /**
    * @brief The finite bremsstrahlung correction \f$f_{710}(\hat s)\f$ from @cite Asatrian:2002va.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{710}(\hat s)\f$
    */
    double f710(double sh);
    
    /**
    * @brief The finite bremsstrahlung correction \f$f_{910}(\hat s)\f$ from @cite Asatrian:2002va.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{910}(\hat s)\f$
    */
    double f910(double sh);
    
    /**
    * @brief The finite bremsstrahlung correction \f$f_{810}(\hat s)\f$ from @cite Asatrian:2003yk.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{810}(\hat s)\f$
    */
    double t810(double sh);
    
    /**
    * @brief The fit of finite bremsstrahlung correction \f$t_{210}(\hat s)\f$ from @cite Asatrian:2003yk.
    * @param[in] sh \f$q^2/m_b^2\f$ of the decay
    * @return \f$tau_{210}(\hat s)\f$
    */
    double t210fit(double sh);
    
    /**
    * @brief Auxiliary function that matches orders_qed to an integer
    * @param[in] ord_qed order to be returned
    */
    unsigned int int_qed(orders_qed order_qed);
    
    /**
    * @brief Auxiliary function that performs the multiplication of Wilson coefficients and matrix elements
    * @param[in] Hij matrix element related to the one of the angular observables of \f$B\to X_s\ell\ell\f$
    */
    double CCH_multiplication(std::vector< gslpp::matrix<gslpp::complex> >& Hij);
    
    /**
    * @brief Temporary method to test Wilson coefficients with C10_OS1 matching and HeffDF1 evolution 
    */
    void Test_WC_DF1();
    
    friend double gslpp_special_functions::dilog(double x);
};
#endif	/* BXqLL_H */
