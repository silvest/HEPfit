/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDB2_H
#define	AMPDB2_H

class StandardModel;
#include "StandardModel.h"
#include "OrderScheme.h"
#include "gslpp.h"

/**
 * @addtogroup Flavour
 * @brief A module for all the flavour observables implemented in HEPfit.
 * @details This module has several flavour physics observables which include
 * quark flavour violation in the beauty, charm and strange sectors. This includes 
 * the evolutors, Hamiltonians and low energy observables.
 * @{
 */


/**
 * @class AmpDB2
 * @ingroup Flavour
 * @brief \f$ | \Delta B = 2 | \f$ Amplitude Class
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is related to the calculation of the \f$ B_{d,s}-\bar{B}_{d,s}\f$
 * mixing.
 *
 */

class AmpDB2 {
public:
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    AmpDB2(const StandardModel& SM_i, bool flag_RI = false);

    /**
    * @brief The value of @f$M_{21}^{bd}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$M_{21}^{bd}@f$
    */
    gslpp::complex getM21_Bd(orders order){
        return M21_Bd(order);
    }

    /**
    * @brief The value of @f$M_{21}^{bs}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$M_{21}^{bs}@f$
    */
    gslpp::complex getM21_Bs(orders order){
        return M21_Bs(order);
    }
    
    //mass schemes used in 2205.07907
    enum mass_schemes {pole, MSbar, PS, only1overmb};
    
    /**
    * @brief The value of @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$ from Gerlach (2205.07907 and thesis)
    * @param[in] order the %QCD order of the computation
    * @param[in] mass_scheme the scheme for the bottom quark mass
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$
    */
    gslpp::complex getGamma21overM21_Bd(orders order, mass_schemes mass_scheme = MSbar){
        return Gamma21overM21_Bd(order, mass_scheme);
    }

    /**
    * @brief The value of @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$ from Gerlach (2205.07907 and thesis)
    * @param[in] order the %QCD order of the computation
    * @param[in] mass_scheme the scheme for the bottom quark mass
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$
    */
    gslpp::complex getGamma21overM21_Bs(orders order, mass_schemes mass_scheme = MSbar){
        return Gamma21overM21_Bs(order, mass_scheme);
    }
    
    /**
    * @brief The value of @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$ in the traditional basis (hep-ph/0308029v2)
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$
    */
    gslpp::complex Gamma21overM21_BdFULLNLO_tradBasis();
    
    /**
    * @brief The value of @f$\frac{\Gamma_{21},M_{21}}^{bD}@f$ in the traditional basis (hep-ph/0308029v2)
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$
    */
    gslpp::complex Gamma21overM21_BsFULLNLO_tradBasis();

    /**
    * @brief The value of @f$\frac{\Gamma_{21},M_{21}}^{bD}@f$ in the traditional basis (hep-ph/0308029v2)
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$
    */
    gslpp::complex Gamma21overM21_BsLO_tradBasis();    

    gslpp::complex getPBd(){
        return PBd();
    }

    gslpp::complex getPBs(){
        return PBs();
    }
    
protected:
    /**
    * @brief A method to compute @f$M_{21}^{bd}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$M_{21}^{bd}@f$
    */
    gslpp::complex M21_Bd(orders order);
    
    /**
    * @brief A method to compute @f$M_{21}^{bs}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$M_{21}^{bs}@f$
    */
    gslpp::complex M21_Bs(orders order);

    /**
    * @brief A method to compute @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$
    * @detail source: Marvin Gerlach (2205.07907 and thesis) with 1/mb corrections from Lenz (hep-ph/0612167)
    * @param[in] order the %QCD order of the computation
    * @param[in] mass_scheme the scheme for the bottom quark mass
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$
    */
    gslpp::complex Gamma21overM21_Bd(orders order, mass_schemes mass_scheme);
    

    /**
    * @brief A method to compute @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$ 
    * @detail source: Marvin Gerlach (2205.07907 and thesis) with 1/mb corrections from Lenz (hep-ph/0612167)
    * @param[in] order the %QCD order of the computation
    * @param[in] mass_scheme the scheme for the bottom quark mass
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$
    */
    gslpp::complex Gamma21overM21_Bs(orders order, mass_schemes mass_scheme);
    
    
    /**
    * @brief A method to compute the ratio of the absolute value of the $B_s$ mixing amplitude over the Standard Model value.
    * @param[in] order the %QCD order of the computation
    * @return @f$\vert (M_{21}^{bs})_\mathrm{full}/(M_{21}^{bs})_\mathrm{SM}\vert@f$
    */
    gslpp::complex RBs(orders order);
    gslpp::complex PBd();
    gslpp::complex PBs();


private:
    const StandardModel& mySM;/**< Model type */
    
    gslpp::complex C_1_SM;/**<Wilson coeffients H_{eff}^{DF2} @f$C_1@f$*/
    
public:
    enum quark {d,s};   /*quark index i used for $B_i$*/
    enum quarks {cc, cu, uu};   /*combinations of u- and c- quarks in diagrams */
    
private:
    //mathematical constants
    const double zeta3 = gslpp_special_functions::zeta(3);
    const double log2 = log(2);
    const double log3 = log(3);
    const double sqrt3 = sqrt(3);        
    const double sqrt5 = sqrt(5);
    const double log12sqrt52 = log(0.5 + sqrt5/2.);    
    const double t_2 = -0.389012;   /*Gerlach thesis, eq. 3.102*/
    const double Cl2PI3 = 1.014941; /*Gerlach thesis, eq. 6.23*/
    
    double mu_1;        /*matching scale of DB=1 theory for leading order in 1/mb */
    double mu_1_overm;  /*matching scale of DB=1 theory for subleading order in 1/mb */
    double mu_2;        /*matching scale of DB=2 theory */
    double mu_b;        /*scale the running MSbar mass of the bottom quark */
    
    gslpp::vector<double> me = gslpp::vector<double>(5, 0.); /*DB=2 matrix elements in SUSY basis (arXiv:1907.01025v2) */
    gslpp::vector<double> me_R = gslpp::vector<double>(5, 0.); /*subleading DB=2 matrix elements R_0 to R_3 (Gerlach thesis) and R_4 (hep-ph/0308029v2) */
    gslpp::vector<double> me_Rtilde = gslpp::vector<double>(3, 0.); /*subleading DB=2 matrix elements R_1 to R_3 (Gerlach thesis) */
 
    
    
    //resummation to use z_bar instead of z and and eliminate z ln z terms (hep-ph/0612167)
    bool flag_resumz;
    
    //transformation matrix to switch to the RI scheme for the 5 matrix elements (hep-ph/0606197 eq. 5.10)
    gslpp::matrix<double> meMStoRI;
    //transformation matrix to switch to the RI scheme for the three DB=2 Wilson coefficients (hep-ph/0606197 eq. 5.10)
    gslpp::matrix<double> coeffsMStoRI;
    bool flag_RI; //flag to signal if transformation to RI is applied
    
    
    
    /**
    * @brief A method to compute CKM elements, quark masses and alpha_s
    * @param[in] order NNLO, use NLO only for tradBasis
    * @param[in] mass_scheme the scheme for the bottom quark mass
    */
    void computeCKMandMasses(orders order = NNLO, mass_schemes mass_scheme = MSbar);
 
    /**
    * @brief A method to compute all DB=2 Wilson coefficients (me, me_R, me_Rtilde)
    * @param[in] quark index of the neutral B mesons
    * @param[in] order %QCD orders used for R_0
    * @detail requires computeCKMandMasses() before use
    */
    void compute_matrixelements(quark q, orders order);


    //returns position in our array parameterization of the corresponding coefficient function
    int index_deltas(quarks qq, quark q);

    
    //often used values
    double Gf2;
    double z;
    double z2; //z^2
    double z3; //z^3
    double z4; //z^4
    double logz; //log(z)
    double log1minusz; //log(1-z)
    double log1minus4z; //log(1-4z)
    double oneminusz2; //(1 - z)^2
    double sqrt1minus4z; //(1-4z)^(1/2)
    double sigma; //(1 - sqrt1minus4z)/(1 + sqrt1minus4z)
    double logsigma; //log(sigma)
    double log2sigma; //log^2(sigma)
    double x_1; //mu_1/Mb
    double x_2; //mu_2/Mb
    double logx_1; //log(x_1)
    double logx_2; //log(x_2)
    double Dilogz; //Li_2(z)
    double Dilogsigma; //Li_2(sigma)
    double Dilogsigma2; //Li_2(sigma^2)
    const double M_PI2 = M_PI * M_PI;
    double as_4pi_mu1; //[alpha_s/(4Pi)](mu_1)
    double as_4pi_mu2; //[alpha_s/(4Pi)](mu_2)    
    double as_4pi; //[alpha_s/(4Pi)](mb(mb))
    
    // z for 1/mb corrections
    double z_1overm;
    double z_1overm2;  //z^2
    double oneminusz_1overm2;  //(1 - z)^2
    double sqrt1minus4z_1overm;  //(1-4z)^(1/2)
    
    double Md; //mass of the down quark in GeV
    double Ms; //mass of the strange quark in GeV
    double Mc; //mass of the charm quark in GeV   
    double Mb; //mass of the bottom quark in GeV
    double MB; //mass of the $B_d$ meson in GeV
    double MB_s; //mass of the $B_s$ meson in GeV
    double Mb2_prefactor; //overall Mb^2 prefactor of @f$\Gamma_{21}@f$
    double Mb_Mb; //MSbar mass of bottom
    double Mb_pole; //pole mass of bottom
    double Mb_PS; //PS mass of bottom 
    
    //parameters to calculate the bottom quark mass in the PS scheme (hep-ph/9804241)
    double mu_f = 2.;
    double K = 13.44 - 1.04 * 4.;
    double a1 = 31./3. - 10./9. * 4.;
    double b0 = 11. - 2. * 4. / 3.;
    
/**********                 DB=1 Wilson coefficients               ************/

    
    //Method to compute the DB=1 Wilson coefficients in the Buras basis to NNLO (arXiv:0401041)
    void computeWilsonCoeffsDB1bsg();
    
    //Method to compute the DB=1 Wilson coefficients in the Misiak basis to NNLO (hep-ph/9711280)    
    void computeWilsonCoeffsMisiak();
    
    gslpp::complex cacheC[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.};      /*FULLNNLO DB=1 Wilson coefficients C_i, i=1-6,8 */
    gslpp::complex cacheC_LO[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.};   /*LO DB=1 Wilson coefficients C_i, i=1-6,8 */
    gslpp::complex cacheC_NLO[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.};  /*NLO DB=1 Wilson coefficients C_i, i=1-6,8 */
    gslpp::complex cacheC_NNLO[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.}; /*NNLO DB=1 Wilson coefficients C_i, i=1-6,8 */
    gslpp::complex C(int i);    /*Value of the FULLNNLO DB=1 Wilson coefficients C_i, i=1-6,8 */
    
    
    
 /*******************************************************************************
 *  @f$\Gamma_{21}@f$ in NLO from Ciuchini (hep-ph/0308029v2)                   * 
 * ****************************************************************************/  
    
    //Values of the products of CKM elements
    gslpp::complex VtbVtd;
    gslpp::complex VtbVts;
    gslpp::complex VtbVtd2;
    gslpp::complex VtbVts2;
    gslpp::complex VcbVcd;
    gslpp::complex VcbVcs;
    gslpp::complex VcbVcd2;
    gslpp::complex VcbVcs2;

    
    //values of coefficient functions needed for DB=2 Wilson coefficients (hep-ph/0308029v2)
    double F0(quarks qq, int k, int i, int j);
    double F1(quarks qq, int k, int i, int j);
    double F(quarks qq, int k, int i, int j);
    double P(quarks qq, int k, int i, int j);
    gslpp::complex D(quarks qq, int k);
    
    double cacheF0[24] = { 0. };
    double cacheF1[24] = { 0. };
    double cacheP[84] = { 0. };
    gslpp::complex cacheD[6] = { 0. };
    
    //Methods to compute coefficient functions needed for DB=2 Wilson coefficients (hep-ph/0308029v2)
    void computeF0();
    void computeF1(); //requires "F0"
    void computeP();
    void computeD(); //requires "F" and "P"
    void computeD_LO(); //requires "F" and "P"

    //returns position in our array parameterization of the corresponding coefficient function
    int indexF(quarks qq, int k, int i, int j);
    int indexP(quarks qq, int k, int i, int j);
    int indexD(quarks qq, int k);
    
    /**
    * @brief Values of DB=2 Wilson coefficients (hep-ph/0308029v2)
    * @param[in] quark index of the neutral B mesons
    * @detail requires computeCKMandMasses() and "D"
    */
    gslpp::vector<gslpp::complex> c(quark q);
    
    
    /**
    * @brief Value of 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0308029v2)
    * @param[in] quark index of the neutral B mesons
    * @detail requires computeCKMandMasses() before use
    */
    gslpp::complex delta_1overm_tradBasis(quark q);
    
    //Values of the contributions to the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0308029v2)  
    gslpp::complex deltas_1overm_tradBasis(quarks qq, quark q);  //require computeCKMandMasses

    gslpp::complex cache_deltas_1overm_tradBasis[6] = { 0. };
    
    //Method to compute the contributions to the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0308029v2)  
    void compute_deltas_1overm_tradBasis(quark q); //require Wilson and computeCKMandMasses

    
    
 /*******************************************************************************
 *  @f$\Gamma_{21}@f$ in NNLO from Marvin Gerlach (2205.07907 and thesis)       * 
 * ****************************************************************************/
    
    gslpp::complex lambda_c; /* V_cq* V_cb  q=d,s */
    gslpp::complex lambda_u; /* V_uq* V_ub  q=d,s */
    
    const double M_PI4 = M_PI2 * M_PI2;
    bool orderofp[3] = {true, true, true}; /*signals if LO, NLO and NNLO contributions are used for DB=2 coefficients*/
    gslpp::vector<gslpp::complex> transformation(gslpp::vector< gslpp::complex > result, orders order);
    
    //Values of DB=2 Wilson coefficients (Gerlach thesis)
    gslpp::vector<gslpp::complex> c_H(); //require compute_pp_s and Wilson coefficients in Misiak basis
    gslpp::complex H(quarks qq, orders order); /*Values of contributions to the DB=2 Wilson coefficients for B_d (Gerlach thesis) */
    gslpp::complex H_s(quarks qq, orders order); /*Values of contributions to the DB=2 Wilson coefficients for B_s (Gerlach thesis) */

    //Values of DB=2 Wilson coefficients (Gerlach thesis) separated for
    //C-12-12 (LO, NLO, NNLO), C-12-36 (LO, NLO), C-36-36 (LO, NLO),C-12-8 (LO, NLO), C-36-8 (LO, NLO), C-8-8 (LO)
    gslpp::vector< gslpp::complex >  c_H_partial(int i);
    gslpp::vector<gslpp::complex> H_allpartial(quarks qq); /*Values of partial contributions to the DB=2 Wilson coefficients for B_d (Gerlach thesis) */
    gslpp::vector<gslpp::complex> H_s_allpartial(quarks qq); /*Values of partial contributions to the DB=2 Wilson coefficients for B_s (Gerlach thesis) */
    gslpp::complex H_partial(quarks qq, int i_start, int i_end, int j_start, int j_end, int n);
    gslpp::complex H_s_partial(quarks qq, int i_start, int i_end, int j_start, int j_end, int n); 
    
    //Values of the coefficient functions needed for DB=2 Wilson coefficients (Gerlach thesis)
    double p(quarks qq, int i, int j);
    double p_s(quarks qq, int i, int j);
    double p(quarks qq, int i, int j, int n);
    double p_s(quarks qq, int i, int j, int n);
    double lastInput_compute_pp_s[3] = {NAN, NAN, NAN};
    
    //Values of the coefficient functions needed for DB=2 Wilson coefficients (Gerlach thesis)
    double cache_p[576] = { 0. };
    double cache_ps[576] = { 0. };
    
    //Method to compute coefficient functions needed for DB=2 Wilson coefficients (Gerlach thesis)
    void compute_pp_s();
    
    //returns position in our array parameterization of p and p_s
    int index_p(quarks qq, int i, int j, int n);
    

    //A Method to adapt the DB=2 coefficient functions for the MSbar scheme (2106.05979 eq. (33))
    void poletoMSbar_pp_s();
    //constants from hep-ph/9912391v2  eq. (11)
    double PoletoMS_as1 = 4./3.;                                
    double PoletoMS_as2 = -(4. * (71./144. + M_PI2/18.) - 3019./288. + 1./6. * zeta3 - M_PI2/9. * log2 - M_PI2/3.);
    
    //A Method to adapt the DB=2 coefficient functions for the PS scheme (analog to 2106.05979 eq. (33))    
    void poletoPS_pp_s();
    //constants from hep-ph/9804241v2 eq. (21)
    double PoletoPS_as1;                
    double PoletoPS_as2;
    
    
    
 /*******************************************************************************
 *  1/mb corrections of @f$\Gamma_{21}@f$ from Lenz (hep-ph/0612167)            * 
 * ****************************************************************************/

    /**
    * @brief Value of 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
    * @param[in] quark index of the neutral B mesons
    * @detail requires compute_matrixelements() and Wilson coefficients in Buras basis
    */
    gslpp::complex delta_1overm(quark q);
    
    //Method to compute the contributions to the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
    void compute_deltas_1overm(quark q);
    
    //Values of the contributions to the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
    gslpp::complex deltas_1overm(quarks qq, quark q);
 
    gslpp::complex cache_deltas_1overm[6] = { 0. };
    
    //A method to compute the coefficients for the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
    void compute_g();
    
    //Values of the coefficients for the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
    gslpp::complex g(quarks qq, int i);
    gslpp::complex gtilde(quarks qq, int i);
    gslpp::complex cacheg[12] = { 0. };
    gslpp::complex cachegtilde[12] = { 0. };

    //returns position in our array parameterization of the corresponding coefficient function
    int indexg(quarks qq, int i);
    
    //LO DB=1 Wilson coefficients for 1/mb corrections
    gslpp::complex C_1LO;
    gslpp::complex C_2LO;
    //combinations of LO DB=1 Wilson coefficients
    gslpp::complex K_1; //3*C_1^2 + 2*C_1 * C_2
    gslpp::complex K_2; //C_2^2
};

/**
 * @}
 */

#endif	/* AMPDB2_H */

