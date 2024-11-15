/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDB2_H
#define AMPDB2_H

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

class AmpDB2
{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] BMeson_i an integer to specify the B meson: 0 for Bd and 1 for Bs
     */
    AmpDB2(const StandardModel &SM_i, int BMeson_i, bool flag_fixmub = false, bool flag_RI = false);

    /**
     * @brief The value of @f$M_{21}@f$ for @f$B_{d,s}@f$ mesons.
     * @param[in] order the %QCD order of the computation
     * @return @f$M_{21}@f$
     */
    gslpp::complex getM21(orders order)
    {
        if (BMeson == 0)
        {
            return M21_Bd(order);
        }
        else
        {
            return M21_Bs(order);
        }
    }

    enum mass_schemes {pole, MSbar, PS, only1overmb, MSbar_partialNNLO, PS_partialNNLO, MSbar_partialN3LO, PS_partialN3LO}; //mass schemes
    enum quark {d,s};   /*quark index i used for $B_i$*/
    enum quarks {cc, cu, uu};   /*combinations of u- and c- quarks in diagrams */
    
    /**
     * @brief The value of @f$\frac{\Gamma_{21},M_{21}}@f$ in the traditional basis
     * @param[in] order the %QCD order of the computation
     * @param[in] mass_scheme the scheme for the bottom quark mass
     * @return @f$\frac{\Gamma_{21},M_{21}}@f$
     */
    gslpp::complex getGamma21overM21_tradBasis(orders order)
    {
        if (BMeson == 0)
        {
            return Gamma21overM21_tradBasis(order, d);
        }
        else
        {
            return Gamma21overM21_tradBasis(order, s);
        }
    }

    /**
     * @brief The value of @f$\frac{\Gamma_{21},M_{21}}@f$ from Gerlach (2205.07907 and thesis)
     * @param[in] order the %QCD order of the computation
     * @param[in] mass_scheme the scheme for the bottom quark mass
     * @return @f$\frac{\Gamma_{21},M_{21}}@f$
     */
    gslpp::complex getGamma21overM21(orders order, mass_schemes mass_scheme = MSbar)
    {
        return Gamma21overM21(order, mass_scheme, BMeson);
    }

    gslpp::complex getPB()
    {
        if (BMeson == 0)
        {
            return PBd();
        }
        else
        {
            return PBs();
        }
    }
    gslpp::complex getRB(orders order)
    {
        if (BMeson == 0)
        {
            return RBd(order);
        }
        else
        {
            return RBs(order);
        }
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
     * @brief A method to compute @f$\frac{\Gamma_{21},M_{21}}^{bq}@f$
     * @detail source: Marvin Gerlach (2205.07907 and thesis) with 1/mb corrections from Lenz (hep-ph/0612167)
     * @param[in] order the %QCD order of the computation
     * @param[in] mass_scheme the scheme for the bottom quark mass
     * @param[in] BMeson the B meson index (0 for Bd and 1 for Bs)
     * @return @f$\frac{\Gamma_{21},M_{21}}^{bq}@f$
     */
    gslpp::complex Gamma21overM21(orders order, mass_schemes mass_scheme, int BMeson);

    /**
     * @brief A method to compute @f$\frac{\Gamma_{21},M_{21}}^{bq}@f$ in the traditional basis
     * @detail source: Ciuchini (hep-ph/0308029v2)
     * @param[in] order the %QCD order of the computation
     * @param[in] BMeson the B meson index (0 for Bd and 1 for Bs)
     * @return @f$\frac{\Gamma_{21},M_{21}}^{bq}@f$
     */
    gslpp::complex Gamma21overM21_tradBasis(orders order, quark q);

    /**
     * @brief A method to compute the ratio of the absolute value of the $B_s$ mixing amplitude over the Standard Model value.
     * @param[in] order the %QCD order of the computation
     * @return @f$\vert (M_{21}^{bs})_\mathrm{full}/(M_{21}^{bs})_\mathrm{SM}\vert@f$
     */
    gslpp::complex RBs(orders order);
     /**
     * @brief A method to compute the ratio of the absolute value of the $B_d$ mixing amplitude over the Standard Model value.
     * @param[in] order the %QCD order of the computation
     * @return @f$\vert (M_{21}^{bd})_\mathrm{full}/(M_{21}^{bd})_\mathrm{SM}\vert@f$
     */
    gslpp::complex RBd(orders order);
    gslpp::complex PBd();
    gslpp::complex PBs();

private:
    const StandardModel &mySM; /**< Model type */

    int BMeson; /**< B meson index (0 for Bd and 1 for Bs) */

    gslpp::complex C_1_SM; /**<Wilson coeffients H_{eff}^{DF2} @f$C_1@f$*/

    // mathematical constants
private:
    //mathematical constants
    const double M_PI2 = M_PI * M_PI;
    const double M_PI3 = M_PI2 * M_PI;
    const double M_PI4 = M_PI2 * M_PI2;
    const double zeta2 = gslpp_special_functions::zeta(2);
    const double zeta3 = gslpp_special_functions::zeta(3);
    const double zeta4 = gslpp_special_functions::zeta(4);
    const double zeta5 = gslpp_special_functions::zeta(5);
    const double log2 = log(2);
    const double log2_2 = log2 * log2;
    const double log2_4 = log2_2 * log2_2;
    const double log3 = log(3);
    const double sqrt3 = sqrt(3);
    const double sqrt5 = sqrt(5);
    const double log12sqrt52 = log(0.5 + sqrt5 / 2.);
    const double t_2 = -0.389011713;   // Im(Dilog((3 - i*sqrt(3))/6)
    const double Cl2PI3 = 1.014941606; // Clausen(2, Pi/3)
    const double polylog4_12 = 0.517479062; //PolyLog[4, 1/2]
    
    double mu_1;        /*matching scale of DB=1 theory for leading order in 1/mb */
    double mu_1_overm;  /*matching scale of DB=1 theory for subleading order in 1/mb */
    double mu_2;        /*matching scale of DB=2 theory */
    double mu_b;        /*scale the running MSbar mass of the bottom quark */
    
    gslpp::vector<double> me = gslpp::vector<double>(5, 0.); /*DB=2 matrix elements in SUSY basis (arXiv:1907.01025v2) */
    gslpp::vector<double> meoverme0 = gslpp::vector<double>(3, 0.); /*DB=2 matrix elements me(1),me(2),me(3) */
    gslpp::vector<double> me_R = gslpp::vector<double>(5, 0.);      /*subleading DB=2 matrix elements R_0 to R_3 (Gerlach thesis) and R_4 (hep-ph/0308029v2) */
    gslpp::vector<double> me_Rtilde = gslpp::vector<double>(3, 0.); /*subleading DB=2 matrix elements R_1 to R_3 (Gerlach thesis) */

    // resummation to use z_bar instead of z and and eliminate z ln z terms (hep-ph/0612167)
    bool flag_resumz;

    // transformation matrix to switch to the RI scheme for the 5 matrix elements (hep-ph/0606197 eq. 5.10)
    gslpp::matrix<double> meMStoRI;
    // transformation matrix to switch to the RI scheme for the three DB=2 Wilson coefficients (hep-ph/0606197 eq. 5.10)
    gslpp::matrix<double> coeffsMStoRI;
    bool flag_fixmub; // flag to fix mu_b=mu_c to 4.2 GeV
    bool flag_RI;     // flag to signal if transformation to RI is applied

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

    // returns position in our array parameterization of the corresponding coefficient function
    int index_deltas(quarks qq, quark q);

    // often used values
    double Gf2;
    double z;
    double z2; //z^2
    double z3; //z^3
    double z4; //z^4
    double sqrtz; //sqrt(z)
    double logz; //log(z)
    double log2z; //log(z^2)
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
    double as_4pi_mu1; //[alpha_s/(4Pi)](mu_1)
    double as_4pi_mu2; //[alpha_s/(4Pi)](mu_2)
    double as_4pi;     //[alpha_s/(4Pi)](mb(mb))

    // z for 1/mb corrections
    double z_1overm;
    double z_1overm2;           // z^2
    double oneminusz_1overm2;   //(1 - z)^2
    double sqrt1minus4z_1overm; //(1-4z)^(1/2)

    double Md;                   // mass of the down quark in GeV
    double Ms;                   // mass of the strange quark in GeV
    double Mc;                   // mass of the charm quark in GeV
    double Mb;                   // mass of the bottom quark in GeV
    double MB;                   // mass of the $B_d$ meson in GeV
    double MB_s;                 // mass of the $B_s$ meson in GeV
    double Mb2_prefactor;        // overall Mb^2 prefactor of @f$\Gamma_{21}@f$
    double Mb2_prefactor_1overm; // Mb^2 prefactor of the 1/mb part of @f$\Gamma_{21}@f$
    double Mb_Mb;                // MSbar mass of bottom
    double Mb_pole;              // pole mass of bottom
    double Mb_PS;                // PS mass of bottom

    // parameters to calculate the bottom quark mass in the PS scheme (hep-ph/9804241)
    double mu_f = 2.;
    double nl = 4.;
    double a1 = 31./3. - 10./9. * nl;
    double a2 = (4343./162. + 6.*M_PI2 - M_PI4/4. + 22./3. * zeta3) * 3. * 3. - (1798./81. + 56./3. * zeta3) * 3. * 0.5 * nl
         - (55./3. - 16. * zeta3) * 4./3. * 0.5 * nl + 400./81. * 0.25 * nl * nl;
    double b0 = 11. - 2. * nl/3.;
    double b0_2 = b0 * b0;
    double b1 = 102. - 38. * nl/3.;
    
/**********                 DB=1 Wilson coefficients               ************/

    
    //Method to compute the DB=1 Wilson coefficients in the Buras basis to NNLO (arXiv:0401041)
    void computeWilsonCoeffsBuras();
    
    //Method to compute the DB=1 Wilson coefficients in the Misiak basis to NNLO (hep-ph/9711280)    
    void computeWilsonCoeffsMisiak();
    
    gslpp::complex cacheC[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.};      /*FULLNNLO DB=1 Wilson coefficients C_i, i=1-6,8 */
    gslpp::complex C_Misiak_LO[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.};   /*LO DB=1 Wilson coefficients in Misiak basis C_i, i=1-6,8 */
    gslpp::complex C_Misiak_NLO[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.};  /*NLO DB=1 Wilson coefficients in Misiak basis C_i, i=1-6,8 */
    gslpp::complex C_Misiak_NNLO[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.}; /*NNLO DB=1 Wilson coefficients in Misiak basis C_i, i=1-6,8 */
    gslpp::complex C_Buras_LO[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.};   /*LO DB=1 Wilson coefficients in Buras basis C_i, i=1-6,8 */
    gslpp::complex C_Buras_NLO[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.};  /*NLO DB=1 Wilson coefficients in Buras basis C_i, i=1-6,8 */
    gslpp::complex C_Buras_NNLO[8] = { 0., 0., 0., 0., 0., 0., NAN, 0.}; /*NNLO DB=1 Wilson coefficients in Buras basis C_i, i=1-6,8 */

    
    
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

    // values of coefficient functions needed for DB=2 Wilson coefficients (hep-ph/0308029v2)
    double F0(quarks qq, int k, int i, int j);
    double F1(quarks qq, int k, int i, int j);
    double F(quarks qq, int k, int i, int j);
    double P(quarks qq, int k, int i, int j);
    gslpp::complex D(quarks qq, int k);

    double cacheF0[24] = {0.};
    double cacheF1[24] = {0.};
    double cacheP[84] = {0.};
    gslpp::complex cacheD[6] = {0.};

    // Methods to compute coefficient functions needed for DB=2 Wilson coefficients (hep-ph/0308029v2)
    void computeF0();
    void computeF1(); // requires "F0"
    void computeP();
    void computeD(orders order); // requires "F" and "P"

    // returns position in our array parameterization of the corresponding coefficient function
    int indexF(quarks qq, int k, int i, int j);
    int indexP(quarks qq, int k, int i, int j);
    int indexD(quarks qq, int k);

    /**
    * @brief Values of DB=2 Wilson coefficients
    * from (hep-ph/0308029v2) transformed to the new basis
    * @param[in] quark index of the neutral B mesons
    * @param[in] order the %QCD order for the transformations
    * @detail requires computeCKMandMasses() and "D"
    */
    gslpp::vector<gslpp::complex> c(quark q, orders order);

    /**
     * @brief Value of 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0308029v2)
     * @param[in] quark index of the neutral B mesons
     * @detail requires computeCKMandMasses() before use
     */
    gslpp::complex delta_1overm_tradBasis(quark q);

    // Values of the contributions to the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0308029v2)
    gslpp::complex deltas_1overm_tradBasis(quarks qq, quark q); // require computeCKMandMasses

    gslpp::complex cache_deltas_1overm_tradBasis[6] = {0.};

    // Method to compute the contributions to the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0308029v2)
    void compute_deltas_1overm_tradBasis(quark q); // require Wilson and computeCKMandMasses

    /*******************************************************************************
     *  @f$\Gamma_{21}@f$ in NNLO from Marvin Gerlach (2205.07907 and thesis)       *
     * ****************************************************************************/

    // Values of the products of CKM elements
    gslpp::complex lambda_c_d; /* V_cd* V_cb  */
    gslpp::complex lambda_u_d; /* V_ud* V_ub  */
    gslpp::complex lambda_c_s; /* V_cs* V_cb  */
    gslpp::complex lambda_u_s; /* V_us* V_ub  */
    
    gslpp::vector<gslpp::complex> transformation(gslpp::vector< gslpp::complex > result, orders order);
    
    //Values of DB=2 Wilson coefficients (Gerlach thesis)
    gslpp::vector<gslpp::complex> c_H(quark q, orders order); //require compute_pp_s and Wilson coefficients in Misiak basis
    gslpp::complex H(quarks qq, orders order); /*Values of contributions to the DB=2 Wilson coefficients for B_d (Gerlach thesis) */
    gslpp::complex H_s(quarks qq, orders order); /*Values of contributions to the DB=2 Wilson coefficients for B_s (Gerlach thesis) */

    // Values of DB=2 Wilson coefficients (Gerlach thesis) separated for
    // C-12-12 (LO, NLO, NNLO), C-12-36 (LO, NLO), C-36-36 (LO, NLO),C-12-8 (LO, NLO), C-36-8 (LO, NLO), C-8-8 (LO)
    gslpp::vector<gslpp::complex> c_H_partial(quark q, int i);
    gslpp::vector<gslpp::complex> H_allpartial(quarks qq);   /*Values of partial contributions to the DB=2 Wilson coefficients for B_d (Gerlach thesis) */
    gslpp::vector<gslpp::complex> H_s_allpartial(quarks qq); /*Values of partial contributions to the DB=2 Wilson coefficients for B_s (Gerlach thesis) */
    gslpp::complex H_partial(quarks qq, int i_start, int i_end, int j_start, int j_end, int n);
    gslpp::complex H_s_partial(quarks qq, int i_start, int i_end, int j_start, int j_end, int n);

    // Values of the coefficient functions needed for DB=2 Wilson coefficients (Gerlach thesis)
    double p(quarks qq, int i, int j, int n, bool flag_LOz = false);
    double p_s(quarks qq, int i, int j, int n, bool flag_LOz = false);
    //double lastInput_compute_pp_s[4] = {NAN, NAN, NAN, NAN};
    
    //Values of the coefficient functions needed for DB=2 Wilson coefficients (Gerlach thesis)
    double cache_p[768] = { 0. };
    double cache_ps[768] = { 0. };
    //Values of the coefficient functions in LO in z needed for DB=2 Wilson coefficients (Gerlach thesis)
    bool flag_LOz = true;
    double cache_p_LO[576] = {0.};
    double cache_ps_LO[576] = {0.};

    // Method to compute coefficient functions needed for DB=2 Wilson coefficients (Gerlach thesis)
    void compute_pp_s();

    // returns position in our array parameterization of p and p_s
    int index_p(quarks qq, int i, int j, int n);

    // A Method to adapt the DB=2 coefficient functions for the MSbar scheme (2106.05979 eq. (33))
    void poletoMSbar_pp_s();
    void poletoMSbar_pp_s_partialN3LO();
    //constants from hep-ph/9912391v2  eq. (11)
    double PoletoMS_as1;                                
    double PoletoMS_as2;
    double PoletoMS_as3;
    double PoletoMS_as2_z0; //0th order in z
    double PoletoMS_as2_z1; //1st order in z
    
    //A Method to adapt the DB=2 coefficient functions for the PS scheme (analog to 2106.05979 eq. (33))    
    void poletoPS_pp_s();
    void poletoPS_pp_s_partialN3LO();
    //constants from hep-ph/9804241v2 eq. (21)
    double PoletoPS_as1;                
    double PoletoPS_as2;
    double PoletoPS_as3;

    //A Method to discard true NNLO contributions from the DB=1 and DB=2 Wilson coefficients (like with partialN3LO)
    void compute_partialNNLO();
    
    
 /*******************************************************************************
 *  1/mb corrections of @f$\Gamma_{21}@f$ from Lenz (hep-ph/0612167)            * 
 * ****************************************************************************/

    /**
     * @brief Value of 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
     * @param[in] quark index of the neutral B mesons
     * @detail requires compute_matrixelements() and Wilson coefficients in Buras basis
     */
    gslpp::complex delta_1overm(quark q);

    // Method to compute the contributions to the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
    void compute_deltas_1overm(quark q);

    // Values of the contributions to the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
    gslpp::complex deltas_1overm(quarks qq, quark q);

    gslpp::complex cache_deltas_1overm[6] = {0.};

    // A method to compute the coefficients for the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
    void compute_g();

    // Values of the coefficients for the 1/mb corrections of @f$\Gamma_{21}@f$ (hep-ph/0612167)
    gslpp::complex g(quarks qq, int i);
    gslpp::complex gtilde(quarks qq, int i);
    gslpp::complex cacheg[12] = {0.};
    gslpp::complex cachegtilde[12] = {0.};

    // returns position in our array parameterization of the corresponding coefficient function
    int indexg(quarks qq, int i);
    
    //LO DB=1 Wilson coefficients for 1/mb corrections
    gslpp::complex C1_LO_1overm;
    gslpp::complex C2_LO_1overm;
    //combinations of LO DB=1 Wilson coefficients
    gslpp::complex K_1; //3*C_1^2 + 2*C_1 * C_2
    gslpp::complex K_2; //C_2^2
};

/**
 * @}
 */

#endif /* AMPDB2_H */
