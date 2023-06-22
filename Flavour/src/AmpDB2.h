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
    enum mass_schemes {pole, MSbar, PS};
    
    /**
    * @brief The value of @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$
    */
    gslpp::complex getGamma21overM21_Bd(orders order, mass_schemes mass_scheme = MSbar){
        return Gamma21overM21_Bd(order, mass_scheme);
    }

    /**
    * @brief The value of @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$
    */
    gslpp::complex getGamma21overM21_Bs(orders order, mass_schemes mass_scheme = MSbar){
        return Gamma21overM21_Bs(order, mass_scheme);
    }
    
    gslpp::complex Gamma21overM21_BdFULLNLO1();
    gslpp::complex Gamma21overM21_BsFULLNLO1();
    

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
    * @brief A method to compute @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$
    */
    gslpp::complex Gamma21overM21_Bd(orders order, mass_schemes mass_scheme);
    

    /**
    * @brief A method to compute @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$.
    * @param[in] order the %QCD order of the computation
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
    enum quark {d,s};
    enum quarks {cc, cu, uu};
    
private:
    double zeta3 = gslpp_special_functions::zeta(3);
    double log2 = log(2);
    double log3 = log(3);
    double sqrt3 = sqrt(3);        
    double sqrt5 = sqrt(5);
    double log12sqrt52 = log(0.5 + sqrt5/2.);    
    //equations (3.102, 6.23)
    double t_2 = -0.389012;
    double Cl2PI3 = 1.014941;
    
    double mu_1;
    double mu_1_overm;
    double mu_2;    
    gslpp::vector<gslpp::complex> c(quark q); //requires computeCKMandMasses(); before use
    gslpp::complex delta_1overm(quark q); //requires computeCKMandMasses(); before use
    gslpp::complex delta_1overm_NLO1(quark q); //requires computeCKMandMasses(); before use
    
//resummation to use z_bar instead of z and  eliminate z ln z terms (hep-ph/0612167)
    bool flag_resumz;
    
//access calculated function values
    double F0(quarks qq, int k, int i, int j);
    double F1(quarks qq, int k, int i, int j);
    double F(quarks qq, int k, int i, int j);
    double P(quarks qq, int k, int i, int j);
    gslpp::complex D(quarks qq, int k);
    gslpp::complex deltas_1overm_NLO1(quarks qq, quark q);  //require computeCKMandMasses
    gslpp::vector<double> me = gslpp::vector<double>(5, 0.);
    gslpp::vector<double> me_R = gslpp::vector<double>(5, 0.); //R0 to R4
    gslpp::vector<double> me_Rtilde = gslpp::vector<double>(3, 0.);
    
//calculate function values
    void computeF0();
    void computeF1();   
    void computeP();
    void computeD(); //requires F and P
    void compute_deltas_1overm(quark q); //require Wilson and computeCKMandMasses
    void compute_deltas_1overm_NLO1(quark q); //require Wilson and computeCKMandMasses
    void compute_matrixelements(quark q); //require computeCKMandMasses

//array for caching function values
    double cacheF0[24] = { 0. };
    double cacheF1[24] = { 0. };
    double cacheP[84] = { 0. };
    gslpp::complex cacheD[6] = { 0. };
    gslpp::complex cache_deltas_1overm_NLO1[6] = { 0. };
    gslpp::complex cache_deltas_1overm[6] = { 0. };

//returns position in the corresponding array
    int indexF(quarks qq, int k, int i, int j);
    int indexP(quarks qq, int k, int i, int j);
    int indexD(quarks qq, int k);
    int index_deltas(quarks qq, quark q);

    //CKM elements
    void computeCKMandMasses(orders order, mass_schemes mass_scheme = MSbar);

    gslpp::complex VtbVtd;
    gslpp::complex VtbVts;
    gslpp::complex VtbVtd2;
    gslpp::complex VtbVts2;
    gslpp::complex VcbVcd;
    gslpp::complex VcbVcs;
    gslpp::complex VcbVcd2;
    gslpp::complex VcbVcs2;
    
    //cache for often used values
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
    double z_1overm2;
    double oneminusz_1overm2;
    double sqrt1minus4z_1overm;
    
    double Md;
    double Ms;
    double Mc;
    double Mb;
    double MB;
    double MB_s;
    double Mb2;
    double Mb2_prefactor; //overall Mb^2 prefactor used in the NNLO implementation
    double MB2;
    double Mb_Mb; //MSbar mass of bottom
    double Mb_pole; //pole mass of bottom
    double Mb_PS; //PS mass of bottom 
    //parameter to calculate PS
    double mu_f = 2.;
    double K = 13.44 - 1.04 * 5.;
    double a1 = 10.33 - 1.11 * 5.;
    double b0 = 11. - 2. * 5. / 3.;
    
    //Buras basis pdf/hep-ph/9512380v1
    void computeWilsonCoeffs();
    double lastInput_computeWilsonCoeffs = NAN;
    void computeWilsonCoeffsDB1bsg();
    double lastInput_computeWilsonCoeffsDB1bsg = NAN;
    gslpp::complex cacheC[6] = { 0. };
    gslpp::complex C_8G;
    gslpp::complex C(int i);
    
    //combinations of Wilson coeffients arxiv.org/abs/hep-ph/0202010
    gslpp::complex K_1;
    gslpp::complex K_2;
    
    
    //NNLO
    void computeWilsonCoeffsMisiak();
    double lastInput_computeWilsonCoeffsMisiak = NAN;
    void compute_pp_s();
    void poletoMSbar_pp_s();
    //hep-ph/9912391v2  eq. (11)
    double PoletoMS_as1 = 4./3.;                                
    double PoletoMS_as2 = -(4. * (71./144. + M_PI2/18.) - 3019./288. + 1./6. * zeta3 - M_PI2/9. * log2 - M_PI2/3.);
    void poletoPS_pp_s();   
    double PoletoPS_as1;                
    double PoletoPS_as2;
    
    int index_p(quarks qq, int i, int j, int n);
    double cache_p[576] = { 0. };
    double cache_ps[576] = { 0. };
    gslpp::vector<gslpp::complex> c_H();
    gslpp::complex H(quarks qq);
    gslpp::complex H_s(quarks qq);
    double p(quarks qq, int i, int j);
    double p_s(quarks qq, int i, int j);
    double p(quarks qq, int i, int j, int n);
    double p_s(quarks qq, int i, int j, int n);
    gslpp::complex lambda_c;
    gslpp::complex lambda_u;
    double lastInput_compute_pp_s[3] = {NAN, NAN, NAN};
    
    const double M_PI4 = M_PI2 * M_PI2;
    bool orderofp[3] = {true, true, true};
    
    //RI: hep-lat/0110091
    gslpp::matrix<double> meMStoRI;
    gslpp::matrix<double> coeffsMStoRI;
    bool flag_RI;
    
    //hep-ph/0612167 
    gslpp::complex deltas_1overm(quarks qq, quark q);
    void compute_g();
    gslpp::complex g(quarks qq, int i);
    gslpp::complex gtilde(quarks qq, int i);
    gslpp::complex cacheg[12] = { 0. };
    gslpp::complex cachegtilde[12] = { 0. };
    int indexg(quarks qq, int i);
    gslpp::complex C_1LO;
    gslpp::complex C_2LO;
};

/**
 * @}
 */

#endif	/* AMPDB2_H */

