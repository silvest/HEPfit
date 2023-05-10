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
    AmpDB2(const StandardModel& SM_i);

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
    
    /**
    * @brief The value of @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bd}@f$
    */
    gslpp::complex getGamma21overM21_Bd(orders order){
        return Gamma21overM21_Bd(order);
    }

    /**
    * @brief The value of @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$
    */
    gslpp::complex getGamma21overM21_Bs(orders order){
        return Gamma21overM21_Bs(order);
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
    gslpp::complex Gamma21overM21_Bd(orders order);
    

    /**
    * @brief A method to compute @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$.
    * @param[in] order the %QCD order of the computation
    * @return @f$\frac{\Gamma_{21},M_{21}}^{bs}@f$
    */
    gslpp::complex Gamma21overM21_Bs(orders order);
    
    
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
    double mu_1;
    double mu_2;    
    gslpp::vector<gslpp::complex> c(quark q); //requires computeCKMandMasses(); before use
    gslpp::complex delta_1overm(quark q); //requires computeCKMandMasses(); before use
    
//resummation to use z_bar instead of z and  eliminate z ln z terms (hep-ph/0612167)
    bool flag_resumz = true;
    
//access calculated function values
    double F0(quarks qq, int k, int i, int j);
    double F1(quarks qq, int k, int i, int j);
    double F(quarks qq, int k, int i, int j);
    double P(quarks qq, int k, int i, int j);
    gslpp::complex D(quarks qq, int k);
    gslpp::complex deltas_1overm(quarks qq, quark q);  //require computeCKMandMasses
    gslpp::vector<double> me = gslpp::vector<double>(5, 0.);
    double me1tilde;
    gslpp::vector<double> me_R = gslpp::vector<double>(4, 0.);
    
//calculate function values
    void computeF0();
    void computeF1();   
    void computeP();
    void computeD(); //requires F and P
    void compute_deltas_1overm(quark q); //require Wilson and computeCKMandMasses
    void compute_matrixelements(quark q); //require computeCKMandMasses

//array for caching function values
    double cacheF0[24];
    double cacheF1[24];
    double cacheP[84];
    gslpp::complex cacheD[6];
    gslpp::complex cache_deltas_1overm[6];

//returns position in the corresponding array
    int indexF(quarks qq, int k, int i, int j);
    int indexP(quarks qq, int k, int i, int j);
    int indexD(quarks qq, int k);
    int index_deltas(quarks qq, quark q);

    //CKM elements
    void computeCKMandMasses(orders order);

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
    double as_4pi; //alpha_s/(4Pi)
    
    double Md;
    double Ms;
    double Mc;
    double Mb;
    double MB;
    double MB_s;
    double Mb2;
    double MB2;
    
    //Buras basis pdf/hep-ph/9512380v1
    void computeWilsonCoeffs();
    double lastInput_computeWilsonCoeffs = NAN;
    void computeWilsonCoeffsDB1bsg();
    double lastInput_computeWilsonCoeffsDB1bsg = NAN;
    gslpp::complex cacheC[6];
    gslpp::complex C_8G;
    gslpp::complex C(int i);
    
    //combinations of Wilson coeffients arxiv.org/abs/hep-ph/0202010
    gslpp::complex K_1;
    gslpp::complex K_2;
    
    
    //NNLO
    void computeWilsonCoeffsMisiak();
    double lastInput_computeWilsonCoeffsMisiak = NAN;
    void compute_pp_s();
    int index_p(quarks qq, int i, int j, int n);
    double cache_p[576];
    double cache_ps[576];
    gslpp::complex H();
    gslpp::complex H_s();
    gslpp::complex H(quarks qq);
    gslpp::complex H_s(quarks qq);
    double p(quarks qq, int i, int j);
    double p_s(quarks qq, int i, int j);
    double p(quarks qq, int i, int j, int n);
    double p_s(quarks qq, int i, int j, int n);
    gslpp::complex lambda_c;
    gslpp::complex lambda_u;
    //for checking cache;
    double lastInput_compute_pp_s[3] = {NAN, NAN, NAN};
    
    const double M_PI4 = M_PI2 * M_PI2;
    bool orderofp[3] = {true, true, true};
        
    //RI
    //arXiv:hep-ph/0606197v1, scheme dependent different from Mathematica
    gslpp::matrix<double> meMStoRI;
    gslpp::matrix<double> coeffsMStoRI;
};

/**
 * @}
 */

#endif	/* AMPDB2_H */

