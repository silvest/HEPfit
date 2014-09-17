/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BKSTARLL_H
#define	BKSTARLL_H

#include <math.h>
#include "Flavour.h"
#include <StandardModel.h>
#include <ThObservable.h>
#include <gsl/gsl_integration.h>
#include <assert.h>


#define CUTOFF 10    //cutoff between LCSR and lattice values for Form Factors, in GeV^2

/*******************************************************************************
 * GSL Function Conversion BEGIN                                                  *
 * ****************************************************************************/

// Option 1. To be called with:
//gsl_function_pp Fp( boost::bind(&Class::member_function, &class, _1) );
//gsl_function *F = static_cast<gsl_function*>(&Fp);

class gsl_function_pp : public gsl_function
{
public:
    gsl_function_pp(std::function<double(double)> const& func) : _func(func){
        function=&gsl_function_pp::invoke;
        params=this;
    }
private:
    std::function<double(double)> _func;
    static double invoke(double x, void *params) {
        return static_cast<gsl_function_pp*>(params)->_func(x);
    }
};


//Option 2. To be used with:
//gslFunction gslF = convertToGslFunction( boost::bind( &Class::member_function, &class, _1 ) );

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

//Option 3. To be used with:
//Class* ptr2 = class;
//auto ptr = [=](double x)->double{return ptr2->member_function(x);};
//gsl_function_p<decltype(ptr)> Fp(ptr);
//gsl_function *F = static_cast<gsl_function*>(&Fp);

template< typename F >
class gsl_function_p : public gsl_function {
public:
    gsl_function_p(const F& func) : _func(func) {
        function = &gsl_function_p::invoke;
        params=this;
    }
private:
    const F& _func;
    static double invoke(double x, void *params) {
        return static_cast<gsl_function_p*>(params)->_func(x);
    }
};

/*******************************************************************************
 * GSL Function conversion END                                                     *
 * ****************************************************************************/

/**
 * @class BKstarll
 * @ingroup flavour
 * @brief A class for the decay B -> K^*ll. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class BKstarll : public ThObservable {
public:
    BKstarll(const StandardModel& SM_i, StandardModel::lepton lep_i = StandardModel::MU);
    virtual ~BKstarll();
    void updateParameters();
    virtual double computeThValue()=0;
    
    double GF;            //Fermi constant
    double ale;           //alpha electromagnetic
    double Mm;            //muon mass
    double MB;            //B meson mass
    double MKstar;        //K star meson mass
    double Mb;            //b quark mass
    double Ms;            //s quark mass
    double MW;            //W boson mass
    complex lambda_t;     //Vckm factor
    double b;             //BF of the decay K^* -> K pi
    double h[3];          //parameter that contains the contribution from the hadronic hamiltonian  
    double q2;            //q^2 of the decay
    double qmin, qmax;
    
    /*lattice fit parameters*/
    double a_0V, a_1V, c_01V, c_01sV, dmV;
    double a_0A0, a_1A0, c_01A0, c_01sA0, dmA0;
    double a_0A1, a_1A1, c_01A1, c_01sA1, dmA1;
    double a_0A12, a_1A12, c_01A12, c_01sA12, dmA12;
    double a_0T1, a_1T1, c_01T1, c_01sT1, dmT1;
    double a_0T2, a_1T2, c_01T2, c_01sT2, dmT2;
    double a_0T23, a_1T23, c_01T23, c_01sT23, dmT23;
    
    /*LCSR fit parameters*/
    double r_1V, r_2V, m_RV, m_fit2V;
    double r_1A0, r_2A0, m_RA0, m_fit2A0;
    double r_2A1, m_fit2A1;
    double r_1A2, r_2A2, m_fit2A2;
    double r_1T1, r_2T1, m_RT1, m_fit2T1;
    double r_2T2, m_fit2T2;
    double r_1T3t, r_2T3t, m_fit2T3t;

    vector<complex> ** allcoeff;
    vector<complex> ** allcoeffprime;
    
    
    
    /**
    * @brief \f$ LCSR_fit1 \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] r_1 fit parameter
    * @param[in] r_2 fit parameter
    * @param[in] m_R2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return return the first fit function from arXiv:hep-ph/0412079v1
    */
    double LCSR_fit1(double q2, double r_1, double r_2, double m_R2, double m_fit2);
    
    
    /**
    * @brief \f$ LCSR_fit2 \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] r_1 fit parameter
    * @param[in] r_2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return return the second fit function from arXiv:hep-ph/0412079v1
    */
    double LCSR_fit2(double q2, double r_1, double r_2, double m_fit2);
    
    
    /**
    * @brief \f$ LCSR_fit3 \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] r_2 fit parameter
    * @param[in] m_fit2 fit parameter
    * @return return the third fit function from arXiv:hep-ph/0412079v1
    */
    double LCSR_fit3(double q2, double r_2, double m_fit2);
    
    
    /**
    * @brief \f$ z \f$
    * @param[in] q2 q^2 of the decay
    * @return return the lattice parameter z from arXiv:1310.3722v3
    */
    double z(double q2);
    
    
    /**
    * @brief \f$ lat_fit(q^2) \f$
    * @param[in] q2 q^2 of the decay
    * @param[in] a_0 fit parameter
    * @param[in] a_1 fit parameter
    * @param[in] c_01 fit parameter
    * @param[in] c_01s fit parameter
    * @param[in] q2 q^2 of the decay
    * @return return the lattice pole factor P(q^2,dm) from arXiv:1310.3722v3
    */
    double lat_fit(double q2, double a_0, double a_1, double c_01, double c_01s, double dm);
    
    
    /**
    * @brief \f$ V \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor V(q^2)
    */
    double V(double q2);

    
    /**
    * @brief \f$ A_0 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor A_0(q^2)
    */
    double A_0(double q2);

    
    /**
    * @brief \f$ A_1 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor A_1(q^2)
    */
    double A_1(double q2);

    
    /**
    * @brief \f$ A_2 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor A_2(q^2)
    */
    double A_2(double q2);

    
    /**
    * @brief \f$ T_1 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor T_1(q^2)
    */
    double T_1(double q2);

    
    /**
    * @brief \f$ V \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor V(q^2)
    */
    double T_2(double q2);

    
    /**
    * @brief \f$ T_3tilde \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor T_3tilde(q^2)
    */
    double T_3tilde(double q2);

    
    /**
    * @brief \f$ T_3 \f$
    * @param[in] q2 q^2 of the decay
    * @return return the transverse form factor T_3(q^2)
    */
    double T_3(double q2);
    
    
    /**
    * polarization   i
    *      0         0
    *      +         1
    *      -         2
    */

    
    /**
    * @brief \f$ V_L \f$
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor V_L(lambda)
    */
    double V_L(int i, double q2);

    
    /**
    * @brief \f$ V_R \f$
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor V_R(lambda)
    */
    double V_R(int i, double q2);


    /**
    * @brief \f$ T_L \f$
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor T_L(lambda)
    */
    double T_L(int i, double q2);


    /**
    * @brief \f$ T_R \f$ 
    * @param[in] i polarization
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor T_R(lambda)
    */
    double T_R(int i, double q2);


    /**
    * 
    * @brief \f$ S_L \f$
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor S_L
    */
    double S_L(double q2);


    /**
    * 
    * @brief \f$ S_R \f$
    * @param[in] q2 q^2 of the decay
    * @return return the helicity form factor S_R
    */
    double S_R(double q2);
    
    
    /**
    * @brief \f$ N \f$ 
    * @return return the helicity amplitude normalization factor N
    */
    
    complex N();
    
    
    /**
    * @brief \f$ H_V \f$ 
    * @param[in] i polarization lambda
    * @param[in] q2 q^2 of the decay
    * @return return the helicity amplitude H_V(lambda)
    */
    gslpp::complex H_V(int i, double q2);


    /**
    * @brief \f$ H_A \f$ 
    * @param[in] i polarization lambda
    * @param[in] q2 q^2 of the decay
    * @return return the helicity amplitude H_A(lambda)
    */
    gslpp::complex H_A(int i, double q2);


    /**
    * @brief \f$ H_S \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the helicity amplitude H_S
    */
    gslpp::complex H_S(double q2);


    /**
    * @brief \f$ H_P \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the helicity amplitude H_P
    */
    gslpp::complex H_P(double q2);
    
    
    /**
    * @brief \f$ k^2 \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the square of the 3-momentum of the recoiling meson in the B rest frame
    */
    double k2 (double q2);
    
    
    /**
    * @brief \f$ beta \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the factor beta used in the angular coefficients I_i
    */
    double beta (double q2);
    
    
    /**
    * @brief \f$ lambda \f$ 
    * @param[in] q2 q^2 of the decay
    * @return return the factor lambda used in the angular coefficients I_i
    */
    double lambda(double q2);

    
    /**
    * @brief \f$ F \f$ 
    * @param[in] q2 q^2 of the decay
    * @param[in] b BF of the decay K* -> K pi
    * @return return the factor F used in the angular coefficients I_i
    */
    double F(double q2, double b);
    
    
    /**
    * i values:
    * 0 = 1c
    * 1 = 1s
    * 2 = 2c
    * 3 = 2s
    * 4 = 3
    * 5 = 4
    * 6 = 5
    * 7 = 6s
    * 8 = 6c
    * 9 = 7
    * 10 = 8
    * 11 = 9
    */
    
    
    /**
    * @brief \f$ I_{i} \f$ 
    * @param[in] i index of the angular coefficient
    * @param[in] q2 q^2 of the decay
    * @return return the angular coefficient I_i
    */
    double  I(int i, double q2);
    
    
    /**
    * @brief \f$ I_{i} \f$ 
    * @param[in] i index of the angular coefficient
    * @param[in] q2 q^2 of the decay
    * @return return the angular coefficient \bar{I}_i
    */
    double I_bar(int i, double q2);
    
    
    /**
    * @brief \f$ Sigma_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q2 q^2 of the decay
    * @return return the CP average Sigma_i
    */
    double Sigma(int i, double q2);
    
    
    /**
    * @brief \f$ Sigma_{i} \f$ 
    * @param[in] i index of the angular coefficient I_i
    * @param[in] q2 q^2 of the decay
    * @return return the CP asymmetry Delta_i
    */
    double Delta(int i, double q2);
    
    
    double getSigma3(double q2){
        return Sigma(3, q2);
    };
    
    double getSigma4(double q2){
        return Sigma(4, q2);
    };


private:
    const StandardModel& mySM;
    StandardModel::lepton lep;

};


/**
 * @class P_1
 * @ingroup flavour
 * @brief A class for the clean observable P_1. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_1 : public BKstarll{
public:
    
    /**
    * @brief \f$ P_{1} \f$ 
    */
    P_1(const StandardModel& SM_i, StandardModel::lepton lep_i = StandardModel::MU);
    
    
    /**
    * @return return the clean observable P_1
    */
    double computeThValue ();
    
private:
    gsl_function F1, F2;
};


/**
 * @class P_2
 * @ingroup flavour
 * @brief A class for the clean observable P_2. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_2 : public BKstarll {
public:
    
    /**
    * @brief \f$ P_{2} \f$ 
    */
    P_2(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable P_2
    */
    double computeThValue ();
};


/**
 * @class P_3
 * @ingroup flavour
 * @brief A class for the clean observable P_3. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_3 : public BKstarll{
public:
    
    /**
    * @brief \f$ P_{3} \f$ 
    */
    P_3(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable P_3
    */
    double computeThValue ();
};


/**
 * @class P_4Prime
 * @ingroup flavour
 * @brief A class for the clean observable P'_4. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_4Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ P'_{4} \f$ 
    */
    P_4Prime(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable P'_4
    */
    double computeThValue ();
};


/**
 * @class P_5Prime
 * @ingroup flavour
 * @brief A class for the clean observable P'_5. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_5Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ P'_{5} \f$ 
    */

    P_5Prime(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable P'_5
    */
    double computeThValue ();
};


/**
 * @class P_6Prime
 * @ingroup flavour
 * @brief A class for the clean observable P'_6. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P_6Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ P'_{6} \f$ 
    */
    P_6Prime(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);

    /**
    * @return return the clean observable P'_6
    */
    double computeThValue ();
};


/**
 * @class GammaPrime
 * @ingroup flavour
 * @brief A class for the clean observable Gamma'. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class GammaPrime : public BKstarll{
public:
    
    /**
    * @brief \f$ Gamma' \f$ 
    */
    GammaPrime(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable Gamma'
    */
    double computeThValue ();
};


/**
 * @class ACP
 * @ingroup flavour
 * @brief A class for the clean observable Gamma'. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class ACP : public BKstarll{
public:
    
    /**
    * @brief \f$ A_{CP} \f$ 
    */
    ACP(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable ACP
    */
    double computeThValue ();
    
    const StandardModel& mySM;
    StandardModel::lepton lep;
};


/**
 * @class P3CP
 * @ingroup flavour
 * @brief A class for the clean observable Gamma'. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class P3CP : public BKstarll{
public:
    
    /**
    * @brief \f$ P_3^{CP} \f$ 
    */
    P3CP(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable P3CP
    */
    double computeThValue ();
};


/**
 * @class F_L
 * @ingroup flavour
 * @brief A class for the clean observable F_L. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class F_L : public BKstarll{
public:
    
    /**
    * @brief \f$ F_L \f$ 
    */
    F_L(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);

    
    /**
    * @return return the clean observable F_L
    */
    double computeThValue ();
    
    const StandardModel& mySM;
    StandardModel::lepton lep;
};


/**
 * @class M_1Prime
 * @ingroup flavour
 * @brief A class for the clean observable M'_1. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class M_1Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ M'_1 \f$ 
    */
    M_1Prime(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable M'_1
    */
    double computeThValue ();
};


/**
 * @class M_2Prime
 * @ingroup flavour
 * @brief A class for the clean observable M'_2. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class M_2Prime : public BKstarll{
public:
    
    /**
    * @brief \f$ M'_2 \f$ 
    */
    M_2Prime(const StandardModel& SM_i, double q2, StandardModel::lepton lep_i = StandardModel::MU);
    
    /**
    * @return return the clean observable M'_2
    */
    double computeThValue ();
};


#endif	/* BKSTARLL_H */

    
