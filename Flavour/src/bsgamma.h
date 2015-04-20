/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BSGAMMA_H
#define	BSGAMMA_H

#include <ThObservable.h>
#include "Flavour.h"
#include <StandardModel.h>
#include <Polylogarithms.h>
#include <ClausenFunctions.h>

/**
 * @class Bsgamma
 * @ingroup flavour
 * @brief A class for the decay b -> s gamma. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Bsgamma : public ThObservable {
public: 
    /**
    * @brief constructor
    * @param[in] SM_i model
    * @param[in] obsFlag flag to choose which observable to compute
    */
    Bsgamma(const StandardModel& SM_i, int obsFlag);
    
    
    /**
    * @brief \f$ H1 \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function H^(1) obtained integrating the equations from arXiv:hep-ph/0607316v1v2
    */
    //double H1(double E0, double Mu);
    
    
    /**
    * @brief \f$ f \f$
    * @param[in] rho squared ratio of charm mass over bottom mass
    * @return return the function f from arXiv:hep-ph/0611123v1
    */
    //double f(double rho);
    
    
    /**
    * @brief \f$ H2a \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function H^(2,a) obtained integrating the equations from arXiv:hep-ph/0607316v1v2
    */
    //double H2a(double E0, double Mu);
    
    
    /**
    * @brief \f$ H2na \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function H^(2,na) obtained integrating the equations from arXiv:hep-ph/0607316v1v2
    */
    //double H2na(double E0, double Mu);
    
    
    /**
    * @brief \f$ H2NH \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function H^(2,NH) obtained integrating the equations from arXiv:hep-ph/0607316v1v2
    */
    //double H2NH(double E0, double Mu);
    
    
    /**
    * @brief \f$ H2NV \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function H^(2,NV) obtained integrating the equations from arXiv:hep-ph/0611123v1
    */
    //double H2NV(double E0, double Mu);
    
    
    /**
    * @brief \f$ H2NL \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function H^(2,NL) obtained integrating the equations from arXiv:hep-ph/0607316v1v2
    */
    //double H2NL(double E0, double Mu);
    
    
    /**
    * @brief \f$ H2 \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function H^(2) obtained integrating the equations from arXiv:hep-ph/0607316v1v2 and arXiv:hep-ph/0611123v1
    */
    //double H2(double E0, double Mu);
    
    
    /**
    * @brief \f$ G77 \f$
    * @param[in] order perturbation theory order
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function G_{77} obtained integrating the equations from arXiv:hep-ph/0607316v1v2 and arXiv:hep-ph/0611123v1
    */
    //double G77(orders order, double E0, double Mu);
    
    
    /**
    * @brief \f$ dY1 \f$
    * @param[in] E0 energy cutoff
    * @return return the function \delta Y^1 from arXiv:0805.3911
    */
    //double dY1(double E0);
    
    
    /**
    * @brief \f$ Y1 \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(1) from arXiv:0805.3911
    */
    //double Y1(double E0, double Mu);
    
    
    /**
    * @brief \f$ Phi1 \f$
    * @param[in] rho squared ratio of charm mass over bottom mass
    * @return return the function \Phi_1 from arXiv:1005.5587
    */
    //double Phi1(double rho);
    
    
    /**
    * @brief \f$ Phi2 \f$
    * @param[in] rho squared ratio of charm mass over bottom mass
    * @return return the function \Phi_2 from arXiv:1005.5587
    */
    //double Phi2(double rho);
    
    
    /**
    * @brief \f$ Phi3 \f$
    * @param[in] rho squared ratio of charm mass over bottom mass
    * @return return the function \Phi_3 from arXiv:1005.5587
    */
    //double Phi3(double rho);
    
    
    /**
    * @brief \f$ Phi4 \f$
    * @param[in] rho squared ratio of charm mass over bottom mass
    * @return return the function \Phi_4 from arXiv:1005.5587
    */
    //double Phi4(double rho);
    
    
    /**
    * @brief \f$ Y2CF \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,CF) from arXiv:1005.5587
    */
    //double Y2CF(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2CA \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,CA) from arXiv:1005.5587
    */
    //double Y2CA(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2NH \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,NH) from arXiv:0805.3911
    */
    //double Y2NH(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2NV \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,NV) from arXiv:0805.3911
    */
    //double Y2NV(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2NL \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2,NL) from arXiv:0805.3911
    */
    //double Y2NL(double E0, double Mu);
    
    
    /**
    * @brief \f$ Y2 \f$
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function Y^(2) from arXiv:0805.3911 and arXiv:1005.5587
    */
    //double Y2(double E0, double Mu);
    
    
    /**
    * @brief \f$ G78 \f$
    * @param[in] order perturbation theory order
    * @param[in] E0 energy cutoff
    * @param[in] Mu mass scale
    * @return return the function G_{78} from arXiv:0805.3911 and arXiv:1005.5587
    */
    //double G78(orders order, double E0, double Mu);
    
    
    /**
    * @brief \f$ delta \f$
    * @param[in] E0 cutoff energy
    * @return return the function delta(E0)
    */
    double delta(double E0);
    
    
    /**
    * @brief \f$ zeta \f$
    * @param[in] E0 cutoff energy
    * @return return the squared ratio between m_c and m_b^{1s}
    */
    double zeta();
    
    
    /**
    * @brief \f$ a(z) \f$
    * @param[in] z squared ratio between m_c and m_b^{1s}
    * @return return the funcion a(z) as defined in hep-ph/0203135
    */
    gslpp::complex a(double z);
    
    
    /**
    * @brief \f$ b(z) \f$
    * @param[in] z squared ratio between m_c and m_b^{1s}
    * @return return the funcion b(z) as defined in hep-ph/0203135
    */
    gslpp::complex b(double z);
    
    
    /**
    * @brief \f$ r_i(z) \f$
    * @param[in] i function index
    * @param[in] z squared ratio between m_c and m_b^{1s}
    * @return return the funcion r_i(z) as defined in hep-ph/0203135
    */
    gslpp::complex r(int i, double z);
    
    
    /**
    * @brief \f$ Gamma \f$
    * @param[in] t dummy variable to be integrated out
    * @return return the function Gamma(t) from hep-ph/0104034v2
    */
    gslpp::complex Gamma_t(double t);
    
    
    /**
    * @brief \f$ GetPhi221 \f$
    * @param[in] t dummy variable to be integrated out
    * @return return the square of the absolute value of the function Gamma(t)/t + 1/2
    */
    double getPhi221(double t){
        return (Gamma_t(t)/t + 1./2.).abs2();
    };
    
    
    /**
    * @brief \f$ GetPhi221_t \f$
    * @param[in] t dummy variable to be integrated out
    * @return return the square of the absolute value of the function Gamma(t)/t + 1/2 times t
    */
    double getPhi221_t(double t){
        return t*(Gamma_t(t)/t + 1./2.).abs2();
    };
    
    
    /**
    * @brief \f$ GetPhi221_t2 \f$
    * @param[in] t dummy variable to be integrated out
    * @return return the square of the absolute value of the function Gamma(t)/t + 1/2 times t^2
    */
    double getPhi221_t2(double t){
        return t*t*(Gamma_t(t)/t + 1./2.).abs2();
    };
    
    
    /**
    * @brief \f$ GetPhi271 \f$
    * @param[in] t dummy variable to be integrated out
    * @return return the square of the real value of the function Gamma(t) + t/2
    */
    double getPhi271(double t){
        return (Gamma_t(t) + t/2.).real();
    };
    
    
    /**
    * @brief \f$ GetPhi271_t \f$
    * @param[in] t dummy variable to be integrated out
    * @return return the square of the real value of the function Gamma(t) + t/2 times t
    */
    double getPhi271_t(double t){
        return t*(Gamma_t(t) + t/2.).real();
    };
    
    
    /**
    * @brief \f$ Phi_{11}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi22(1) function from  from hep-ph/0104034v2
    */
    double Phi11_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{12}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi22(1) function from  from hep-ph/0104034v2
    */
    double Phi12_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{17}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi27(1) function from  from hep-ph/0104034v2
    */
    double Phi17_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{18}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi27(1) function from  from hep-ph/0104034v2
    */
    double Phi18_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{22}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi22(1) function from  from hep-ph/0104034v2
    */
    double Phi22_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{27}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi27(1) function from  from hep-ph/0104034v2
    */
    double Phi27_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{28}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi27(1) function from  from hep-ph/0104034v2
    */
    double Phi28_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{47}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi42(1) function from  from 1005:1173
    */
    double Phi47_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{77}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi22(1) function from  from hep-ph/0104034v2
    */
    double Phi77_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{78}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi27(1) function from  from hep-ph/0104034v2
    */
    double Phi78_1(double E0);
    
    
    /**
    * @brief \f$ Phi_{88}^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @return return the phi27(1) function from  from hep-ph/0104034v2
    */
    double Phi88_1(double E0);
    
    
    /**
    * @brief \f$ K_{ij}^{(1)} \f$
    * @param[in] i first index
    * @param[in] j second index
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return return the phi27(1) function from  from hep-ph/0104034v2
    */
    double Kij_1(int i, int j, double E0, double mu);
    
    
    /**
    * @brief \f$ computeCoeff \f$
    * @param[in] mu low scale of the decay
    * @return compute the Wilson Coefficient
    */
    void computeCoeff(double mu);
    
    
    /**
    * @brief \f$ P_2^{(1)} \f$
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return return the perturbative part P_2^(1) of the BR
    */
    double P21(double E0, double mu);
    
    
    /**
    * @brief \f$ P_2^{(2)} \f$
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return return the perturbative part P_2^(2) of the BR
    */
    double P22(double E0, double mu);
    
    
    /**
    * @brief \f$ P_3^{(2)} \f$
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return return the perturbative part P_3^(2) of the BR
    */
    double P32(double E0, double mu);
    
    
    /**
    * @brief \f$ P \f$
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @param[in] order perturbation theory order
    * @return return the perturbative part of the BR
    */
    double P(double E0, double mu, orders order);
    
    
    /**
    * @brief \f$ computeBR \f$
    * @param[in] order perturbation theory order
    * @return compute the branching fraction BR
    */
    void computeBR(orders order);
    
    
    /**
    * @brief \f$ computeThValue \f$
    * @return compute the BR and A_{CP}
    */
    double computeThValue();
    
    
private:
    const StandardModel& mySM;
    
    double ale;
    double E0;
    double mu_b;
    double Mb1s;
    double Mc;
    double Ms;
    double BRsl;
    double C;
    gslpp::complex lambda_t;
    gslpp::complex V_cb;
    
    gslpp::vector<gslpp::complex> ** allcoeff;
    
    gslpp::complex C1_0;
    gslpp::complex C2_0;
    gslpp::complex C3_0;
    gslpp::complex C4_0;
    gslpp::complex C5_0;
    gslpp::complex C6_0;
    gslpp::complex C7_0;
    gslpp::complex C8_0;
    
    gslpp::complex C1_1;
    gslpp::complex C2_1;
    gslpp::complex C3_1;
    gslpp::complex C4_1;
    gslpp::complex C5_1;
    gslpp::complex C6_1;
    gslpp::complex C7_1;
    gslpp::complex C8_1;
    
    gslpp::complex C7_2;
    
    int obs;
    
    double BR;
    double BR_conj;
    
    double avaPhi221_1;
    double avaPhi221_2;
    double avaPhi221_3;
    double avaPhi221_4;
    double avaPhi221_5;

    double avaPhi271_1;
    double avaPhi271_2;
    double avaPhi271_3;
    
    double errPhi221_1;
    double errPhi221_2;
    double errPhi221_3;
    double errPhi221_4;
    double errPhi221_5;
    
    double errPhi271_1;
    double errPhi271_2;
    double errPhi271_3;
    
    gsl_function FPhi221_1;
    gsl_function FPhi221_2;
    gsl_function FPhi221_3;
    gsl_function FPhi221_4;
    gsl_function FPhi221_5;
    
    gsl_function FPhi271_1;
    gsl_function FPhi271_2;
    gsl_function FPhi271_3;
    
    gsl_integration_workspace * w_Phi221_1;
    gsl_integration_workspace * w_Phi221_2;
    gsl_integration_workspace * w_Phi221_3;
    gsl_integration_workspace * w_Phi221_4;
    gsl_integration_workspace * w_Phi221_5;
    
    gsl_integration_workspace * w_Phi271_1;
    gsl_integration_workspace * w_Phi271_2;
    gsl_integration_workspace * w_Phi271_3;
    
};

#endif	/* BSGAMMA_H */