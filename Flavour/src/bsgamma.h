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
 * @ingroup Flavour
 * @brief A class for the @f$b \to s \gamma@f$ decay. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to 
 * compute the observables relative to the @f$b \to s \gamma@f$ decay.
 */
class Bsgamma : public ThObservable {
public: 
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    * @param[in] obsFlag flag to choose which observable to compute
    */
    Bsgamma(const StandardModel& SM_i, int obsFlag);
    
    
    /**
    * @brief The cutoff energy function \f$ \delta \f$.
    * @param[in] E0 cutoff energy
    * @return \f$ \delta(E0) \f$ 
    */
    double delta(double E0);
    
    
    /**
    * @brief The squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$, \f$ z \f$.
    * @return \f$ z \f$
    */
    double zeta();
    
    
    /**
    * @brief The funcion \f$ a(z) \f$ as defined in hep-ph/0203135.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ a(z) \f$ 
    */
    gslpp::complex a(double z);
    
    
    /**
    * @brief The funcion \f$ b(z) \f$ as defined in hep-ph/0203135.
    * @param[in] z squared ratio between \f$m_c\f$ and \f$m_b^{1s}\f$
    * @return \f$ b(z) \f$ 
    */
    gslpp::complex b(double z);
    
    
    /**
    * @brief The funcion \f$ r_i(z) \f$ as defined in hep-ph/0203135.
    * @param[in] i function index
    * @param[in] z squared ratio between m_c and m_b^{1s}
    * @return \f$ r_i(z) \f$
    */
    gslpp::complex r(int i, double z);
    
    
    /**
    * @brief The function \f$ \Gamma \f$ as defined in hep-ph/0104034v2.
    * @param[in] t dummy variable to be integrated out
    * @return \f$ \Gamma \f$ 
    */
    gslpp::complex Gamma_t(double t);
    
    
    /**
    * @brief The square of the absolute value of the function \f$\frac{\Gamma(t)}{t} + \frac{1}{2}\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$|\frac{\Gamma(t)}{t} + \frac{1}{2}|^2\f$ 
    */
    double getPhi221(double t){
        return (Gamma_t(t)/t + 1./2.).abs2();
    };
    
    
    /**
    * @brief The square of the absolute value of the function \f$\frac{\Gamma(t)}{t} + \frac{1}{2}\f$ times \f$t\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$t|\frac{\Gamma(t)}{t} + \frac{1}{2}|^2\f$ 
    */
    double getPhi221_t(double t){
        return t*(Gamma_t(t)/t + 1./2.).abs2();
    };
    
    
    /**
    * @brief The square of the absolute value of the function \f$\frac{\Gamma(t)}{t} + \frac{1}{2}\f$ times \f$t^2\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$t^2|\frac{\Gamma(t)}{t} + \frac{1}{2}|^2\f$ 
    */
    double getPhi221_t2(double t){
        return t*t*(Gamma_t(t)/t + 1./2.).abs2();
    };
    
    
    /**
    * @brief The square of the real part of the function \f$\Gamma(t) + \frac{t}{2}\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$\mathrm{RE}({\Gamma(t) + \frac{t}{2})\f$ 
    */
    double getPhi271(double t){
        return (Gamma_t(t) + t/2.).real();
    };
    
    
    /**
    * @brief The square of the real part of the function \f$\Gamma(t) + \frac{t}{2}\f$ times \f$t\f$.
    * @param[in] t dummy variable to be integrated out
    * @return \f$t\mathrm{RE}({\Gamma(t) + \frac{t}{2})\f$ 
    */
    double getPhi271_t(double t){
        return t*(Gamma_t(t) + t/2.).real();
    };
    
    
    /**
    * @brief The \f$ \Phi_{11}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{11}^{(1)} \f$
    */
    double Phi11_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{12}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{12}^{(1)} \f$
    */
    double Phi12_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{17}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{17}^{(1)} \f$
    */
    double Phi17_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{18}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{18}^{(1)} \f$
    */
    double Phi18_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{22}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{22}^{(1)} \f$
    */
    double Phi22_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{27}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{27}^{(1)} \f$
    */
    double Phi27_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{28}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{28}^{(1)} \f$
    */
    double Phi28_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{47}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{47}^{(1)} \f$
    */
    double Phi47_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{77}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{77}^{(1)} \f$
    */
    double Phi77_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{78}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{78}^{(1)} \f$
    */
    double Phi78_1(double E0);
    
    
    /**
    * @brief The \f$ \Phi_{88}^{(1)} \f$ function from hep-ph/0104034v2.
    * @param[in] E0 energy cutoff
    * @return \f$ \Phi_{88}^{(1)} \f$
    */
    double Phi88_1(double E0);
    
    
    /**
    * @brief The \f$ K_{ij}^{(1)} \f$ function from arXiv:1005.1173.
    * @param[in] i first index
    * @param[in] j second index
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ K_{ij}^{(1)} \f$
    */
    double Kij_1(int i, int j, double E0, double mu);
    
    
    /**
    * @brief Compute the Wilson Coefficient.
    * @param[in] mu low scale of the decay
    */
    void computeCoeff(double mu);
    
    
    /**
    * @brief The perturbative part \f$ P_2^{(1)} \f$ of the BR as defined in arXiv:1005.1173.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ P_2^{(1)} \f$
    */
    double P21(double E0, double mu);
    
    
    /**
    * @brief The perturbative part \f$ P_2^{(2)} \f$ of the BR as defined in arXiv:1005.1173.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ P_2^{(2)} \f$
    */
    double P22(double E0, double mu);
    
    
    /**
    * @brief The perturbative part \f$ P_3^{(2)} \f$ of the BR as defined in arXiv:1005.1173.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @return \f$ P_3^{(2)} \f$
    */
    double P32(double E0, double mu);
    
    
    /**
    * @brief The perturbative part of the \f$BR\f$ as defined in arXiv:1005.1173, \f$P\f$.
    * @param[in] E0 energy cutoff
    * @param[in] mu low scale of the decay
    * @param[in] order perturbation theory order
    * @return \f$P\f$
    */
    double P(double E0, double mu, orders order);
    
    
    /**
    * @brief The \f$BR\f$ as defined in arXiv:1005.1173.
    * @param[in] order perturbation theory order
    */
    void computeBR(orders order);
    
    
    /**
    * @brief The \f$BR\f$ as defined in arXiv:1005.1173.
    * @return \f$BR\f$
    */
    double computeThValue();
    
    
private:
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