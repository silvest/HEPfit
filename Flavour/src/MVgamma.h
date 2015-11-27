/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MVGAMMA_H
#define	MVGAMMA_H

#include <StandardModel.h>
#include <ThObservable.h>



/**
 * @class MVgamma
 * @ingroup Flavour
 * @brief A class for the @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute all the functions needed in order to 
 * compute the observables relative to the @f$M \to V \gamma@f$ decay. After the
 * parameters are updated in updateParameters() and the form factor @f$ T_1 @f$
 * is computed in T_1() following @cite Straub:2015ica, the QCDF corrections to the Wilson coefficient @f$ C_7 @f$
 * is computed in the functions G1(), G8(), H1() and H8() @cite Bosch:2001gv. The helicity amplitudes
 * @f$H_V^{(+,-)},\overline{H}_V^{(+,-)}@f$ are build in H_V_p(), H_V_m(), H_V_p_bar() 
 * and H_V_m_bar() following @cite Jager:2012uw, in order to be further used to build the observables.
 */
class MVgamma  : public ThObservable {
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~MVgamma();
    
    /**
     * @brief The update parameter method for MVgamma.
     */
    void updateParameters();
    
    double GF;            /**<Fermi constant */
    double ale;           /**<alpha electromagnetic */
    double MM;            /**<initial meson mass */
    double MM2;           /**<square of the initial meson mass */
    double MV;            /**<final vector meson mass */
    double Mb;            /**<b quark mass */
    double Mc;            /**<c quark mass */
    double mu_b;          /**<b mass scale */
    double mu_h;          /**<sqrt(mu_b*lambda_QCD) */
    double width;         /**<initial meson width */
    double Ms;            /**<s quark mass */
    double MW;            /**<W boson mass */
    gslpp::complex lambda_t;     /**<Vckm factor */
    gslpp::complex h[2];         /**<parameter that contains the contribution from the hadronic hamiltonian */
    double lambda;        /**<cinematic parameter */
    
    double a_0T1;/**<LCSR fit parameter */
    double a_1T1;/**<LCSR fit parameter */
    double a_2T1;/**<LCSR fit parameter */
    double MRT1_2;/**<LCSR fit parameter */
    
    gslpp::vector<gslpp::complex> ** allcoeff;/**<vector that contains the Wilson coeffients at mub*/
    gslpp::vector<gslpp::complex> ** allcoeffh;/**<vector that contains the Wilson coeffients at muh*/
    gslpp::vector<gslpp::complex> ** allcoeffprime;/**<vector that contains the primed Wilson coeffients at mub*/
    
    gslpp::complex C_7;/**<Wilson coeffients @f$C_7@f$*/
    gslpp::complex C_7p;/**<Wilson coeffients @f$C_7'@f$*/
    gslpp::complex C_2;/**<Wilson coeffients @f$C_2(mu_b)@f$*/
    gslpp::complex C_8;/**<Wilson coeffients @f$C_8(mu_b)@f$*/
    gslpp::complex C_2h;/**<Wilson coeffients @f$C_2(mu_h)@f$*/
    gslpp::complex C_8h;/**<Wilson coeffients @f$C_8(mu_h)@f$*/
    
    
    /**
    * @brief The transverse form factor @f$ T_1 @f$.
    * @return @f$ T_1 @f$ 
    */
    double T_1();

    /**
    * @brief The helicity amplitude @f$ H_V^+ @f$.
    * @return @f$ H_V^+ @f$ 
    */
    gslpp::complex H_V_p();
    
    /**
    * @brief The helicity amplitude @f$ H_V^- @f$.
    * @return @f$ H_V^- @f$ 
    */
    gslpp::complex H_V_m();
    
    /**
    * @brief The helicity amplitude @f$ \bar{H}_V^+ @f$.
    * @return @f$ \bar{H}_V^+ @f$ 
    */
    gslpp::complex H_V_p_bar();
    
    /**
    * @brief The helicity amplitude @f$ \bar{H}_V^- @f$.
    * @return @f$ \bar{H}_V^- @f$ 
    */
    gslpp::complex H_V_m_bar();

    /**
     * @brief Correction in eq. (42) of @cite Bosch:2001gv.
     * @param s @f$ m_c^2/m_b^2  @f$
     * @return @f$ G_1(s) @f$
     */
    gslpp::complex G1(double s);
    
    /**
     * @brief Correction in eq. (42) of @cite Bosch:2001gv.
     * @return @f$ G_8 @f$
     */
    gslpp::complex G8();

    /**
     * @brief Correction in eq. (42) of @cite Bosch:2001gv.
     * @param s @f$ m_c^2/m_b^2  @f$
     * @return @f$ H_1(s) @f$
     */
    gslpp::complex H1(double s);

    /**
     * @brief Correction in eq. (42) of @cite Bosch:2001gv.
     * @return  @f$ H_8 @f$
     */
    gslpp::complex H8();
    
private:
    StandardModel::meson meson;
    StandardModel::meson vectorM;
};



/**
 * @class BR_MVgamma
 * @ingroup Flavour
 * @brief A class for the @f$BR@f$ in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$BR@f$ in @f$M \to V \gamma@f$ 
 * in terms of the helicity amplitudes @f$H_V^{(+,-)},\overline{H}_V^{(+,-)}@f$, 
 * computed in the MVgamma class:
 * @f[
 * BR = \frac {\alpha_e G_F^2 M_b^2 M_M \lambda}{(4\pi)^2 4 W_M} ( |H_V^+|^2 + |H_A^+|^2 +|\overline{H}_V^-|^2 + |\overline{H}_A^-|^2) \,.
 * @f]
 */
class BR_MVgamma : public MVgamma{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    BR_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i);
    
    /**
    * @brief The @f$BR@f$ in @f$M \to V \gamma@f$.
    * @return @f$BR@f$
    */
    double computeThValue ();

private:
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};



/**
 * @class ACP_MVgamma
 * @ingroup Flavour
 * @brief A class for the @f$A_{CP}@f$ in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$A_{CP}@f$ in @f$M \to V \gamma@f$ 
 * in terms of the helicity amplitudes @f$H_V^{(+,-)},\overline{H}_V^{(+,-)}@f$, 
 * computed in the MVgamma class:
 * @f[
 * A_{CP} = \frac {( |H_V^+|^2 + |H_A^+|^2 - |\overline{H}_V^-|^2 - |\overline{H}_A^-|^2)}{( |H_V^+|^2 + |H_A^+|^2 +|\overline{H}_V^-|^2 + |\overline{H}_A^-|^2)} \,.
 * @f]
 */
class ACP_MVgamma : public MVgamma{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    ACP_MVgamma(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i);
    
    /**
    * @brief The @f$A_{CP}@f$ in @f$M \to V \gamma@f$.
    * @return @f$A_{CP}@f$
    */
    double computeThValue ();

private:
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};



/**
 * @class DC7_1
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^1@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^1@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^1 = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_- - h_+| \,.
 * @f]
 */
class DC7_1 : public MVgamma{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    DC7_1(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^1@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^1@f$
    */
    double computeThValue();

private:
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};



/**
 * @class DC7_2
 * @ingroup Flavour
 * @brief A class for the @f$\Delta C_7^2@f$ correction in @f$M \to V \gamma@f$ decay. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$ 
 * due to the hadronic parameters  @f$h_{+,-}@f$, computed in the MVgamma class:
 * @f[
 * \Delta C_7^2 = \frac {8 \pi^2 M_M^3}{\lambda m_b T_1(0)}|h_- + h_+| \,.
 * @f]
 */
class DC7_2 : public MVgamma{
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     * @param[in] meson_i initial meson of the decay
     * @param[in] vector_i final vector meson of the decay
     */
    DC7_2(const StandardModel& SM_i, StandardModel::meson meson_i, StandardModel::meson vector_i);

    /**
    * @brief The @f$\Delta C_7^2@f$ in @f$M \to V \gamma@f$.
    * @return @f$\Delta C_7^2@f$
    */
    double computeThValue();

private:
    StandardModel::meson meson; /**< Initial meson type. */
    StandardModel::meson vectorM; /**< Final vector meson type. */
};
#endif	/* MVLL_H */

