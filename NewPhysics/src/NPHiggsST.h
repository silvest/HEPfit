/* 
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPHIGGSST_H
#define	NPHIGGSST_H

#include "NPbase.h"

/**
 * @class NPHiggsST
 * @brief A model class for new physics with a light Higgs-like scalar. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This is a Model class describing new physics models based on the
 * electroweak chiral Lagrangian with a light neutral scalar \f$h\f$. 
 * The effective Lagrangian expansion is truncated at the level of two-derivatives
 * \cite Contino:2010mh,
 * 
 * \f{eqnarray}{
 *   {\cal L} &=&
 *   \frac {1}{2} (\partial_\mu h)^2 - V(h)
 *   + \frac{v^2}{4}\mathrm{Tr}(D_\mu \Sigma^\dagger D^\mu \Sigma)
 *     \left(1+2a\frac{h}{v}+b\frac{h^2}{v^2}+\cdots\right)
 *  \\
 *  &&
 *   - \left[
 *     m_{u,i}\overline{Q}_{L,i}\Sigma \left(\begin{array}{c}u_{R,i}\\0\end{array}\right)\left(1+c_u \frac{h}{v}+\cdots\right)
 *     + m_{d,i}\overline{Q}_{L,i}\Sigma \left(\begin{array}{c}0\\d_{R,i}\end{array}\right)\left(1+c_d \frac{h}{v}+\cdots\right)
 *     + m_{e,i}\overline{L}_{L,i}\Sigma \left(\begin{array}{c}0\\e_{R,i}\end{array}\right)\left(1+c_e \frac{h}{v}+\cdots\right)
 *     + \mathrm{h.c.}
 * \right],
 * \f}
 * where @f$V(h)@f$ is given by 
 * \f[ V(h)
 *   = \frac{m_h}{2}h^2
 *   + \frac{d_3}{6}\left(\frac{3m_h^2}{v}\right)\,h^3
 *   + \frac{d_4}{24}\left(\frac{3m_h^2}{v^2}\right)\,h^4 + \cdots\,.
 * \f]
 *
 * The SM corresponds to the choice @f$a=b=c_u=c_d=c_e=d_3=d_4=1@f$.
 * The dominant deviations from the SM in the electroweak precision observables
 * are induced by the non-standard @f$HVV@f$ coupling @f$a\neq 1@f$. 
 * This generates extra contributions to the oblique parameters
 * @cite Barbieri:2007bh ,
 * @f{eqnarray}{
 * S &=& \frac{1}{12\pi} (1 - a^2)
 *        \ln\bigg(\frac{\Lambda^2}{m_h^2}\bigg)\,,\\
 * T &=& - \frac{3}{16\pi c_W^2} (1 - a^2)
 *        \ln\bigg(\frac{\Lambda^2}{m_h^2}\bigg)\,,
 * @f}
 * where @f$\Lambda@f$ is the cutoff scale of the model.
 *
 *
 * @anchor NPHiggsSTInitialization
 * <h3>Initialization</h3>
 *
 * After creating an instance of the current class with the constructor 
 * NPHiggsST(), it is required to call the initialization method 
 * InitializeModel(). 
 * In a Monte Carlo run, the constructor as well as the initialization
 * method are called in InputParser::ReadParameters().
 *
 * 
 * @anchor NPHiggsSTParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of %NPHiggsST are summarized below:
 * <table class="model">
 * <tr>
 *   <th>Label</th>
 *   <th>LaTeX symbol</th>
 *   <th>Description</th>
 * </tr>
 * <tr>
 *   <td class="mod_name">%a </td>
 *   <td class="mod_symb">\f$a \f$</td>
 *   <td class="mod_desc">The \f$hVV\f$ coupling \f$a\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%b </td>
 *   <td class="mod_symb">\f$b \f$</td>
 *   <td class="mod_desc">The \f$hhVV\f$ coupling \f$b\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%c_u </td>
 *   <td class="mod_symb">\f$c_u \f$</td>
 *   <td class="mod_desc">The \f$h\overline{u}u\f$ coupling \f$c_u\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%c_d</td>
 *   <td class="mod_symb">\f$c_d \f$</td>
 *   <td class="mod_desc">The \f$h\overline{d}d\f$ coupling \f$c_d\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%c_e </td>
 *   <td class="mod_symb">\f$c_e \f$</td>
 *   <td class="mod_desc">The \f$h\overline{e}e\f$ coupling \f$c_e\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%d_3 </td>
 *   <td class="mod_symb">\f$d_3\f$</td>
 *   <td class="mod_desc">The \f$hhh\f$ coupling \f$d_3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%d_4 </td>
 *   <td class="mod_symb">\f$d_4\f$</td>
 *   <td class="mod_desc">The \f$hhhh\f$ coupling \f$d_4\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%LambdaNP</td>
 *   <td class="mod_symb">\f$\Lambda\f$</td>
 *   <td class="mod_desc">The cutoff scale of the model. When this parameter is 
 *   set to 0 in the model configuration file, the cutoff scale is taken to be
 *   @f$\Lambda = 4\pi v/\sqrt{|1-a^2|}@f$, where the @f$W_LW_L@f$ scattering
 *   becomes non-perturbative. </td>
 * </tr>
 * </table>
 *
 *
 * @anchor NPHiggsSTFlags
 * <h3>%Model flags</h3>
 *
 * There is no model flag in the current class.
 *
 *
 * @anchor NPHiggsSTReimplemented
 * <h3>Reimplemented quantities</h3>
 *
 * Compared to the base classes NPbase and StandardModel, the methods for the
 * following quantities are reimplemented in the current class:
 *
 * @li @f$S@f$, @f$T@f$ and @f$U@f$&nbsp;&nbsp;
 * (with obliqueS(), obliqueT() and obliqueU()),
 * @li @f$M_W@f$, @f$c_W^2@f$ and @f$s_W^2@f$&nbsp;&nbsp;
 * (with Mw(), cW2() and sW2()),
 * @li @f$\Gamma_W@f$&nbsp;&nbsp;
 * (with GammaW()). 
 *
 */
class NPHiggsST : public NPbase {
public:
    /**
     * @brief The number of the model parameters in %NPHiggsST.
     */
    static const int NNPHIGGSSTvars = 8;
    /**
     * @brief A string array containing the labels of the model parameters in %NPHiggsST.
     */
    static const std::string NPHIGGSSTvars[NNPHIGGSSTvars];

    /**
     * @brief The default constructor.
     */
    NPHiggsST();

    /**
     * @brief A method to fetch the name of %NPHiggsST.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const 
    {
        return "NPHiggsST";
    }

    /**
     * @brief A method to initialize %NPHiggsST.
     * @return a boolean that is true if model initialization is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief A method to initialize the model parameters.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);

    /**
     * @brief The update method for %NPHiggsST.
     * @details This method updates all the model parameters with giving DPars.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to check if all the mandatory parameters for %NPHiggsST
     * have been provided in the model configuration file.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A method to set a flag of %NPHiggsST.
     * @param[in] name name of a model flag
     * @param[in] value the boolean to be assigned to the flag specified by name
     * @return a boolean that is true if the execution is successful
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     * @brief A method to check the sanity of the set of model flags.
     * @return a boolean that is true if the set of model flags is sane
     */
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief The oblique parameter \f$S\f$.
     * @details
     * @f[ 
     *  S = \frac{1}{12\pi} (1 - a^2) \ln\bigg(\frac{\Lambda^2}{m_h^2}\bigg)\,.
     * @f]
     *
     * See @cite Barbieri:2007bh.
     * @return \f$S\f$
     */
    virtual double obliqueS() const;
        
    /**
     * @brief The oblique parameter \f$T\f$.
     * @details
     * @f[
     *  T = - \frac{3}{16\pi c_W^2} (1 - a^2)
     *        \ln\bigg(\frac{\Lambda^2}{m_h^2}\bigg)\,.
     * @f]
     *
     * See @cite Barbieri:2007bh.
     * @return \f$T\f$
     */
    virtual double obliqueT() const;
    
    /**
     * @brief The oblique parameter \f$U\f$.
     * @return \f$U=0\f$
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The @f$W@f$-boson mass in the on-shell scheme, @f$M_W@f$.
     * @details
     * @f[
     * M_W = M_{W,\mathrm{SM}}
     * \left[
     * 1 - \frac{\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     * \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     * \right].
     * @f]
     *
     * See, e.g., @cite Ciuchini:2013pca.
     * @return @f$M_W@f$ in GeV.
     */
    virtual double Mw() const;

    /**
     * @brief The square of the cosine of the weak mixing angle
     * in the on-shell scheme, denoted as @f$c_W^2@f$.
     * @return @f$c_W^2=\cos^2{\theta_W}=M_W^2/M_Z^2@f$
     */
    virtual double cW2() const;

    /**
     * @brief The square of the sine of the weak mixing angle
     * in the on-shell scheme, denoted as @f$s_W^2@f$.
     * @return @f$s_W^2=\sin^2{\theta_W}=1-M_W^2/M_Z^2@f$
     */
    virtual double sW2() const;

    /**
     * @brief The total width of the W boson, @f$\Gamma_W@f$.
     * @details
     * @f[
     * \Gamma_W = \Gamma_{W,\mathrm{SM}}
     * \left[ 1
     * - \frac{3\alpha(M_Z^2)}{4(c_W^2-s_W^2)}
     *  \left( S - 2c_W^2\,T - \frac{c_W^2-s_W^2}{2s_W^2}\,U \right)
     * \right].
     * @f]
     *
     * See, e.g., @cite Ciuchini:2013pca.
     * @return @f$\Gamma_W@f$ in GeV
     */
    virtual double GammaW() const;
    

    ////////////////////////////////////////////////////////////////////////    
protected:    

    double a;///< The \f$hVV\f$ coupling \f$a\f$.
    double b;///< The \f$hhVV\f$ coupling \f$b\f$.
    double c_u;///< The \f$h\overline{u}u\f$ coupling \f$c_u\f$.
    double c_d;///< The \f$h\overline{d}d\f$ coupling \f$c_d\f$. 
    double c_e;///< The \f$h\overline{e}e\f$ coupling \f$c_e\f$.
    double d_3;///< The \f$hhh\f$ coupling \f$d_3\f$. 
    double d_4;///< The \f$hhhh\f$ coupling \f$d_4\f$. 
    double LambdaNP_in;///< The new physics scale \f$\Lambda\f$

    /**
     * @brief A method to set the value of a parameter of %NPHiggsST.
     * @param[in] name name of a model parameter
     * @param[in] value the value to be assigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPHIGGSST_H */

