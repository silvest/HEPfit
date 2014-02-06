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
 * @brief A class for new physics with a non-standard @f$HVV@f$ couplings. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to describe new physics models based on the electroweak chiral
 * Lagrangian plus all the interactions involving a light neutral scalar \f$H\f$. The
 * effective Lagrangian expansion is truncated at the level of two-derivatives \cite Contino:2010mh,
 * 
 * \f${\cal L}=\frac {1}{2} (\partial_\mu h)^2 -V(h)+\frac{v^2}{4}\mathrm{Tr}(D_\mu \Sigma^\dagger D^mu \Sigma)
 * \left(1+2a\frac{h}{v}+b\frac{h^2}{v^2}+\cdots\right)-\left[
 * m_{u,i}\overline{Q}_{L,i}\Sigma \left(\begin{array}{c}u_{R,i}\\0\end{array}\right)\left(1+c_u \frac{h}{v}+\cdots\right)
 * +m_{d,i}\overline{Q}_{L,i}\Sigma \left(\begin{array}{c}0\\d_{R,i}\end{array}\right)\left(1+c_d \frac{h}{v}+\cdots\right)
 * +m_{\ell,i}\overline{L}_{L,i}\Sigma \left(\begin{array}{c}0\\\ell_{R,i}\end{array}\right)\left(1+c_e \frac{h}{v}+\cdots\right)+h.c.
 * \right],\f$
 * 
 * \f$ V(h)=\frac{m_h}{2}h^2+\frac{d_3}{6}\left(\frac{3m_h^2}{v}\right)h^3+\frac{d_4}{24}\left(\frac{3m_h^2}{v^2}\right)h^4+\cdots\f$
 * 
 * Currently, only the dominant effects in EWPD, coming from non-standard @f$HVV@f$
 * interactions, are considered.
 * 
 * @anchor NPHiggsSTParameters
 * <h3>%Model parameters</h3>
 *
 * The model parameters of NPHiggsST are summarized below: 
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
 *   <td class="mod_name">%cu </td>
 *   <td class="mod_symb">\f$c_u \f$</td>
 *   <td class="mod_desc">The \f$h\overline{u}u\f$ coupling \f$c_u\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%cd</td>
 *   <td class="mod_symb">\f$c_d \f$</td>
 *   <td class="mod_desc">The \f$h\overline{d}d\f$ coupling \f$c_d\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%ce </td>
 *   <td class="mod_symb">\f$c_e \f$</td>
 *   <td class="mod_desc">The \f$h\overline{\ell}\ell\f$ coupling \f$c_e\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%d3 </td>
 *   <td class="mod_symb">\f$d_3\f$</td>
 *   <td class="mod_desc">The \f$hhh\f$ coupling \f$d_3\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%d4 </td>
 *   <td class="mod_symb">\f$d_4\f$</td>
 *   <td class="mod_desc">The \f$hhhh\f$ coupling \f$d_4\f$. </td>
 * </tr>
 * <tr>
 *   <td class="mod_name">%LambdaNP_in</td>
 *   <td class="mod_symb">\f$\Lambda_{NP}\f$</td>
 *   <td class="mod_desc">The new physics scale. </td>
 * </tr>
 * </table>
 * 
 */
class NPHiggsST : public NPbase {
public:
    /**
     * @brief The number of new physics parameters in the model.
     */
    static const int NNPHIGGSSTvars = 8;
    /**
     * @brief A string array with the names of the new physics parameters in the model.
     */
    static const std::string NPHIGGSSTvars[NNPHIGGSSTvars];

    /**
     * @brief Constructor.
     */
    NPHiggsST();

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const 
    {
        return "NPHiggsST";
    }

    /**
     * @brief A method to initialize the model.
     * @return true is model initialization is successful
     */
    virtual bool InitializeModel();
    
    /**
     * @brief A method to initialize the model.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * @return true is model initialization is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);

    /**
     * @brief The update method for the model class.
     * @details This method updates all the parameters of the model every time a
     * new set of parameters is generated.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful.
     */    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * @brief A method to check if all the mandatory parameters for the model have been
     * provided in the model configuration file.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A set method to fix the flags for the model.
     * @param[in] name the name of the flag
     * @param[in] value the value of the flag that can be true or false
     * @return a boolean to designate the success or failure of this procedure
     */
    virtual bool setFlag(const std::string name, const bool value);
    
    /**
     * @brief A method to check the sanity of the set of flags.
     * @return true if the set of flags is sane.
     */
    virtual bool CheckFlags() const;

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief The oblique parameter \f$S\f$.
     * @return the contribution to the oblique parameter \f$S\f$ in this model
     */
    virtual double obliqueS() const;
        
    /**
     * @brief The oblique parameter \f$T\f$.
     * @return the contribution to the oblique parameter \f$T\f$ in this model
     */
    virtual double obliqueT() const;
    
    /**
     * @brief The oblique parameter \f$U\f$.
     * @return the contribution to the oblique parameter \f$U\f$ in this model 
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The \f$W\f$ boson mass.
     * @return the \f$W\f$-boson mass in GeV
     */
    virtual double Mw() const;

    /**
     * @brief The (square of the) cosine of the weak angle \f$\cos^2{\theta_W}\f$.
     * @return the value of \f$\cos^2{\theta_W}\f$ in the On-mass-shell renormalization scheme,
     *  \f$\cos^2{\theta_W}=\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double cW2() const;

    /**
     * @brief The (square of the) sine of the weak angle \f$\sin^2{\theta_W}\f$.
     * @return the value of \f$\sin^2{\theta_W}\f$ in the On-mass-shell renormalization scheme,
     *  \f$\sin^2{\theta_W}=1-\frac{M_W^2}{M_Z^2}\f$
     */
    virtual double sW2() const;

    /**
     * @brief The \f$W\f$ decay width \f$\Gamma_W\f$.
     * @return the total width of the \f$W\f$ boson in GeV
     */
    virtual double GammaW() const;
    

    ////////////////////////////////////////////////////////////////////////    
protected:    
    double a;///< The \f$hVV\f$ coupling \f$a\f$.
    double b;///< The \f$hhVV\f$ coupling \f$b\f$.
    double c_u;///< The \f$h\overline{u}u\f$ coupling \f$c_u\f$.
    double c_d;///< The \f$h\overline{d}d\f$ coupling \f$c_d\f$. 
    double c_e;///< The \f$h\overline{\ell}\ell\f$ coupling \f$c_e\f$. 
    double d_3;///< The \f$hhh\f$ coupling \f$d_3\f$. 
    double d_4;///< The \f$hhhh\f$ coupling \f$d_4\f$. 
    double LambdaNP_in;///< The new physics scale \f$\Lambda_{NP}\f$
     /**
     * @brief A set method to fix the parameters of the model.
     * @param[in] name a string with the parameter name
     * @param[in] value the value to be asigned to the parameter specified by name
     */
    virtual void setParameter(const std::string name, const double& value);

};

#endif	/* NPHIGGSST_H */

