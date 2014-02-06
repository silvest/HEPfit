/* 
 * Copyright (C) 2013-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPEFFECTIVE_H
#define	NPEFFECTIVE_H

#include <stdexcept>
#include "NPbase.h"

/**
 * @class NPEffective
 * @brief A base class for new physics in the form of contributions to the 
 * dimension-six effective Lagrangian.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class contains the necessary functions to compute new physics 
 * tree-level corrections to electroweak precision observables, in the form of
 * contributions to a general class of SM gauge-invariant dimension-six operators
 * in the effective Lagrangian,
 * 
 * \f${\cal L}_{Eff}={\cal L}_{SM}+\frac{1}{\Lambda_{NP}^2}{\cal L}_6,\f$
 * 
 * with
 * 
 * \f${\cal L}_6=c_{WB} \left(H^\dagger\sigma_aH\right)W_{\mu\nu}^a B^{\mu\nu}+ c_H \left|H^\dagger D_\mu H\right|^2+
 * c_{LL} \frac 12\left(\overline{L}\gamma^\mu\sigma_a L\right)\left(\overline{L}\gamma_\mu\sigma_a L\right)+
 * c_{HL}^\prime i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L}\sigma_a\gamma^\mu L\right)+
 * c_{HQ}^\prime i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q}\sigma_a\gamma^\mu Q\right)\f$
 * 
 *\f$\phantom{{\cal L}_6=} +c_{HL} i\left(H^\dagger D_\mu H\right)\left(\overline{L}\gamma^\mu L\right)+
 * c_{HQ} i\left(H^\dagger D_\mu H\right)\left(\overline{Q}\gamma^\mu Q\right)+
 * c_{HE} i\left(H^\dagger D_\mu H\right)\left(\overline{E}\gamma^\mu E\right)+
 * c_{HU} i\left(H^\dagger D_\mu H\right)\left(\overline{U}\gamma^\mu U\right)+
 * c_{HD} i\left(H^\dagger D_\mu H\right)\left(\overline{D}\gamma^\mu D\right)+h.c..\f$
 * 
 * The new physics corrections are parameterized in terms of the operator 
 * contributions to \f$M_W\f$, and to \f$Z\f$-pole observables through the corrections
 * to the different neutral-current effective couplings to leptons and quarks, 
 * \f$\delta g_{L,R}^f\f$ or \f$\delta g_{V,A}^f\f$. This class also contains get
 * methods to retrieve the value of all the included dimension six operators.
 * 
 */
class NPEffective : public NPbase {
public:

    /**
     * @brief Constructor.
     */
    NPEffective();

    /**
     * @brief The name of the model.
     * @return the name of the model as a string
     */
    virtual std::string ModelName() const
    {
        return "NPEffective";
    }

    /**
     * @brief A method to initialize the model.
     * @return true is model initialization is successful
     */
    virtual bool InitializeModel();
    
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
     * @brief A method to initialize the model.
     * @param[in] Dpars a map of parameters that are being updated in the Monte Carlo run
     * @return true is model initialization is successful
     */
    virtual bool Init(const std::map<std::string, double>& DPars);  
   
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
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HL}^\prime\right)_{11}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}^\prime\right)_{11}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L^1}\sigma_a\gamma^\mu L^1\right)\f$
     */
    double getCHL1p() const 
    {
        return cHL1p;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HL}^\prime\right)_{22}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}^\prime\right)_{22}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L^2}\sigma_a\gamma^\mu L^2\right)\f$
     */
    double getCHL2p() const 
    {
        return cHL2p;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HL}^\prime\right)_{33}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}^\prime\right)_{33}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L^3}\sigma_a\gamma^\mu L^3\right) \f$
     */
    double getCHL3p() const 
    {
        return cHL3p;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HL}\right)_{11}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{L^1}\gamma^\mu L^1\right) \f$
     */
    double getCHL1() const 
    {
        return cHL1;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HL}\right)_{22}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{L^2}\gamma^\mu L^2\right) \f$
     */
    double getCHL2() const 
    {
        return cHL2;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HL}\right)_{33}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{L^3}\gamma^\mu L^3\right) \f$
     */
    double getCHL3() const 
    {
        return cHL3;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HQ}^\prime\right)_{11}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}^\prime\right)_{11}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q^1}\sigma_a\gamma^\mu Q^1\right)\f$
     */
    double getCHQ1p() const
    {
        return cHQ1p;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HQ}^\prime\right)_{22}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}^\prime\right)_{22}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q^2}\sigma_a\gamma^\mu Q^2\right)\f$
     */
    double getCHQ2p() const
    {
        return cHQ2p;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HQ}^\prime\right)_{33}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}^\prime\right)_{33}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q^3}\sigma_a\gamma^\mu Q^3\right)\f$
     */
    double getCHQ3p() const
    {
        return cHQ3p;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HQ}\right)_{11}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q^1}\gamma^\mu Q^1\right) \f$
     */
    double getCHQ1() const
    {
        return cHQ1;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HQ}\right)_{22}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q^2}\gamma^\mu Q^2\right) \f$
     */    
    double getCHQ2() const
    {
        return cHQ2;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HQ}\right)_{33}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q^3}\gamma^\mu Q^3\right) \f$
     */
    double getCHQ3() const
    {
        return cHQ3;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HE}\right)_{11}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HE}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{E^1}\gamma^\mu E^1\right) \f$
     */
    double getCHE1() const 
    {
        return cHE1;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HE}\right)_{22}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HE}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{E^2}\gamma^\mu E^2\right) \f$
     */
    double getCHE2() const 
    {
        return cHE2;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HE}\right)_{33}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HE}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{E^3}\gamma^\mu E^3\right) \f$
     */
    double getCHE3() const 
    {
        return cHE3;
    }
  
    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HU}\right)_{11}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HU}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{U^1}\gamma^\mu U^1\right) \f$
     */  
    double getCHU1() const 
    {
        return cHU1;
    }
    
    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HU}\right)_{22}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HU}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{U^2}\gamma^\mu U^2\right) \f$
     */  
    double getCHU2() const 
    {
        return cHU2;
    }
    
    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HU}\right)_{33}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HU}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{U^3}\gamma^\mu U^3\right) \f$
     */  
    double getCHU3() const 
    {
        return cHU3;
    }

    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HD}\right)_{11}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HD}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{D^1}\gamma^\mu D^1\right) \f$
     */
    double getCHD1() const 
    {
        return cHD1;
    }
    
    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HD}\right)_{22}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HD}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{D^2}\gamma^\mu D^2\right) \f$
     */
    double getCHD2() const 
    {
        return cHD2;
    }
    
    /**
     * @brief A get method to retrieve the value of the operator coefficient \f$\left(c_{HD}\right)_{33}\f$.
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HD}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{D^3}\gamma^\mu D^3\right) \f$
     */
    double getCHD3() const 
    {
        return cHD3;
    }

    /**
     * @brief A get method to retrieve the value of the new physics scale \f$\Lambda_\mathrm{NP}\f$.
     * @return the value of the new physics scale \f$\Lambda_\mathrm{NP}\f$ [GeV] 
     */
    double getLambdaNP() const
    {
        return LambdaNP;
    }


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The SM Higgs vev \f$v\f$.
     * @return the SM value of the Higgs vev, as extracted from the experimental value of \f$G_F\f$,
     * \f$v=\frac{1}{\sqrt{\sqrt{2} G_F}}\f$ 
     */
    virtual double v() const;

    /**
     * @brief The SM tree-level \f$W\f$ mass.
     * @return the SM tree-level value of the \f$W\f$ mass [GeV] 
     */
    virtual double Mw_tree() const;

    /**
     * @brief The new physics corrections to the Fermi constant.
     * @return the new physics correction to the Fermi constant, 
     * \f$ G_F\equiv G_F^\mathrm{SM}(1+\Delta G_F)\f$ 
     */
    virtual double DeltaGF() const;


    ////////////////////////////////////////////////////////////////////////     

    /**
     * @brief The oblique parameter \f$S\f$.
     * @return the value of the oblique parameter \f$S\f$
     */
    virtual double obliqueS() const;

    /**
     * @brief The oblique parameter \f$T\f$.
     * @return the value of the oblique parameter \f$T\f$
     */
    virtual double obliqueT() const;

    /**
     * @brief The oblique parameter \f$U\f$.
     * @return The value of the oblique parameter \f$U\f$
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @brief The new physics correction to @f$g_L^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current left-handed coupling @f$g_L^l@f$
     */
    double deltaGLl(StandardModel::lepton l) const;

    /**
     * @brief The new physics correction to @f$g_L^q@f$.
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current left-handed coupling @f$g_L^q@f$
     */
    double deltaGLq(StandardModel::quark q) const;

    /**
     * @brief The new physics correction to @f$g_R^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current right-handed coupling @f$g_R^l@f$
     */    
    double deltaGRl(StandardModel::lepton l) const;
    
    /**
     * @brief The new physics correction to @f$g_R^q@f$.
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current right-handed coupling @f$g_R^q@f$
     */
    double deltaGRq(StandardModel::quark q) const;

    /**
     * @brief The new physics correction to @f$g_V^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current vector coupling @f$g_A^l@f$
     */
    virtual double deltaGVl(StandardModel::lepton l) const;

    /**
     * @brief The new physics correction to @f$g_V^q@f$.
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current vector coupling @f$g_V^q@f$
     */    
    virtual double deltaGVq(StandardModel::quark q) const;

    /**
     * @brief The new physics correction to @f$g_A^l@f$.
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^l@f$
     */    
    virtual double deltaGAl(StandardModel::lepton l) const;

    /**
     * @brief The new physics correction to @f$g_A^q@f$.
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^q@f$
     */    
    virtual double deltaGAq(StandardModel::quark q) const;


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
    double cWB;///< The dimension-6 operator coefficient \f$c_{WB}\f$.
    double cH;///< The dimension-6 operator coefficient \f$c_{H}\f$.
    double cL1L1;///< The dimension-6 operator coefficient \f$(c_{LL})_{1111}\f$.
    double cL1L2;///< The dimension-6 operator coefficient \f$(c_{LL})_{1122}\f$.
    double cL1L3;///< The dimension-6 operator coefficient \f$(c_{LL})_{1133}\f$.
    double cL2L2;///< The dimension-6 operator coefficient \f$(c_{LL})_{2222}\f$.
    double cL2L3;///< The dimension-6 operator coefficient \f$(c_{LL})_{2233}\f$.
    double cL3L3;///< The dimension-6 operator coefficient \f$(c_{LL})_{3333}\f$.
    double cHL1p;///< The dimension-6 operator coefficient \f$(c_{HL}^\prime)_{11}\f$.
    double cHL2p;///< The dimension-6 operator coefficient \f$(c_{HL}^\prime)_{22}\f$.
    double cHL3p;///< The dimension-6 operator coefficient \f$(c_{HL}^\prime)_{33}\f$.
    double cHQ1p;///< The dimension-6 operator coefficient \f$(c_{HQ}^\prime)_{11}\f$.
    double cHQ2p;///< The dimension-6 operator coefficient \f$(c_{HQ}^\prime)_{22}\f$.
    double cHQ3p;///< The dimension-6 operator coefficient \f$(c_{HQ}^\prime)_{33}\f$.
    double cHL1;///< The dimension-6 operator coefficient \f$(c_{HL})_{11}\f$.
    double cHL2;///< The dimension-6 operator coefficient \f$(c_{HL})_{22}\f$.
    double cHL3;///< The dimension-6 operator coefficient \f$(c_{HL})_{33}\f$.
    double cHQ1;///< The dimension-6 operator coefficient \f$(c_{HQ})_{11}\f$.
    double cHQ2;///< The dimension-6 operator coefficient \f$(c_{HQ})_{22}\f$.
    double cHQ3;///< The dimension-6 operator coefficient \f$(c_{HQ})_{33}\f$.
    double cHE1;///< The dimension-6 operator coefficient \f$(c_{HE})_{11}\f$.
    double cHE2;///< The dimension-6 operator coefficient \f$(c_{HE})_{22}\f$.
    double cHE3;///< The dimension-6 operator coefficient \f$(c_{HE})_{33}\f$.
    double cHU1;///< The dimension-6 operator coefficient \f$(c_{HU})_{11}\f$.
    double cHU2;///< The dimension-6 operator coefficient \f$(c_{HU})_{22}\f$.
    double cHU3;///< The dimension-6 operator coefficient \f$(c_{HU})_{33}\f$.
    double cHD1;///< The dimension-6 operator coefficient \f$(c_{HD})_{11}\f$.
    double cHD2;///< The dimension-6 operator coefficient \f$(c_{HD})_{22}\f$.
    double cHD3;///< The dimension-6 operator coefficient \f$(c_{HD})_{33}\f$.
    double LambdaNP;///< The new physics scale \f$\Lambda_{NP}\f$.
    

    ////////////////////////////////////////////////////////////////////////   
private:

};

#endif	/* NPEFFECTIVE_H */

