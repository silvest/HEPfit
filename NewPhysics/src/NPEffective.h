/* 
 * Copyright (C) 2013 SusyFit Collaboration
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
 * @brief A base class for new physics with the effective Lagrangian approach.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPEffective : public NPbase {
public:

    /**
     * @brief NPEffective constructor. 
     */
    NPEffective();

    virtual std::string ModelName() const
    {
        return "NPEffective";
    }

    virtual bool InitializeModel();
    virtual void setEWSMflags(EWSM& myEWSM);

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    virtual bool setFlag(const std::string, const bool&); 
    virtual bool CheckFlags() const;
    

    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}^\prime\right)_{11}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L^1}\sigma_a\gamma^\mu L^1\right)\f$
     */
    double getCHL1p() const {
        return cHL1p;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}^\prime\right)_{22}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L^2}\sigma_a\gamma^\mu L^2\right)\f$
     */
    double getCHL2p() const {
        return cHL2p;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}^\prime\right)_{33}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{L^3}\sigma_a\gamma^\mu L^3\right) \f$
     */
    double getCHL3p() const {
        return cHL3p;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{L^1}\gamma^\mu L^1\right) \f$
     */
    double getCHL1() const {
        return cHL1;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{L^2}\gamma^\mu L^2\right) \f$
     */
    double getCHL2() const {
        return cHL2;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HL}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{L^3}\gamma^\mu L^3\right) \f$
     */
    double getCHL3() const {
        return cHL3;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}^\prime\right)_{11}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q^1}\sigma_a\gamma^\mu Q^1\right)\f$
     */
    double getCHQ1p() const
    {
        return cHQ1p;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}^\prime\right)_{22}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q^2}\sigma_a\gamma^\mu Q^2\right)\f$
     */
    double getCHQ2p() const
    {
        return cHQ2p;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}^\prime\right)_{33}=i\left(H^\dagger\sigma_a D_\mu H\right)\left(\overline{Q^3}\sigma_a\gamma^\mu Q^3\right)\f$
     */
    double getCHQ3p() const
    {
        return cHQ3p;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q^1}\gamma^\mu Q^1\right) \f$
     */
    double getCHQ1() const
    {
        return cHQ1;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q^2}\gamma^\mu Q^2\right) \f$
     */    double getCHQ2() const
    {
        return cHQ2;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HQ}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{Q^3}\gamma^\mu Q^3\right) \f$
     */
    double getCHQ3() const
    {
        return cHQ3;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HE}\right)_{11}=i\left(H^\dagger D_\mu H\right)\left(\overline{E^1}\gamma^\mu E^1\right) \f$
     */
    double getCHE1() const {
        return cHE1;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HE}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{E^2}\gamma^\mu E^2\right) \f$
     */
    double getCHE2() const {
        return cHE2;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HE}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{E^3}\gamma^\mu E^3\right) \f$
     */
    double getCHE3() const {
        return cHE3;
    }
  
    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HU}\right)_{22}=i\left(H^\dagger D_\mu H\right)\left(\overline{U^2}\gamma^\mu U^2\right) \f$
     */  
    double getCHU2() const {
        return cHU2;
    }

    /**
     * @return the coefficient of the dimension-six effective operator 
     * \f$\left({\cal O}_{HD}\right)_{33}=i\left(H^\dagger D_\mu H\right)\left(\overline{D^3}\gamma^\mu D^3\right) \f$
     */
    double getCHD3() const {
        return cHD3;
    }

    /**
     * @return the value of the new physics scale \f$\Lambda_\mathrm{NP}\f$ [GeV] 
     */
    double getLambdaNP() const
    {
        return LambdaNP;
    }


    ////////////////////////////////////////////////////////////////////////

    virtual double v() const;

    virtual double Mw_tree() const;

    virtual double DeltaGF() const;


    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return the oblique parameter S
     */
    virtual double obliqueS() const;

    /**
     * @return the oblique parameter T
     */
    virtual double obliqueT() const;

    /**
     * @return the oblique parameter U
     */
    virtual double obliqueU() const;


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current left-handed coupling @f$g_L^l@f$
     */
    double deltaGLl(StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current left-handed coupling @f$g_L^q@f$
     */
    double deltaGLq(StandardModel::quark q) const;

    /**
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current right-handed coupling @f$g_R^l@f$
     */    
    double deltaGRl(StandardModel::lepton l) const;
    
   /**
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current right-handed coupling @f$g_R^q@f$
     */
    double deltaGRq(StandardModel::quark q) const;

  /**
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current vector coupling @f$g_A^l@f$
     */
    virtual double deltaGVl(StandardModel::lepton l) const;

   /**
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current vector coupling @f$g_V^q@f$
     */    
    virtual double deltaGVq(StandardModel::quark q) const;

   /**
     * @param[in] l name of a lepton
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^l@f$
     */    
    virtual double deltaGAl(StandardModel::lepton l) const;

   /**
     * @param[in] q name of a quark
     * @return the new physics correction to the neutral-current axial-vector coupling @f$g_A^q@f$
     */    
    virtual double deltaGAq(StandardModel::quark q) const;


    ////////////////////////////////////////////////////////////////////////

    virtual double epsilon1() const;

    virtual double epsilon2() const;

    virtual double epsilon3() const;

    virtual double epsilonb() const;


    ////////////////////////////////////////////////////////////////////////

    /**
     * @return the \f$W\f$-boson mass in GeV
     */
    virtual double Mw() const;

    /**
     * @return the (square of the) cosine of the weak angle in the On-mass-shell renormalization scheme,
     *  \f$\cos^2{\theta_W}=\frac{M_W^2}{M_Z^2}
     */
    virtual double cW2() const;

    /**
     * @return the (square of the) sine of the weak angle in the On-mass-shell renormalization scheme,
     *  \f$\sin^2{\theta_W}=1-\frac{M_W^2}{M_Z^2}
     */
    virtual double sW2() const;

    /**
     * @return the total width of the \f$W\f$ boson in GeV
     */
    virtual double GammaW() const;
    
    
    ////////////////////////////////////////////////////////////////////////    
protected:
    double cWB, cH;
    double cL1L1, cL1L2, cL1L3, cL2L2, cL2L3, cL3L3;
    double cHL1p, cHL2p, cHL3p;
    double cHQ1p, cHQ2p, cHQ3p;
    double cHL1, cHL2, cHL3;
    double cHQ1, cHQ2, cHQ3;
    double cHE1, cHE2, cHE3;
    double cHU1, cHU2, cHU3;
    double cHD1, cHD2, cHD3;
    double LambdaNP;
    

    ////////////////////////////////////////////////////////////////////////   
private:

};

#endif	/* NPEFFECTIVE_H */

