/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPHIGGSST_H
#define	NPHIGGSST_H

#include "NPHiggs.h"

/**
 * @class NPHiggsST
 * @brief A class for new physics with a non-standard @f$HVV@f$ coupling. 
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class NPHiggsST : public NPHiggs {
public:
    
    NPHiggsST()
    : NPHiggs()
    {
    }

    virtual std::string ModelName() const 
    {
        return "NPHiggsST";
    }

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return NP contribution to oblique parameter S
     */
    virtual double obliqueS() const;
        
    /**
     * @return NP contribution to oblique parameter T
     */
    virtual double obliqueT() const;
    
    /**
     * @return NP contribution to oblique parameter U
     */
    virtual double obliqueU() const;
    
    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return epsilon_1
     */
    virtual double epsilon1() const;

    /**
     * @return epsilon_2
     */
    virtual double epsilon2() const;
    
    /**
     * @return epsilon_3
     */
    virtual double epsilon3() const;
    
    /**
     * @return epsilon_b
     */
    virtual double epsilonb() const;    

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return the W boson mass
     */
    virtual double Mw() const;

    /**
     * @return Mw^2/Mz^2
     */
    virtual double cW2() const;
    
    /**
     * @return 1-Mw^2/Mz^2
     */
    virtual double sW2() const;
    
    /**
     * @brief effective coupling rho_Z^l
     * @param[in] l lepton
     * @return rho_Z^l
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;
    
    /**
     * @brief effective coupling rho_Z^q
     * @param[in] q quark
     * @return rho_Z^q
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;

    /**
     * @brief effective coupling kappa_Z^l
     * @param[in] l name of lepton
     * @return kappa_Z^l in the SM
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @brief effective coupling kappa_Z^q
     * @param[in] q name of quark
     * @return kappa_Z^q in the SM
     */
    virtual complex kappaZ_q(const StandardModel::quark q) const;       
    
    /**
     * @brief vector effective coupling for neutral-current interactions
     * @param[in] l lepton
     * @return g_V^l
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @brief vector effective coupling for neutral-current interactions
     * @param[in] q quark
     * @return g_V^q
     */
    virtual complex gVq(const StandardModel::quark q) const;

    /**
     * @brief axial-vector effective coupling for neutral-current interactions
     * @param[in] l lepton
     * @return g_A^l
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @brief axial-vector effective coupling for neutral-current interactions
     * @param[in] q quark
     * @return g_A^q
     */
    virtual complex gAq(const StandardModel::quark q) const; 
    
    /**
     * @return the total width of the W boson
     */
    virtual double GammaW() const;    

    
    ////////////////////////////////////////////////////////////////////////     
    
protected:    
    
private:

};

#endif	/* NPHIGGSST_H */

