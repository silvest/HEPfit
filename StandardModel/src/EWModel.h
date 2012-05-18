/* 
 * File:   EWModel.h
 * Author: mishima
 */

#ifndef EWMODEL_H
#define	EWMODEL_H

#include <gslpp.h>
#include "StandardModel.h"
using namespace gslpp;


class EWModel {
public:
    
    /**
     * @return the W boson mass
     */
    virtual double Mw() const = 0;

    /**
     * @return Mw^2/Mz^2
     */
    virtual double cW2() const = 0;
    
    /**
     * @return 1-Mw^2/Mz^2
     */
    virtual double sW2() const = 0;
    
    /**
     * @brief effective coupling rho_Z^l
     * @param[in] l lepton
     * @return rho_Z^l
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const = 0;
    
    /**
     * @brief effective coupling rho_Z^q
     * @param[in] q quark
     * @return rho_Z^q
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const = 0;
    
    /**
     * @brief vector effective coupling for neutral-current interactions
     * @param[in] l lepton
     * @return g_V^l
     */
    virtual complex gVl(const StandardModel::lepton l) const = 0;

    /**
     * @brief vector effective coupling for neutral-current interactions
     * @param[in] q quark
     * @return g_V^q
     */
    virtual complex gVq(const StandardModel::quark q) const = 0;

    /**
     * @brief axial-vector effective coupling for neutral-current interactions
     * @param[in] l lepton
     * @return g_A^l
     */
    virtual complex gAl(const StandardModel::lepton l) const = 0;

    /**
     * @brief axial-vector effective coupling for neutral-current interactions
     * @param[in] q quark
     * @return g_A^q
     */
    virtual complex gAq(const StandardModel::quark q) const = 0;    
    
    /**
     * @return the total width of the W boson
     */
    virtual double GammaW() const = 0;
    
    /**
     * @return NP contribution to oblique parameter S
     */
    virtual double obliqueS() const = 0;
        
    /**
     * @return NP contribution to oblique parameter T
     */
    virtual double obliqueT() const = 0;
    
    /**
     * @return NP contribution to oblique parameter U
     */
    virtual double obliqueU() const = 0;
    
    
};

#endif	/* EWMODEL_H */

