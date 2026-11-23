/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2OBLIQUE_H
#define	LEP2OBLIQUE_H

#include "StandardModel.h"

/**
 * @class LEP2oblique
 * @ingroup EW
 * @brief A class for NP analyses of LEP-II observables with the extended oblique parameters. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Linear interference of the universal neutral-current corrections of
 * Barbieri et al. (hep-ph/0405040) with the SM Born amplitude at LEP2.
 * The reference uses the on-shell SM angle, the running photon coupling and
 * the G_mu normalization of the Z current. MWinput selects which of alpha
 * and M_W is held fixed in the input shifts. The SM width is held fixed.
 * QED radiation and weak loop corrections to the NP term are beyond this order.
 */
class LEP2oblique {
public:

    enum Oblique {Shat=0, That, Uhat, V, W, X, Y};
    
    /**
     * @brief LEP2oblique constructor
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    LEP2oblique(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////
    
    double sigma_l_LEP2_NP(const QCD::lepton l, const double s,
                           const double ml, const double ObParam_i[]) const;
    double sigma_q_LEP2_NP(const QCD::quark q, const double s,
                           const double mq, const double ObParam_i[]) const;

    double AFB_l_LEP2_NP(const QCD::lepton l, const double s, 
                         const double ml, const double ObParam_i[]) const;
    double AFB_q_LEP2_NP(const QCD::quark q, const double s, 
                         const double mq, const double ObParam_i[]) const;
    
    double R_q_LEP2_NP(const QCD::quark q, const double s, 
                       const double mq, const double ObParam_i[]) const;
        
    
    ////////////////////////////////////////////////////////////////////////    
private:
    const StandardModel& SM;

    struct Rates {
        double sigma, deltaSigma, forwardBackward, deltaForwardBackward;
        double deltaAFB() const;
    };

    Rates rates(const Particle& finalState, double s, double mass, double colors,
                const double oblique[]) const;
};

#endif /* LEP2OBLIQUE_H */
