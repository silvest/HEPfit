/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWEPSILONS_H
#define	EWEPSILONS_H

#include "StandardModel.h"

/**
 * @class EWepsilons
 * @ingroup StandardModel
 * @brief A class for the @f$W@f$-boson mass and the @f$Zf\bar{f}@f$ effective
 * couplings with the epsilon parameters.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class EWepsilons {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EWepsilons(const StandardModel& SM_i) 
    : SM(SM_i) 
    {
    };


    ////////////////////////////////////////////////////////////////////////     
    
    double Mw(const double eps1, const double eps2, const double eps3) const;

    complex rhoZ_l(const StandardModel::lepton l, const double eps1) const;
    complex rhoZ_q(const StandardModel::quark q, const double eps1) const;
    complex kappaZ_l(const StandardModel::lepton l, 
                     const double eps1, const double eps3) const;
    complex kappaZ_q(const StandardModel::quark q, 
                     const double eps1, const double eps3) const;       
    complex gVl(const StandardModel::lepton l, 
                const double eps1, const double eps3) const;
    complex gVq(const StandardModel::quark q, 
                const double eps1, const double eps3) const;
    complex gAl(const StandardModel::lepton l, const double eps1) const;
    complex gAq(const StandardModel::quark q, const double eps1) const; 

    complex rhoZ_b(const double eps1, const double epsb) const;
    complex kappaZ_b(const double eps1, const double eps3, const double epsb) const;       
    complex gVb(const double eps1, const double eps3, const double epsb) const; 
    complex gAb(const double eps1, const double epsb) const; 
        
    
    ////////////////////////////////////////////////////////////////////////     
protected:
    double Delta_rW(const double eps1, const double eps2, const double eps3) const;
    double Delta_kappaPrime(const double eps1, const double eps3) const;
 
    complex rhoZ_e(const double eps1) const;
    complex kappaZ_e(const double eps1, const double eps3) const;
    complex gVe(const double eps1, const double eps3) const;
    complex gAe(const double eps1) const;    
    
    
    ////////////////////////////////////////////////////////////////////////     
private:
    const StandardModel& SM;///< A reference to an object of type StandardModel.

    
};

#endif	/* EWEPSILONS_H */

