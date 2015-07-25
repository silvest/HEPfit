/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALEPTON_H
#define	ALEPTON_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Alepton
 * @ingroup EW 
 * @brief An observable class for the left-right asymmetry in
 * @f$e^+e^-\to Z\to \ell^+ \ell^-@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the left-right asymmetry for an 
 * @f$e^+e^-\to Z\to \ell^+ \ell^-@f$ process at the @f$Z@f$ pole, where
 * @f$\ell@f$ denotes a charged lepton, and lepton-flavour universality
 * is assumed.
 *
 * @sa EW_NPZff::Alepton() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Alepton : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Alepton(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The left-right asymmetry for @f$e^+e^-\to Z\to \ell^+ \ell^-@f$
     * at the @f$Z@f$ pole, @f$\mathcal{A}_\ell@f$.
     * @return @f$\mathcal{A}_\ell@f$
     */
    double computeThValue();

    
private:


};


#endif	/* ALEPTON_H */

