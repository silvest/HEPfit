/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BETASJPSIPHI_H
#define	BETASJPSIPHI_H

#include "ThObservable.h"

/**
 * @class Phis_JPsiPhi
 * @ingroup Flavour
 * @brief A class for the @f$B_s@f$ mixing phase @f$\phi_s@f$ measured in @f$B_s\to J/\psi\phi@f$
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * @f$\phi_s = -\mathrm{arg}(M_{21}^{bs}) + 2\,\mathrm{arg}(V_{cs}V_{cb}^*) - 2\phi_{B_s}@f$
 * (reduced to @f$[-\pi,\pi]@f$), with @f$\phi_s^\mathrm{SM} = -2\beta_s@f$.
 * The new-physics phase StandardModel::getPhiBs() therefore enters as
 * @f$\phi_s = \phi_s^\mathrm{SM} - 2\phi_{B_s}@f$, corresponding to
 * @f$(M_{12}^{bs})_\mathrm{full} = C_{B_s}\, e^{-2i\phi_{B_s}} (M_{12}^{bs})_\mathrm{SM}@f$;
 * this is the opposite sign convention w.r.t. @f$\phi_{B_d}@f$ in the @f$B_d@f$ observables.
 */
class Phis_JPsiPhi : public ThObservable{
public:
    
   /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    Phis_JPsiPhi(const StandardModel& SM_i);
    
    /**
     *
     * @return theoretical value of @f$\phi_s@f$
     */
    virtual double computeThValue();

};

#endif	/* BETASJPSIPHI_H */

