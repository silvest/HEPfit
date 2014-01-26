/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MW_H
#define	MW_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class Mw 
 * @ingroup EW 
 * @brief A class for the \f$W\f$-boson mass.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the \f$W\f$-boson mass. The computation
 * includes the EW 2-loop corrections of \f${\cal O}(\alpha^2)\f$ as well as
 * the leading \f${\cal O}(G_\mu^2\alpha_s m_t^4)\f$ and \f${\cal O}(G_\mu^3m_t^6)\f$
 * contributions. We use the approximate formula in Ref. \cite  .
 */
class Mw : public ThObservable {
public:

    /**
     * A constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base class of 
     * the electroweak precision observables.
     */
    Mw(const EW& EW_i) 
    : ThObservable(EW_i)//, myEW(EW_i)
    {
    };

    /**
     * @return the \f$W\f$-boson mass
     */
    double computeThValue();

    
private:

    /**
     * A reference to an object of EW class, which is the base class of the electroweak 
     * precision observables.
     */
    //const EW& myEW;
};

#endif	/* MW_H */

