/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RBOTTOM_H
#define	RBOTTOM_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Rbottom
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_b^0=\Gamma(Z\to b\bar{b})/\Gamma(Z\to {\rm hadrons})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to b\bar{b}@f$ 
 * width to the @f$Z@f$-boson hadronic width:
 * @f[
 * R_b = \frac{\Gamma_b}{\Gamma_h}\,.
 * @f]
 *
 * @sa EW_NPZff::Rbottom() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Rbottom : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Rbottom(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_b^0=\Gamma(Z\to b\bar{b})/\Gamma(Z\to {\rm hadrons})@f$.
     * @return @f$R_b^0@f$
     */
    double computeThValue();

    
private:


};

#endif	/* RBOTTOM_H */

