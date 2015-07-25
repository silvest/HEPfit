/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RCHARM_H
#define	RCHARM_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Rcharm
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_c^0=\Gamma(Z\to c\bar{c})/\Gamma(Z\to {\rm hadrons})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to c\bar{c}@f$ 
 * width to the @f$Z@f$-boson hadronic width:
 * @f[
 * R_c = \frac{\Gamma_c}{\Gamma_h}\,.
 * @f]
 *
 * @sa EW_NPZff::Rcharm() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Rcharm : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Rcharm(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_c^0=\Gamma(Z\to c\bar{c})/\Gamma(Z\to {\rm hadrons})@f$.
     * @return @f$R_c^0@f$
     */
    double computeThValue();

    
private:


};

#endif	/* RCHARM_H */

