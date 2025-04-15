/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RSTRANGE_H
#define	RSTRANGE_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Rstrange
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_s^0=\Gamma(Z\to s\bar{s})/\Gamma(Z\to {\rm hadrons})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to s\bar{s}@f$ 
 * width to the @f$Z@f$-boson hadronic width:
 * @f[
 * R_s = \frac{\Gamma_s}{\Gamma_h}\,.
 * @f]
 *
 * @sa EW_NPZff::Rstrange() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Rstrange : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Rstrange(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_s^0=\Gamma(Z\to s\bar{s})/\Gamma(Z\to {\rm hadrons})@f$.
     * @return @f$R_s^0@f$
     */
    double computeThValue();

    
private:


};

#endif	/* RSTRANGE_H */

