/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RWC_H
#define	RWC_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class RWc
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_{Wc}=\Gamma(W\to cs)/\Gamma(W\to {\rm hadrons})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to c\bar{c}@f$ 
 * width to the @f$Z@f$-boson hadronic width:
 * @f[
 * R_{Wc}=\Gamma(W\to cs)/\Gamma(W\to {\rm hadrons})\,.
 * @f]
 *
 * @sa EW_NPZff::RWc() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class RWc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    RWc(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_{Wc}=\Gamma(W\to cs)/\Gamma(W\to {\rm hadrons})@f$.
     * @return @f$R_{Wc}@f$
     */
    double computeThValue();

    
private:


};

#endif	/* RWC_H */

