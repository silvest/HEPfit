/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RUC_H
#define	RUC_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Ruc
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_{uc}=(\Gamma(Z\to\ u\bar{u})+\Gamma(Z\to\ c\bar{c}))/2\Gamma(Z\to {\rm hadrons})@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z\to u\bar{u} + Z\to c\bar{c}@f$ 
 * width to the @f$Z@f$-boson hadronic width:
 * @f[
 * R_{uc} = \frac{\Gamma_u+\Gamma_c}{\Gamma_h}\,.
 * @f]
 *
 * @sa EW_NPZff::Ruc() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Ruc : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Ruc(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio @f$R_{uc}^0=(\Gamma(Z\to\ u\bar{u})+\Gamma(Z\to\ c\bar{c}))/2\Gamma(Z\to {\rm hadrons})@f$.
     * @return @f$R_{uc}^0@f$
     */
    double computeThValue();

    
private:


};

#endif	/* RUC_H */

