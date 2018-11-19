/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef N_NU_H
#define N_NU_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Nneutrinos
 * @ingroup EW 
 * @brief An observable class for the number of neutrinos obtained indirectly
 * from the measurements at the Z pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the number of neutrinos obtained indirectly
 * from the measurements at the Z pole,
 * @f[
 * N_{\nu}
 * = \frac{\Gamma_\ell}{\Gamma_{\nu}}\left(\sqrt{\frac{12\pi R_\ell}{M_Z^2 \sigma_\mathrm{had}^0}}-R_\ell - 3\right)\,.
 * @f]
 *
 * 
 */
class Nneutrinos : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Nneutrinos(const StandardModel& SM_i)
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The number of neutrinos obtained indirectly
     * from the measurements at the Z pole, @f$N_{\nu}@f$.
     * @return @f$N_{\nu}@f$
     */
    double computeThValue();

    
private:


};

#endif	/* N_NU_H */

