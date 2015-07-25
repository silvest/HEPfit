/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SIGMAHADRON_H
#define	SIGMAHADRON_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class sigmaHadron 
 * @ingroup EW 
 * @brief An observable class for the cross section of 
 * @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$ at the @f$Z@f$ pole
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the cross section for
 * @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$ at the @f$Z@f$ pole, 
 * in terms of @f$\Gamma_e@f$, @f$\Gamma_h@f$ and @f$\Gamma_Z@f$ in the EW class: 
 * @f[
 * \sigma_h^0 = \frac{12\pi}{M_Z^2}\frac{\Gamma_e\Gamma_h}{\Gamma_Z^2}\,.
 * @f]
 *
 * @sa EW_NPZff::sigmaHadron() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class sigmaHadron : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    sigmaHadron(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The cross section for the process @f$e^+ e^-\to Z\to \mathrm{hadrons}@f$
     * at the @f$Z@f$ pole, @f$\sigma_h^0@f$, in units of nb.
     * @return @f$\sigma_h^0@f$ in units of nb. 
     */
    double computeThValue();


private:


};

#endif	/* SIGMAHADRON_H */

