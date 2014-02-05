/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SIGMAHADRON_H
#define	SIGMAHADRON_H

#include <stdexcept>
#include <ThObservable.h>
#include "EW.h"

/**
 * @class sigmaHadron 
 * @ingroup EW 
 * @brief A class for the hadronic cross section of the @f$Z@f$ boson at the 
 * @f$Z@f$ pole (in units of nb)
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$e^+ e^-\rightarrow Z\rightarrow \mathrm{hadrons}@f$
 * cross section at the @f$Z@f$ pole (in units of nb), in terms of \f$\Gamma_e\f$, \f$\Gamma_h\f$ and 
 * \f$\Gamma_Z\f$ in the \b EW class (for the SM contribution):
 * \f[\sigma_h^0=\frac{12\pi}{M_Z^2}\frac{\Gamma_e\Gamma_h}{\Gamma_Z^2}.\f] 
 */
class sigmaHadron : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] EW_i A reference to an object of EW class, which is the base 
     * class of the electroweak precision observables.
     */
    sigmaHadron(const EW& EW_i) 
    : ThObservable(EW_i), myEW(EW_i)
    {
    };

    /**
     * @brief The theory prediction for the hadronic cross section at the \f$Z\f$ pole.
     * @return The cross section for the process \f$e^+ e^-\to Z\to \mathrm{hadrons}\f$
     * at the @f$Z@f$ pole in GeV\f$^{-2}\f$. 
     */
    double computeThValue();


private:

    /**
     * A reference to an object of EW class, which is the base class of the 
     * electroweak precision observables.
     */
    const EW& myEW;///< A reference to an object of the EW class.
};

#endif	/* SIGMAHADRON_H */

