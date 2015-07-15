/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RLEPTON_H
#define	RLEPTON_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Rlepton
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_\ell^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to \ell^+ \ell^-)@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z@f$-boson 
 * hadronic width to the @f$Z\to \ell^+ \ell^-@f$ width:
 * @f[
 * R_\ell = \frac{\Gamma_h}{\Gamma_\ell}\,,
 * @f]
 * where @f$\ell@f$ denotes a charged lepton, and lepton-flavour universality
 * is assumed.
 *
 * @sa EW_NPZff::Rlepton() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Rlepton : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Rlepton(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The ratio @f$R_\ell^0=\Gamma(Z\to {\rm hadrons})/\Gamma(Z\to \ell^+ \ell^-)@f$.
     * @return @f$R_\ell^0@f$
     */
    double computeThValue();


private:


};

#endif	/* RLEPTON_H */

