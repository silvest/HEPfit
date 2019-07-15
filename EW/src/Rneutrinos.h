/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RNEUTRINOS_H
#define	RNEUTRINOS_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Rneutrinos
 * @ingroup EW 
 * @brief An observable class for
 * @f$R_\nu^0=\sum \Gamma(Z\to \nu \bar{\nu})/\Gamma(Z\to {\rm hadrons})/@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the ratio of the @f$Z@f$-boson 
 * invisible and hadronic widths:
 * @f[
 * R_\nu = \sum \Gamma(Z\to \nu \bar{\nu})/\Gamma(Z\to {\rm hadrons})\,,
 * @f]
 * where @f$\nu@f$ denotes a neutrino.
 *
 * @sa EW_NPZff::Rneutrinos() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class Rneutrinos : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Rneutrinos(const StandardModel& SM_i)
    : ThObservable(SM_i)
    {
    };

    /**
     * @brief The ratio @f$R_\nu^0=\sum \Gamma(Z\to \nu \bar{\nu})/\Gamma(Z\to {\rm hadrons})@f$.
     * @return @f$R_\nu^0@f$
     */
    double computeThValue();


private:


};

#endif	/* RNEUTRINOS_H */

