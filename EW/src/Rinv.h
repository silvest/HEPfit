/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef R_INV_H
#define R_INV_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class Rinv
 * @ingroup EW 
 * @brief An observable class for the ratio of the invisible decay width and 
 * the leptonic (electron) decay width of the @f$Z@f$ boson.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute ratio of the invisible decay width and
 * the leptonic (electron) decay width of the @f$Z@f$
 * boson,
 * @f[
 * R_{Inv}
 * = \frac{3\,\Gamma_\nu}{\Gamma_e}\,.
 * @f]
 *
 * @sa EW_NPZff::Rinv() and the detailed description of EW class
 * for the inclusion of new physics contribution
 * 
 */
class Rinv : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Rinv(const StandardModel& SM_i)
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The ratio of the invisible and leptonic (electron) decay widths
     * of the @f$Z@f$ boson, @f$\R_{inv}@f$.
     * @return @f$R_{inv}@f$
     */
    double computeThValue();

    
private:


};

#endif	/* R_INV_H */

