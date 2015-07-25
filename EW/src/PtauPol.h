/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PTAUPOL_H
#define	PTAUPOL_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class PtauPol
 * @ingroup EW 
 * @brief An observable class for the @f$\tau@f$ polarization in 
 * @f$e^+e^-\to Z\to \tau^+\tau^-@f$ at the @f$Z@f$ pole.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the @f$\tau@f$ polarization in the 
 * @f$e^+e^-\to Z\to \tau^+\tau^-@f$ process at the @f$Z@f$ pole,
 * @f[
 * P_\tau^{\mathrm{pol}}=\mathrm{A}_\tau\,.
 * @f]
 * 
 * @sa EW_NPZff::PtauPol() and the detailed description of EW class
 * for the inclusion of new physics contribution
 *
 */
class PtauPol : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    PtauPol(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The @f$\tau@f$ polarization in @f$e^+e^-\to Z\to \tau^+\tau^-@f$
     * at the @f$Z@f$ pole, @f$P_\tau^{\mathrm{pol}}@f$.
     * @return @f$P_\tau^{\mathrm{pol}}@f$
     */
    double computeThValue();

    
private:


};

#endif	/* PTAUPOL_H */

