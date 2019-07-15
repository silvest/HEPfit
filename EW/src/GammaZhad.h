/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GAMMAZHAD_H
#define	GAMMAZHAD_H

#include <stdexcept>
#include <ThObservable.h>

/**
 * @class GammaZhad
 * @ingroup EW 
 * @brief An observable class for the hadronic decay width of the @f$Z@f$ boson. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the hadronic decay width of the @f$Z@f$
 * boson,
 * @f$\Gamma_{Z,had}=\sum_{q\neq t}\Gamma_q@f$.
 *
 * @sa EW_NPZff::GammaZhad() and the detailed description of EW class
 * for the inclusion of new physics contribution
 * 
 */
class GammaZhad : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    GammaZhad(const StandardModel& SM_i) 
    : ThObservable(SM_i) 
    {
    };

    /**
     * @brief The hadronic decay width of the @f$Z@f$ boson, @f$\Gamma_{Z,had}@f$,
     * in units of GeV.
     * @return @f$\Gamma_{Z,had}@f$ in units of GeV
     */
    double computeThValue();

    
private:


};

#endif	/* GAMMAZHAD_H */

