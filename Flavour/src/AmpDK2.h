/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDK2_H
#define	AMPDK2_H

class StandardModel;
#include "gslpp.h"
#include "OrderScheme.h"
#include "gslpp.h"


/**
 * @class AmpDK2
 * @ingroup Flavour
 * @brief A class for calculating the amplitudes contributing to @f$\epsilon_K@f$
 * and @f$\Delta M_K@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the theoretical value of
 * amplitudes contributing to @f$\epsilon_K@f$ and @f$\Delta M_K@f$. The
 * hadronic matrix elements are defined for the operators @f$O_1, \ldots, O_5@f$
 * in the chiral limit as can be found in ...
 *
 */
class AmpDK2 {
public:
    /**
    * @brief Constructor.
    * @param[in] SM_i a reference to an object of type StandardModel
    */
    AmpDK2(const StandardModel& SM_i);

protected:
    /**
     * 
     * @brief compute the amplitude for kaon oscillations
    * @param[in] order the %QCD order of the computation
     */
    gslpp::complex AmpDK(orders order);
    
    /**
     * 
     * @brief compute the NP part of the amplitude for kaon oscillations
    * @param[in] order the %QCD order of the computation
     */
    gslpp::complex AmpDMKNP(orders order);

private:
    
    const StandardModel& mySM;

};

#endif	/* AMPDK2_H */
