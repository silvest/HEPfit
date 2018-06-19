/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMEWPO_H
#define	GENERALTHDMEWPO_H

#include "ThObservable.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"
#include "StandardModel.h"
//#include "gslpp.h"

/**
 * @class GeneralTHDMEWPO
 * @ingroup GeneralTHDM
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */

/**
 * @class Rb0
 * @ingroup GeneralTHDM
 * @brief 
 */
class  Rb0GTHDM: public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Rb0GTHDM(const StandardModel& SM_i);

    /**
     * @return Rb0GTHDM
     */
    double computeThValue ();
private:
    const GeneralTHDM& myGTHDM;
};

#endif	/* GENERALTHDMEWPO_H */
