/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GMUNITARITY_H
#define	GMUNITARITY_H

#include "ThObservable.h"
#include "GeorgiMachacek.h"
#include "GMcache.h"
//#include <gslpp.h>

/**
 * @class GMunitarityLO
 * @ingroup GeorgiMachacek
 * @brief An observable class for the requirement of perturbative unitarity at leading order.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to require unitarity for all the tree level 
 * scalar-scalar scattering amplitudes.
 * The eigenvalues of the S-matrix can be found in @cite ?.
 * They should be smaller than ? in magnitude to preserve the unitarity of the S-matrix.
 */
class GMunitarityLO: public ThObservable {
public:

    /**
     * @brief GMunitarityLO constructor.
     */
    GMunitarityLO(const StandardModel& SM_i, unsigned int index_i);

    /**
     * @brief Destructor.
     */
    virtual ~GMunitarityLO();

    /**
     * @return Unitarity eigenvalues
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
    unsigned int index;
};

#endif	/* GMUNITARITY_H */
