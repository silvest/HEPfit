/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDB2_H
#define	AMPDB2_H

#include <gslpp_complex.h>
#include "Flavour.h"
#include <StandardModel.h>

using namespace gslpp;

/**
 * @addtogroup Flavour
 * @brief A module for flavour observables.
 * @{
 */

class AmpDB2 {
public:
    AmpDB2(const StandardModel& SM_i);

protected:
    complex AmpBd(orders order);
    complex AmpBs(orders order);

private:
    
    const StandardModel& mySM;

};

/**
 * @}
 */

#endif	/* AMPDB2_H */

