/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef AMPDB2_H
#define	AMPDB2_H

#include <gslpp_complex.h>
#include "Flavour.h"

using namespace gslpp;

/**
 * @addtogroup Flavour
 * @brief A module for flavour observables.
 * @{
 */

class AmpDB2 {
public:
    AmpDB2(Flavour& Flavour);

protected:
    complex AmpBd(orders order);
    complex AmpBs(orders order);

private:
    Flavour& myFlavour;

};

/**
 * @}
 */

#endif	/* AMPDB2_H */

