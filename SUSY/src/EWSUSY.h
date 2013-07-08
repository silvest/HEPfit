/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSUSY_H
#define	EWSUSY_H

#include <gslpp.h>
#include "SUSY.h"

using namespace gslpp;

/**
 * @class EWSUSY
 * @ingroup SUSY
 * @brief A class for SUSY contributions to the EW precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class EWSUSY {
public:

    /**
     * @brief A constructor.
     * @param[in] SUSY_in A reference to a SUSY object. 
     */
    EWSUSY(const SUSY& SUSY_in);

    /**
     * @brief Sets parameters in Rosiek's notation. 
     */
    void SetRosiekParameters();

private:

    /* A reference to the SUSY object passed to the constructor. */
    const SUSY& mySUSY;

    /* Yukawa couplings in Janusz Rosiek's notation */
    matrix<complex> Yu_JR, Yd_JR, Yl_JR;

    /* Trilinear couplings in Janusz Rosiek's notation */
    matrix<complex> Au_JR, Ad_JR, Al_JR;

};

#endif	/* EWSUSY_H */

