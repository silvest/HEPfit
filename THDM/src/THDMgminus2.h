/*
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMGMINUS2_H
#define	THDMGMINUS2_H

#include "ThObservable.h"

/**
 * @class THDMgminus2_mu
 * @ingroup THDM
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class THDMgminus2_mu : public ThObservable {
public:

    /**
     * @brief Constructor of the class THDMgminus2_mu
     */
    THDMgminus2_mu(const StandardModel& SM_i);

    /**
     * @return value of \f$ (g-2)_{\mu} \f$
     */
    double computeThValue();

private:

};

#endif	/* THDMGMINUS2_H */
