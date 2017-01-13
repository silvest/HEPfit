/*
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMGMINUS2_H
#define	GENERALTHDMGMINUS2_H

#include "gslpp.h"
#include "ThObservable.h"
#include "LeptonFlavour.h"

/**
 * @class GeneralTHDMgminus2_mu
 * @ingroup GeneralTHDM
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class GeneralTHDMgminus2_mu : public ThObservable {
public:
    
    /**
     * @brief Constructor of the class GTHDMgminus2_mu
     */
    GeneralTHDMgminus2_mu(const StandardModel& SM_i);
    
    /**
     * @return value of \f$ (g-2)_{\mu} \f$
     */
    double computeThValue();
    
private:
    /**
     * @brief Constructor containing the Wilson coefficient 
     */
    const StandardModel& mySM;

};

#endif	/* GENERALTHDMGMINUS2_H */
