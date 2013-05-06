/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ALSMZ_H
#define	ALSMZ_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "QCD.h"

/**
 * @addtogroup StandardModel 
 * @brief A project for Standard Model. 
 * @{
 */

/**
 * @class alsMz
 * @brief A class for an interface to the input parameter @f$\alpha_s(M_Z^2)@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class alsMz : public ThObservable {
public:

    alsMz(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.getAlsMz();
    };

private:

};

/** 
 * @}
 */

#endif	/* ALSMZ_H */

