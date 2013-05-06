/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MTPOLE_H
#define	MTPOLE_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "QCD.h"

/**
 * @class mtpole
 * @ingroup StandardModel
 * @brief A class for an interface to the input parameter @f$m_t@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class mtpole : public ThObservable {
public:

    mtpole(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.getMtpole();
    };

private:

};

#endif	/* MTPOLE_H */

