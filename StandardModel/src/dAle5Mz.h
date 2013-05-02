/* 
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef DALE5MZ_H
#define	DALE5MZ_H

#include <ThObservable.h>
#include <ThObsType.h>
#include "StandardModel.h"

/**
 * @class dAle5Mz
 * @ingroup StandardModel
 * @brief A class for an interface to the input parameter @f$\Delta\alpha_{\rm had}^{(5)}(M_Z^2)@f$.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class dAle5Mz : public ThObservable {
public:

    dAle5Mz(const ThObsType& ObsType) 
    : ThObservable(ObsType) 
    {
    };
    
    double getThValue()
    {
        return SM.getDAle5Mz();
    };

private:

};

#endif	/* DALE5MZ_H */

