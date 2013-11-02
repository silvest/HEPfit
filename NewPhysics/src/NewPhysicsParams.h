/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NEWPHYSICSPARAMS_H
#define	NEWPHYSICSPARAMS_H

#include <string>
#include <ThObservable.h>
#include <ThObsType.h>

/**
 * @class NewPhysicsParams
 * @ingroup NewPhysics
 * @brief A class for posteriors of model parameters in NewPhysics.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class NewPhysicsParams : public ThObservable {
public:

    NewPhysicsParams(const ThObsType& ObsType, const std::string name_i)
    : ThObservable(ObsType), name(name_i)
    {
    };

    double computeThValue();

private:
    const std::string name;

};

#endif	/* NEWPHYSICSPARAMS_H */

