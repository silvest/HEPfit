/*
 * Copyright (C) 2013-2014 SusyFit Collaboration
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
 * @brief A class for retrieving parameters associated with NewPhysics project.
 * @ingroup NewPhysics
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class NewPhysicsParams : public ThObservable {
public:

    /**
     * @brief Constructor. 
     * @param[in] ObsType
     * @param[in] name_i the name of the parameter to be retrieved
     */
    NewPhysicsParams(const ThObsType& ObsType, const std::string name_i)
    : ThObservable(ObsType), name(name_i)
    {
    };

    /**
     * @brief A method to retrieve the value of the desired new physics parameter.
     * @return the value of the new physics parameter
     */
    double computeThValue();

private:
    const std::string name;

};

#endif	/* NEWPHYSICSPARAMS_H */

