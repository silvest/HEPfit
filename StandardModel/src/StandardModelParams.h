/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef STANDARDMODELPARAMS_H
#define	STANDARDMODELPARAMS_H

#include <string>
#include <ThObservable.h>
#include <ThObsType.h>

/**
 * @class StandardModelParams
 * @ingroup StandardModel
 * @brief A class for posteriors of model parameters in StandardModel.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class StandardModelParams : public ThObservable {
public:

    StandardModelParams(const ThObsType& ObsType, const std::string name_i)
    : ThObservable(ObsType), name(name_i)
    {
    };

    double computeThValue();

private:
    const std::string name;

};

#endif	/* STANDARDMODELPARAMS_H */

