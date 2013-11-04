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
 * @brief A class for retrieving parameters in StandardModel.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class StandardModelParams : public ThObservable {
public:

    /**
     * @brief StandardModelParams constructor. 
     * @param[in] ObsType An object of ModelObservable class.
     * @param[in] name_i The name of a parameter: "AlsMz", "dAle5Mz", "Mz", "mtop",
     * "mHl", "delMw", "delSin2th_l", "delGammaZ", "delRhoZ_nu", "delRhoZ_e",
     * "delRhoZ_u", "delRhoZ_d", "delRhoZ_b".
     */
    StandardModelParams(const ThObsType& ObsType, const std::string name_i)
    : ThObservable(ObsType), name(name_i)
    {
    };

    /**
     * @brief Retrieves the value of the parameter specified in the constructor. 
     * @return The value of the parameter. 
     */
    double computeThValue();

private:
    const std::string name;

};

#endif	/* STANDARDMODELPARAMS_H */

