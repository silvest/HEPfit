/*
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef PARAMOBS_H
#define	PARAMOBS_H

#include "ThObservable.h"

/**
 * @class ParamObs
 * @ingroup Observable
 * @brief A class for setting a parameter as an observable
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class sets a parameter as an observable when it is requested in the config file.
 */
class ParamObs : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     * @param[in] SM a reference to the model
     * @param[in] name the name of the parameter
     */
    ParamObs(const StandardModel& SM, const std::string name);

    /**
     * @brief the method to compute the theory value of the parameter
     * @return the value of the parameter
     */
    double computeThValue();
    
private:
    const double& param;///< The parameter.
};

#endif	/* PARAMOBS_H */

