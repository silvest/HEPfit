/*
 * Copyright (C) 2026 HEPfit Collaboration
 * For the licensing terms see doc/COPYING.
 */

#ifndef SCALEDOBS_H
#define SCALEDOBS_H

#include "ThObservable.h"
#include <memory>
#include <string>

/**
 * @class ScaledObs
 * @ingroup Observable
 * @brief A theory observable multiplied by a model parameter
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Requested in the configuration file with the theory name "<parameter>*<observable>", e.g.
 * "r_bsg*BR_bsgamma": the value is the model parameter times the observable, which receives the bin limits of
 * this one. Used for multiplicative theory-uncertainty factors (nuisance parameters) on a prediction.
 */
class ScaledObs : public ThObservable {
public:

    /**
     * @brief Constructor.
     * @param[in] SM a reference to the model
     * @param[in] parameter the name of the model parameter
     * @param[in] obs the observable to be scaled (owned by this object)
     */
    ScaledObs(const StandardModel& SM, const std::string parameter, ThObservable* obs);

    /**
     * @return the parameter times the value of the observable, computed in the bin of this observable
     */
    double computeThValue();

private:
    const double& scale; ///< The model parameter.
    std::unique_ptr<ThObservable> obs; ///< The observable that is scaled.
};

#endif /* SCALEDOBS_H */
