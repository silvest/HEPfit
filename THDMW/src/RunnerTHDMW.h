/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RUNNERTHDMW_H
#define	RUNNERTHDMW_H

#include "THDMW.h"

/**
 * @addtogroup THDMW
 * @brief A module for the extension of either the Standard Model or the @f$Z_2@f$ symmetric Two-Higgs-Doublet model by a scalar octet.
 * @details The Standard Model extension by a scalar octet can be accessed setting the THDMWmodel flag to "ManoharWise".
 * A custodial limiting case of this model can be obtained setting the flag to "custodialMW".
 * The THDM plus octet can be address by attributing "custodial1" to the THDMWmodel flag.
 * The implemented observables are positivity, unitarity, Higgs signal strengths and electroweak precision measurements.
 * (Not all observables are defined for all model types!)
 * 
 * @{
 */

/**
 * @class RunnerTHDMW
 * @ingroup THDMW
 * @brief An RGE running algorithm for the THDMW parameters.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Renormalization group evolution of the relevant SM and THDMW parameters
 */
class RunnerTHDMW {
public:

    /**
     * @brief RunnerTHDMW constructor.
     */
    RunnerTHDMW(const StandardModel& SM_i);

    /**
     * @brief Runner destructor.
     */
    virtual ~RunnerTHDMW();

    virtual double RGERunnerTHDMW(/*int RGEs, const*/ double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order, double Rpeps, double tNLOuni);
    virtual double RGERunnerMW(/*int RGEs, const*/ double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order, double Rpeps, double tNLOuni);
    virtual double RGERunnercustodialMW(/*int RGEs, const*/ double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order, double Rpeps, double tNLOuni);

    const THDMW * myTHDMW;
};

/**
 * @}
 */

#endif	/* RUNNERTHDMW_H */
