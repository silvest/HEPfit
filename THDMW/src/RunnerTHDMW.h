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

    const THDMW * myTHDMW;
};

#endif	/* RUNNERTHDMW_H */
