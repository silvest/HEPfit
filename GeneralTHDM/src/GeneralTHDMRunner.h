/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMRUNNER_H
#define	GENERALTHDMRUNNER_H

#include "GeneralTHDM.h"

/**
 * @class GeneralTHDMRunner
 * @ingroup GeneralTHDM
 * @brief An RGE running algorithm for the GeneralTHDM parameters.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Renormalization group evolution of the relevant SM and GeneralTHDM parameters
 */
class GeneralTHDMRunner {
public:

    /**
     * @brief GeneralTHDMRunner constructor.
     */
    GeneralTHDMRunner(const StandardModel& SM_i);

    /**
     * @brief Runner destructor.
     */
    virtual ~GeneralTHDMRunner();

    virtual double RGEGeneralTHDMRunner(/*int RGEs, const*/ double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order, double Rpeps, double tNLOuni);

    const GeneralTHDM * myGTHDM;
};

#endif	/* GENERALTHDMRUNNER_H */
