/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RUNNER_H
#define	RUNNER_H

#include "ThObservable.h"
#include "THDM.h"
#include "THDMquantities.h"

/**
 * @class Runner
 * @ingroup THDM 
 * @brief An RGE running algorithm.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Renormalization group evolution of the relevant SM and THDM parameters
 */
class Runner : public ThObservable {
public:

    /**
     * @brief Runner constructor.
     */
    Runner(const StandardModel& SM_i);

    /**
     * @brief Runner destructor.
     */
    ~Runner();

    /**
     * @brief Empty function.
     */
    double computeThValue();

    virtual double RGERunner(/*int RGEs, const*/ double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2, int order, double Rpeps, double tNLOuni);

    const THDM * myTHDM;
    m11_2 * mym11_2;
    m22_2 * mym22_2;
    lambda1 * mylambda1;
    lambda2 * mylambda2;
    lambda3 * mylambda3;
    lambda4 * mylambda4;
    lambda5 * mylambda5;

};

#endif	/* RUNNER_H */
