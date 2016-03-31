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
 * @details ?
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

    virtual double RGERunner(/*int RGEs, const*/ double InitialValues[], unsigned long int NumberOfRGEs, double Q1, double Q2);

    const THDM * myTHDM;
    m11_2 * mym11_2;
    m22_2 * mym22_2;
    lambda1 * mylambda1;
    lambda2 * mylambda2;
    lambda3 * mylambda3;
    lambda4 * mylambda4;
    lambda5 * mylambda5;

private:

//    virtual int RGEs(double t, const double y[], double beta[], void *params);

};

/**
 * @class lambda1atQ
 * @ingroup THDM 
 * @brief lambda1atQ.
 */
class lambda1atQ: public Runner {
public:

    /**
     * @brief lambda1atQ constructor.
     */
    lambda1atQ(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_1(Q)@f$
     */
    double computeThValue();

private:
//    const THDM * myTHDM;
//    Runner * myRunner;
//    m11_2 * mym11_2;
//    m22_2 * mym22_2;
//    lambda1 * mylambda1;
//    lambda2 * mylambda2;
//    lambda3 * mylambda3;
//    lambda4 * mylambda4;
//    lambda5 * mylambda5;
};

/**
 * @class lambda2atQ
 * @ingroup THDM 
 * @brief lambda2atQ.
 */
class lambda2atQ: public ThObservable {
public:

    /**
     * @brief lambda2atQ constructor.
     */
    lambda2atQ(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_2(Q)@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class lambda3atQ
 * @ingroup THDM 
 * @brief lambda3atQ.
 */
class lambda3atQ: public ThObservable {
public:

    /**
     * @brief lambda3atQ constructor.
     */
    lambda3atQ(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_3(Q)@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class lambda4atQ
 * @ingroup THDM 
 * @brief lambda4atQ.
 */
class lambda4atQ: public ThObservable {
public:

    /**
     * @brief lambda4atQ constructor.
     */
    lambda4atQ(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_4(Q)@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class lambda5atQ
 * @ingroup THDM 
 * @brief lambda5atQ.
 */
class lambda5atQ: public ThObservable {
public:

    /**
     * @brief lambda5atQ constructor.
     */
    lambda5atQ(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_5(Q)@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};



#endif	/* RUNNER_H */
