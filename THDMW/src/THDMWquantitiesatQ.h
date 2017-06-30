/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWQUANTITIESATQ_H
#define	THDMWQUANTITIESATQ_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDMW.h"
#include "THDMWcache.h"

/**
 * @class Q_stTHDMW
 * @ingroup THDMW
 * @brief Q_stTHDMW.
 */
class Q_stTHDMW: public ThObservable {
public:

    /**
     * @brief Q_stTHDMW constructor.
     */
    Q_stTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$Q_{\text{stability}}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class DeltaQ_THDMW
 * @ingroup THDMW
 * @brief DeltaQ_THDMW.
 */
class DeltaQ_THDMW: public ThObservable {
public:

    /**
     * @brief DeltaQ_THDMW constructor.
     */
    DeltaQ_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$Q_{\text{THDMW}}-Q_{\text{stability}}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class lambda1atQTHDMW
 * @ingroup THDMW
 * @brief lambda1atQTHDMW.
 */
class lambda1atQTHDMW: public ThObservable {
public:

    /**
     * @brief lambda1atQTHDMW constructor.
     */
    lambda1atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_1(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class lambda2atQTHDMW
 * @ingroup THDMW
 * @brief lambda2atQTHDMW.
 */
class lambda2atQTHDMW: public ThObservable {
public:

    /**
     * @brief lambda2atQTHDMW constructor.
     */
    lambda2atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_2(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class lambda3atQTHDMW
 * @ingroup THDMW
 * @brief lambda3atQTHDMW.
 */
class lambda3atQTHDMW: public ThObservable {
public:

    /**
     * @brief lambda3atQTHDMW constructor.
     */
    lambda3atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_3(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class lambda4atQTHDMW
 * @ingroup THDMW
 * @brief lambda4atQTHDMW.
 */
class lambda4atQTHDMW: public ThObservable {
public:

    /**
     * @brief lambda4atQTHDMW constructor.
     */
    lambda4atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_4(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mu1atQTHDMW
 * @ingroup THDMW
 * @brief mu1atQTHDMW.
 */
class mu1atQTHDMW: public ThObservable {
public:

    /**
     * @brief mu1atQTHDMW constructor.
     */
    mu1atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\mu_1(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mu3atQTHDMW
 * @ingroup THDMW
 * @brief mu3atQTHDMW.
 */
class mu3atQTHDMW: public ThObservable {
public:

    /**
     * @brief mu3atQTHDMW constructor.
     */
    mu3atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\mu_3(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mu4atQTHDMW
 * @ingroup THDMW
 * @brief mu4atQTHDMW.
 */
class mu4atQTHDMW: public ThObservable {
public:

    /**
     * @brief mu4atQTHDMW constructor.
     */
    mu4atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\mu_4(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class nu1atQTHDMW
 * @ingroup THDMW
 * @brief nu1atQTHDMW.
 */
class nu1atQTHDMW: public ThObservable {
public:

    /**
     * @brief nu1atQTHDMW constructor.
     */
    nu1atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\nu_1(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class omega1atQTHDMW
 * @ingroup THDMW
 * @brief omega1atQTHDMW.
 */
class omega1atQTHDMW: public ThObservable {
public:

    /**
     * @brief omega1atQTHDMW constructor.
     */
    omega1atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\omega_1(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class kappa1atQTHDMW
 * @ingroup THDMW
 * @brief kappa1atQTHDMW.
 */
class kappa1atQTHDMW: public ThObservable {
public:

    /**
     * @brief kappa1atQTHDMW constructor.
     */
    kappa1atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\kappa_1(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class nu2atQTHDMW
 * @ingroup THDMW
 * @brief nu2atQTHDMW.
 */
class nu2atQTHDMW: public ThObservable {
public:

    /**
     * @brief nu2atQTHDMW constructor.
     */
    nu2atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\nu_2(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class omega2atQTHDMW
 * @ingroup THDMW
 * @brief omega2atQTHDMW.
 */
class omega2atQTHDMW: public ThObservable {
public:

    /**
     * @brief omega2atQTHDMW constructor.
     */
    omega2atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\omega_2(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class kappa2atQTHDMW
 * @ingroup THDMW
 * @brief kappa2atQTHDMW.
 */
class kappa2atQTHDMW: public ThObservable {
public:

    /**
     * @brief kappa2atQTHDMW constructor.
     */
    kappa2atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\kappa_2(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class nu4atQTHDMW
 * @ingroup THDMW
 * @brief nu4atQTHDMW.
 */
class nu4atQTHDMW: public ThObservable {
public:

    /**
     * @brief nu4atQTHDMW constructor.
     */
    nu4atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\nu_4(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class omega4atQTHDMW
 * @ingroup THDMW
 * @brief omega4atQTHDMW.
 */
class omega4atQTHDMW: public ThObservable {
public:

    /**
     * @brief omega4atQTHDMW constructor.
     */
    omega4atQTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$\omega_4(Q)@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

#endif	/* THDMWQUANTITIESATQ_H */
