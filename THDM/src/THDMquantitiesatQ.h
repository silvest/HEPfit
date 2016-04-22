/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMQUANTITIESATQ_H
#define	THDMQUANTITIESATQ_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"
#include "THDMcache.h"

/**
 * @class Q_st
 * @ingroup THDM
 * @brief Q_st.
 */
class Q_st: public ThObservable {
public:

    /**
     * @brief Q_st constructor.
     */
    Q_st(const StandardModel& SM_i);

    /**
     * @return @f$Q_{\text{stability}}@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class DeltaQ_THDM
 * @ingroup THDM
 * @brief DeltaQ_THDM.
 */
class DeltaQ_THDM: public ThObservable {
public:

    /**
     * @brief DeltaQ_THDM constructor.
     */
    DeltaQ_THDM(const StandardModel& SM_i);

    /**
     * @return @f$Q_{\text{THDM}}-Q_{\text{stability}}@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class g1atQ
 * @ingroup THDM
 * @brief g1atQ.
 */
class g1atQ: public ThObservable {
public:

    /**
     * @brief g1atQ constructor.
     */
    g1atQ(const StandardModel& SM_i);

    /**
     * @return @f$g_1(Q)@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class g2atQ
 * @ingroup THDM
 * @brief g2atQ.
 */
class g2atQ: public ThObservable {
public:

    /**
     * @brief g2atQ constructor.
     */
    g2atQ(const StandardModel& SM_i);

    /**
     * @return @f$g_2(Q)@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class g3atQ
 * @ingroup THDM
 * @brief g3atQ.
 */
class g3atQ: public ThObservable {
public:

    /**
     * @brief g3atQ constructor.
     */
    g3atQ(const StandardModel& SM_i);

    /**
     * @return @f$g_3(Q)@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class YtopatQ
 * @ingroup THDM
 * @brief YtopatQ.
 */
class YtopatQ: public ThObservable {
public:

    /**
     * @brief YtopatQ constructor.
     */
    YtopatQ(const StandardModel& SM_i);

    /**
     * @return @f$Y_t(Q)@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class YbottomatQ
 * @ingroup THDM
 * @brief YbottomatQ.
 */
class YbottomatQ: public ThObservable {
public:

    /**
     * @brief YbottomatQ constructor.
     */
    YbottomatQ(const StandardModel& SM_i);

    /**
     * @return @f$Y_b(Q)@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class YtauatQ
 * @ingroup THDM
 * @brief YtauatQ.
 */
class YtauatQ: public ThObservable {
public:

    /**
     * @brief YtauatQ constructor.
     */
    YtauatQ(const StandardModel& SM_i);

    /**
     * @return @f$Y_\tau(Q)@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class m11_2atQ
 * @ingroup THDM
 * @brief m11_2atQ.
 */
class m11_2atQ: public ThObservable {
public:

    /**
     * @brief m11_2atQ constructor.
     */
    m11_2atQ(const StandardModel& SM_i);

    /**
     * @return @f$m_{11}^2(Q)@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class m22_2atQ
 * @ingroup THDM
 * @brief m22_2atQ.
 */
class m22_2atQ: public ThObservable {
public:

    /**
     * @brief m22_2atQ constructor.
     */
    m22_2atQ(const StandardModel& SM_i);

    /**
     * @return @f$m_{22}^2(Q)@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class m12_2atQ
 * @ingroup THDM
 * @brief m12_2atQ.
 */
class m12_2atQ: public ThObservable {
public:

    /**
     * @brief m12_2atQ constructor.
     */
    m12_2atQ(const StandardModel& SM_i);

    /**
     * @return @f$m_{12}^2(Q)@f$
     */
    double computeThValue();

private:
    const THDM& myTHDM;
};

/**
 * @class lambda1atQ
 * @ingroup THDM
 * @brief lambda1atQ.
 */
class lambda1atQ: public ThObservable {
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
    const THDM& myTHDM;
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

private:
    const THDM& myTHDM;
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

private:
    const THDM& myTHDM;
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

private:
    const THDM& myTHDM;
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

private:
    const THDM& myTHDM;
};

#endif	/* THDMQUANTITIESATQ_H */
