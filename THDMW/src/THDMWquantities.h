/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWQUANTITIES_H
#define	THDMWQUANTITIES_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDMW.h"
#include "THDMWcache.h"

/**
 * @class m12sqTHDMW
 * @ingroup THDMW
 * @brief m12sqTHDMW.
 */
class m12sqTHDMW: public ThObservable {
public:

    /**
     * @brief m12sqTHDMW constructor.
     */
    m12sqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{12}^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class m11sqTHDMW
 * @ingroup THDMW
 * @brief m11sqTHDMW.
 */
class m11sqTHDMW: public ThObservable {
public:

    /**
     * @brief m11sqTHDMW constructor.
     */
    m11sqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{11}^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class m22sqTHDMW
 * @ingroup THDMW
 * @brief m22sqTHDMW.
 */
class m22sqTHDMW: public ThObservable {
public:

    /**
     * @brief m22sqTHDMW constructor.
     */
    m22sqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{22}^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mhsqTHDMW
 * @ingroup THDMW
 * @brief mhsqTHDMW.
 */
class mhsqTHDMW: public ThObservable {
public:

    /**
     * @brief mhsqTHDMW constructor.
     */
    mhsqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_h^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mHHsqTHDMW
 * @ingroup THDMW
 * @brief mHHsqTHDMW.
 */
class mHHsqTHDMW: public ThObservable {
public:

    /**
     * @brief mHHsqTHDMW constructor.
     */
    mHHsqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_H^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mAsqTHDMW
 * @ingroup THDMW
 * @brief mAsqTHDMW.
 */
class mAsqTHDMW: public ThObservable {
public:

    /**
     * @brief mAsqTHDMW constructor.
     */
    mAsqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_A^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mSRsqTHDMW
 * @ingroup THDMW
 * @brief mSRsqTHDMW.
 */
class mSRsqTHDMW: public ThObservable {
public:

    /**
     * @brief mSRsqTHDMW constructor.
     */
    mSRsqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$Re[m_S]^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mSIsqTHDMW
 * @ingroup THDMW
 * @brief mSIsqTHDMW.
 */
class mSIsqTHDMW: public ThObservable {
public:

    /**
     * @brief mSIsqTHDMW constructor.
     */
    mSIsqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$Im[m_S]^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mHpsqTHDMW
 * @ingroup THDMW
 * @brief mHpsqTHDMW.
 */
class mHpsqTHDMW: public ThObservable {
public:

    /**
     * @brief mHpsqTHDMW constructor.
     */
    mHpsqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{H^+}^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mSpsqTHDMW
 * @ingroup THDMW
 * @brief mSpsqTHDMW.
 */
class mSpsqTHDMW: public ThObservable {
public:

    /**
     * @brief mSpsqTHDMW constructor.
     */
    mSpsqTHDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{S^+}^2}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

#endif	/* THDMWQUANTITIES_H */
