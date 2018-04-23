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

/**
 * @class mAmmHH_THDMW
 * @ingroup THDMW
 * @brief mAmmHH_THDMW.
 */
class mAmmHH_THDMW: public ThObservable {
public:

    /**
     * @brief mAmmHH_THDMW constructor.
     */
    mAmmHH_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_A-m_H}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mHHmmA_THDMW
 * @ingroup THDMW
 * @brief mHHmmA_THDMW.
 */
class mHHmmA_THDMW: public ThObservable {
public:

    /**
     * @brief mHHmmA_THDMW constructor.
     */
    mHHmmA_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_H-m_A}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mAmmSR_THDMW
 * @ingroup THDMW
 * @brief mAmmSR_THDMW.
 */
class mAmmSR_THDMW: public ThObservable {
public:

    /**
     * @brief mAmmSR_THDMW constructor.
     */
    mAmmSR_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_A-m_{S_R}}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mSRmmA_THDMW
 * @ingroup THDMW
 * @brief mSRmmA_THDMW.
 */
class mSRmmA_THDMW: public ThObservable {
public:

    /**
     * @brief mSRmmA_THDMW constructor.
     */
    mSRmmA_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{S_R}-m_A}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mAmmSI_THDMW
 * @ingroup THDMW
 * @brief mAmmSI_THDMW.
 */
class mAmmSI_THDMW: public ThObservable {
public:

    /**
     * @brief mAmmSI_THDMW constructor.
     */
    mAmmSI_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_A-m_{S_I}}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mSImmA_THDMW
 * @ingroup THDMW
 * @brief mSImmA_THDMW.
 */
class mSImmA_THDMW: public ThObservable {
public:

    /**
     * @brief mSImmA_THDMW constructor.
     */
    mSImmA_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{S_I}-m_A}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mHHmmSR_THDMW
 * @ingroup THDMW
 * @brief mHHmmSR_THDMW.
 */
class mHHmmSR_THDMW: public ThObservable {
public:

    /**
     * @brief mHHmmSR_THDMW constructor.
     */
    mHHmmSR_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_H-m_{S_R}}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mSRmmHH_THDMW
 * @ingroup THDMW
 * @brief mSRmmHH_THDMW.
 */
class mSRmmHH_THDMW: public ThObservable {
public:

    /**
     * @brief mSRmmHH_THDMW constructor.
     */
    mSRmmHH_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{S_R}-m_H}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mHHmmSI_THDMW
 * @ingroup THDMW
 * @brief mHHmmSI_THDMW.
 */
class mHHmmSI_THDMW: public ThObservable {
public:

    /**
     * @brief mHHmmSI_THDMW constructor.
     */
    mHHmmSI_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_H-m_{S_I}}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mSImmHH_THDMW
 * @ingroup THDMW
 * @brief mSImmHH_THDMW.
 */
class mSImmHH_THDMW: public ThObservable {
public:

    /**
     * @brief mSImmHH_THDMW constructor.
     */
    mSImmHH_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{S_I}-m_H}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mSRmmSI_THDMW
 * @ingroup THDMW
 * @brief mSRmmSI_THDMW.
 */
class mSRmmSI_THDMW: public ThObservable {
public:

    /**
     * @brief mSRmmSI_THDMW constructor.
     */
    mSRmmSI_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{S_R}-m_{S_I}}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

/**
 * @class mSImmSR_THDMW
 * @ingroup THDMW
 * @brief mSImmSR_THDMW.
 */
class mSImmSR_THDMW: public ThObservable {
public:

    /**
     * @brief mSImmSR_THDMW constructor.
     */
    mSImmSR_THDMW(const StandardModel& SM_i);

    /**
     * @return @f$m_{S_I}-m_{S_R}}@f$
     */
    double computeThValue();

private:
    const THDMW& myTHDMW;
};

#endif	/* THDMWQUANTITIES_H */
