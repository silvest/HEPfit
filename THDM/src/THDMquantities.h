/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMQUANTITIES_H
#define	THDMQUANTITIES_H

#include <stdexcept>
#include "ThObservable.h"
#include "THDM.h"

/**
 * @class mass_mHh
 * @ingroup THDM 
 * @brief The mass of the heavy CP-even Higgs state.
 */
class mass_mHh: public ThObservable {
public:

    /**
     * @brief mass_mHh constructor.
     */
    mass_mHh(const StandardModel& SM_i);

    /**
     * @return @f$m_H@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class mass_mA
 * @ingroup THDM 
 * @brief The mass of the CP-odd Higgs state.
 */
class mass_mA: public ThObservable {
public:

    /**
     * @brief mass_mA constructor.
     */
    mass_mA(const StandardModel& SM_i);

    /**
     * @return @f$m_A@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class mass_mHp
 * @ingroup THDM 
 * @brief The mass of the charged Higgs state.
 */
class mass_mHp: public ThObservable {
public:

    /**
     * @brief mass_mHp constructor.
     */
    mass_mHp(const StandardModel& SM_i);

    /**
     * @return @f$m_{H^+}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class masssquare_mA
 * @ingroup THDM 
 * @brief The squared mass of the CP-odd Higgs state.
 */
class masssquare_mA: public ThObservable {
public:

    /**
     * @brief masssquare_mA constructor.
     */
    masssquare_mA(const StandardModel& SM_i);

    /**
     * @return @f$m_A^2@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class masssquare_mHp
 * @ingroup THDM 
 * @brief The squared mass of the charged Higgs state.
 */
class masssquare_mHp: public ThObservable {
public:

    /**
     * @brief masssquare_mHp constructor.
     */
    masssquare_mHp(const StandardModel& SM_i);

    /**
     * @return @f$m_{H^+}^2@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class massdifference_mHhmmA
 * @ingroup THDM 
 * @brief The difference between the masses of the heavy CP-even and the CP-odd Higgs.
 */
class massdifference_mHhmmA: public ThObservable {
public:

    /**
     * @brief massdifference_mHhmmA constructor.
     */
    massdifference_mHhmmA(const StandardModel& SM_i);

    /**
     * @return @f$m_H-m_A@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class massdifference_mAmmHh
 * @ingroup THDM 
 * @brief The difference between the masses of the CP-odd and the heavy CP-even Higgs.
 */
class massdifference_mAmmHh: public ThObservable {
public:

    /**
     * @brief massdifference_mAmmHh constructor.
     */
    massdifference_mAmmHh(const StandardModel& SM_i);

    /**
     * @return @f$m_A-m_H@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class massdifference_mHhmmHp
 * @ingroup THDM 
 * @brief The difference between the masses of the heavy CP-even and the charged Higgs.
 */
class massdifference_mHhmmHp: public ThObservable {
public:

    /**
     * @brief massdifference_mHhmmHp constructor.
     */
    massdifference_mHhmmHp(const StandardModel& SM_i);

    /**
     * @return @f$m_H-m_{H^+}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class massdifference_mHpmmHh
 * @ingroup THDM 
 * @brief The difference between the masses of the charged and the heavy CP-even Higgs.
 */
class massdifference_mHpmmHh: public ThObservable {
public:

    /**
     * @brief massdifference_mHpmmHh constructor.
     */
    massdifference_mHpmmHh(const StandardModel& SM_i);

    /**
     * @return @f$m_{H^+}-m_H@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class massdifference_mAmmHp
 * @ingroup THDM 
 * @brief The difference between the masses of the CP-odd and the charged Higgs.
 */
class massdifference_mAmmHp: public ThObservable {
public:

    /**
     * @brief massdifference_mAmmHp constructor.
     */
    massdifference_mAmmHp(const StandardModel& SM_i);

    /**
     * @return @f$m_A-m_{H^+}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class massdifference_mHpmmA
 * @ingroup THDM 
 * @brief The difference between the masses of the charged and the CP-odd Higgs.
 */
class massdifference_mHpmmA: public ThObservable {
public:

    /**
     * @brief massdifference_mHpmmA constructor.
     */
    massdifference_mHpmmA(const StandardModel& SM_i);

    /**
     * @return @f$m_{H^+}-m_A@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class lambda1 
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_1@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_1@f$.
 */
class lambda1 : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    lambda1(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_1@f$.
     * @return @f$\lambda_1@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class lambda2
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_2@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_2@f$.
 */
class lambda2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    lambda2(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_2@f$.
     * @return @f$\lambda_2@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class lambda3
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_3@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_3@f$.
 */
class lambda3 : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    lambda3(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_3@f$.
     * @return @f$\lambda_3@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class lambda4
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_4@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_4@f$.
 */
class lambda4 : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    lambda4(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_4@f$.
     * @return @f$\lambda_4@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class lambda5
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_5@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_5@f$.
 */
class lambda5 : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    lambda5(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_5@f$.
     * @return @f$\lambda_5@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

#endif	/* THDMQUANTITIES_H */
