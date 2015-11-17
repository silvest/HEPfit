/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMQUANTITIES_H
#define	THDMQUANTITIES_H

#include <stdexcept>
#include <ThObservable.h>
#include "THDM.h"

/**
 * @class THDMquantities
 * @brief .
 */
class THDMquantities : public ThObservable {
public:
    THDMquantities(const StandardModel& SM_i);
    virtual ~THDMquantities();
    double computeThValue();
    const THDM * myTHDM;
    const StandardModel& mySM;
};

class mass_mHh: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    mass_mHh(const StandardModel& SM_i);

    /**
     * @return mHh
     */
    double computeThValue();
};

class mass_mA: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    mass_mA(const StandardModel& SM_i);

    /**
     * @return mA
     */
    double computeThValue();
};

class mass_mHp: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    mass_mHp(const StandardModel& SM_i);

    /**
     * @return mHp
     */
    double computeThValue();
};

class masssquare_mA: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    masssquare_mA(const StandardModel& SM_i);

    /**
     * @return mA2
     */
    double computeThValue();
};

class masssquare_mHp: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    masssquare_mHp(const StandardModel& SM_i);

    /**
     * @return mHp2
     */
    double computeThValue();
};

class massdifference_mHhmmA: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    massdifference_mHhmmA(const StandardModel& SM_i);

    /**
     * @return mHh-mA
     */
    double computeThValue();
};

class massdifference_mAmmHh: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    massdifference_mAmmHh(const StandardModel& SM_i);

    /**
     * @return mA-mHh
     */
    double computeThValue();
};

class massdifference_mHhmmHp: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    massdifference_mHhmmHp(const StandardModel& SM_i);

    /**
     * @return mHh-mHp
     */
    double computeThValue();
};

class massdifference_mHpmmHh: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    massdifference_mHpmmHh(const StandardModel& SM_i);

    /**
     * @return mHp-mHh
     */
    double computeThValue();
};

class massdifference_mAmmHp: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    massdifference_mAmmHp(const StandardModel& SM_i);

    /**
     * @return mA-mHp
     */
    double computeThValue();
};

class massdifference_mHpmmA: public THDMquantities {
public:

    /**
     * @brief Constructor.
     */
    massdifference_mHpmmA(const StandardModel& SM_i);

    /**
     * @return mHp-mA
     */
    double computeThValue();
};

/**
 * @class lambda1 
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_1@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_1@f$.
 */
class lambda1 : public THDMquantities {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    lambda1(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_1@f$.
     * @return @f$\lambda_1@f$
     */
    double computeThValue();
};

/**
 * @class lambda2
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_2@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_2@f$.
 */
class lambda2 : public THDMquantities {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    lambda2(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_2@f$.
     * @return @f$\lambda_2@f$
     */
    double computeThValue();
};

/**
 * @class lambda3
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_3@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_3@f$.
 */
class lambda3 : public THDMquantities {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    lambda3(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_3@f$.
     * @return @f$\lambda_3@f$
     */
    double computeThValue();
};

/**
 * @class lambda4
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_4@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_4@f$.
 */
class lambda4 : public THDMquantities {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    lambda4(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_4@f$.
     * @return @f$\lambda_4@f$
     */
    double computeThValue();
};

/**
 * @class lambda5
 * @ingroup THDM 
 * @brief An observable class for the quartic Higgs potential coupling @f$\lambda_5@f$.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the quartic Higgs potential coupling @f$\lambda_5@f$.
 */
class lambda5 : public THDMquantities {
public:

    /**
     * @brief Constructor.
     * @param[in] ?
     */
    lambda5(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$\lambda_5@f$.
     * @return @f$\lambda_5@f$
     */
    double computeThValue();
};

#endif	/* THDMQUANTITIES_H */