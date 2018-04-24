/*
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWPOSITIVITY_H
#define	THDMWPOSITIVITY_H

#include "ThObservable.h"
#include "THDMW.h"
#include "THDMWcache.h"

/**
 * @class THDMWpositivity
 * @ingroup THDMW 
 * @brief Base class for "boundedness-from-below" constraints
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details Formulae taken from @cite Deshpande:1977rw.
 */

/**
 * @class THDMWpositivity1
 * @brief Controls that the scalar %THDMW potential is bounded from below.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class THDMWpositivity1: public ThObservable {
public:

    /**
     * @brief THDMWpositivity1 constructor.
     */
    THDMWpositivity1(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity2
 * @brief Controls that the scalar %THDMW potential is bounded from below.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class THDMWpositivity2: public ThObservable {
public:

    /**
     * @brief THDMWpositivity2 constructor.
     */
    THDMWpositivity2(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity3
 */
class THDMWpositivity3: public ThObservable {
public:

    /**
     * @brief THDMWpositivity3 constructor.
     */
    THDMWpositivity3(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity4
 */
class THDMWpositivity4: public ThObservable {
public:

    /**
     * @brief THDMWpositivity4 constructor.
     */
    THDMWpositivity4(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity5
 */
class THDMWpositivity5: public ThObservable {
public:

    /**
     * @brief THDMWpositivity5 constructor.
     */
    THDMWpositivity5(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity6
 */
class THDMWpositivity6: public ThObservable {
public:

    /**
     * @brief THDMWpositivity6 constructor.
     */
    THDMWpositivity6(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity7
 */
class THDMWpositivity7: public ThObservable {
public:

    /**
     * @brief THDMWpositivity7 constructor.
     */
    THDMWpositivity7(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity8
 */
class THDMWpositivity8: public ThObservable {
public:

    /**
     * @brief THDMWpositivity8 constructor.
     */
    THDMWpositivity8(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity9
 */
class THDMWpositivity9: public ThObservable {
public:

    /**
     * @brief THDMWpositivity9 constructor.
     */
    THDMWpositivity9(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity10
 */
class THDMWpositivity10: public ThObservable {
public:

    /**
     * @brief THDMWpositivity10 constructor.
     */
    THDMWpositivity10(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity11
 */
class THDMWpositivity11: public ThObservable {
public:

    /**
     * @brief THDMWpositivity11 constructor.
     */
    THDMWpositivity11(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositivity12
 */
class THDMWpositivity12: public ThObservable {
public:

    /**
     * @brief THDMWpositivity12 constructor.
     */
    THDMWpositivity12(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

/**
 * @class THDMWpositiveMassSquares
 */
class THDMWpositiveMassSquares: public ThObservable {
public:

    /**
     * @brief THDMWpositiveMassSquares constructor.
     */
    THDMWpositiveMassSquares(const StandardModel& SM_i);

    /**
     * @return 
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
};

#endif	/* THDMWPOSITIVITY_H */
