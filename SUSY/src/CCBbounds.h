/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CCBBOUNDS_H
#define	CCBBOUNDS_H

//#include <stdexcept>
#include "ThObservable.h"
#include "SUSY.h"


/**
 * @addtogroup SUSY
 * @brief A module for a basis of %SUSY models.
 * @{
 */

/**
 * @class CCBu11
 * @ingroup SUSY
 * @brief CCBu11.
 */
class CCBu11: public ThObservable {
public:

    /**
     * @brief CCBu11 constructor.
     */
    CCBu11(const StandardModel& SM_i);

    /**
     * @return @f$CCBu11@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu22
 * @ingroup SUSY
 * @brief CCBu22.
 */
class CCBu22: public ThObservable {
public:

    /**
     * @brief CCBu22 constructor.
     */
    CCBu22(const StandardModel& SM_i);

    /**
     * @return @f$CCBu22@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu33
 * @ingroup SUSY
 * @brief CCBu33.
 */
class CCBu33: public ThObservable {
public:

    /**
     * @brief CCBu33 constructor.
     */
    CCBu33(const StandardModel& SM_i);

    /**
     * @return @f$CCBu33@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu12
 * @ingroup SUSY
 * @brief CCBu12.
 */
class CCBu12: public ThObservable {
public:

    /**
     * @brief CCBu12 constructor.
     */
    CCBu12(const StandardModel& SM_i);

    /**
     * @return @f$CCBu12@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu13
 * @ingroup SUSY
 * @brief CCBu13.
 */
class CCBu13: public ThObservable {
public:

    /**
     * @brief CCBu13 constructor.
     */
    CCBu13(const StandardModel& SM_i);

    /**
     * @return @f$CCBu13@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu23
 * @ingroup SUSY
 * @brief CCBu23.
 */
class CCBu23: public ThObservable {
public:

    /**
     * @brief CCBu23 constructor.
     */
    CCBu23(const StandardModel& SM_i);

    /**
     * @return @f$CCBu23@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd11
 * @ingroup SUSY
 * @brief CCBd11.
 */
class CCBd11: public ThObservable {
public:

    /**
     * @brief CCBd11 constructor.
     */
    CCBd11(const StandardModel& SM_i);

    /**
     * @return @f$CCBd11@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd22
 * @ingroup SUSY
 * @brief CCBd22.
 */
class CCBd22: public ThObservable {
public:

    /**
     * @brief CCBd22 constructor.
     */
    CCBd22(const StandardModel& SM_i);

    /**
     * @return @f$CCBd22@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd33
 * @ingroup SUSY
 * @brief CCBd33.
 */
class CCBd33: public ThObservable {
public:

    /**
     * @brief CCBd33 constructor.
     */
    CCBd33(const StandardModel& SM_i);

    /**
     * @return @f$CCBd33@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd12
 * @ingroup SUSY
 * @brief CCBd12.
 */
class CCBd12: public ThObservable {
public:

    /**
     * @brief CCBd12 constructor.
     */
    CCBd12(const StandardModel& SM_i);

    /**
     * @return @f$CCBd12@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd13
 * @ingroup SUSY
 * @brief CCBd13.
 */
class CCBd13: public ThObservable {
public:

    /**
     * @brief CCBd13 constructor.
     */
    CCBd13(const StandardModel& SM_i);

    /**
     * @return @f$CCBd13@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd23
 * @ingroup SUSY
 * @brief CCBd23.
 */
class CCBd23: public ThObservable {
public:

    /**
     * @brief CCBd23 constructor.
     */
    CCBd23(const StandardModel& SM_i);

    /**
     * @return @f$CCBd23@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe11
 * @ingroup SUSY
 * @brief CCBe11.
 */
class CCBe11: public ThObservable {
public:

    /**
     * @brief CCBe11 constructor.
     */
    CCBe11(const StandardModel& SM_i);

    /**
     * @return @f$CCBe11@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe22
 * @ingroup SUSY
 * @brief CCBe22.
 */
class CCBe22: public ThObservable {
public:

    /**
     * @brief CCBe22 constructor.
     */
    CCBe22(const StandardModel& SM_i);

    /**
     * @return @f$CCBe22@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe33
 * @ingroup SUSY
 * @brief CCBe33.
 */
class CCBe33: public ThObservable {
public:

    /**
     * @brief CCBe33 constructor.
     */
    CCBe33(const StandardModel& SM_i);

    /**
     * @return @f$CCBe33@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe12
 * @ingroup SUSY
 * @brief CCBe12.
 */
class CCBe12: public ThObservable {
public:

    /**
     * @brief CCBe12 constructor.
     */
    CCBe12(const StandardModel& SM_i);

    /**
     * @return @f$CCBe12@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe13
 * @ingroup SUSY
 * @brief CCBe13.
 */
class CCBe13: public ThObservable {
public:

    /**
     * @brief CCBe13 constructor.
     */
    CCBe13(const StandardModel& SM_i);

    /**
     * @return @f$CCBe13@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe23
 * @ingroup SUSY
 * @brief CCBe23.
 */
class CCBe23: public ThObservable {
public:

    /**
     * @brief CCBe23 constructor.
     */
    CCBe23(const StandardModel& SM_i);

    /**
     * @return @f$CCBe23@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @}
 */

#endif	/* CCBBOUNDS_H */
