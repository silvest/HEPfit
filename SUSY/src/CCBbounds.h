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
 * @class CCBu11
 * @ingroup LeptonFlavour
 * @brief CCBu11.
 */
class CCBu11: public ThObservable {
public:

    /**
     * @brief CCBu11 constructor.
     */
    CCBu11(const StandardModel& SM_i);

    /**
     * @return @CCBu11@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu22
 * @ingroup LeptonFlavour
 * @brief CCBu22.
 */
class CCBu22: public ThObservable {
public:

    /**
     * @brief CCBu22 constructor.
     */
    CCBu22(const StandardModel& SM_i);

    /**
     * @return @CCBu22@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu33
 * @ingroup LeptonFlavour
 * @brief CCBu33.
 */
class CCBu33: public ThObservable {
public:

    /**
     * @brief CCBu33 constructor.
     */
    CCBu33(const StandardModel& SM_i);

    /**
     * @return @CCBu33@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu12
 * @ingroup LeptonFlavour
 * @brief CCBu12.
 */
class CCBu12: public ThObservable {
public:

    /**
     * @brief CCBu12 constructor.
     */
    CCBu12(const StandardModel& SM_i);

    /**
     * @return @CCBu12@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu13
 * @ingroup LeptonFlavour
 * @brief CCBu13.
 */
class CCBu13: public ThObservable {
public:

    /**
     * @brief CCBu13 constructor.
     */
    CCBu13(const StandardModel& SM_i);

    /**
     * @return @CCBu13@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBu23
 * @ingroup LeptonFlavour
 * @brief CCBu23.
 */
class CCBu23: public ThObservable {
public:

    /**
     * @brief CCBu23 constructor.
     */
    CCBu23(const StandardModel& SM_i);

    /**
     * @return @CCBu23@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd11
 * @ingroup LeptonFlavour
 * @brief CCBd11.
 */
class CCBd11: public ThObservable {
public:

    /**
     * @brief CCBd11 constructor.
     */
    CCBd11(const StandardModel& SM_i);

    /**
     * @return @CCBd11@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd22
 * @ingroup LeptonFlavour
 * @brief CCBd22.
 */
class CCBd22: public ThObservable {
public:

    /**
     * @brief CCBd22 constructor.
     */
    CCBd22(const StandardModel& SM_i);

    /**
     * @return @CCBd22@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd33
 * @ingroup LeptonFlavour
 * @brief CCBd33.
 */
class CCBd33: public ThObservable {
public:

    /**
     * @brief CCBd33 constructor.
     */
    CCBd33(const StandardModel& SM_i);

    /**
     * @return @CCBd33@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd12
 * @ingroup LeptonFlavour
 * @brief CCBd12.
 */
class CCBd12: public ThObservable {
public:

    /**
     * @brief CCBd12 constructor.
     */
    CCBd12(const StandardModel& SM_i);

    /**
     * @return @CCBd12@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd13
 * @ingroup LeptonFlavour
 * @brief CCBd13.
 */
class CCBd13: public ThObservable {
public:

    /**
     * @brief CCBd13 constructor.
     */
    CCBd13(const StandardModel& SM_i);

    /**
     * @return @CCBd13@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBd23
 * @ingroup LeptonFlavour
 * @brief CCBd23.
 */
class CCBd23: public ThObservable {
public:

    /**
     * @brief CCBd23 constructor.
     */
    CCBd23(const StandardModel& SM_i);

    /**
     * @return @CCBd23@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe11
 * @ingroup LeptonFlavour
 * @brief CCBe11.
 */
class CCBe11: public ThObservable {
public:

    /**
     * @brief CCBe11 constructor.
     */
    CCBe11(const StandardModel& SM_i);

    /**
     * @return @CCBe11@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe22
 * @ingroup LeptonFlavour
 * @brief CCBe22.
 */
class CCBe22: public ThObservable {
public:

    /**
     * @brief CCBe22 constructor.
     */
    CCBe22(const StandardModel& SM_i);

    /**
     * @return @CCBe22@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe33
 * @ingroup LeptonFlavour
 * @brief CCBe33.
 */
class CCBe33: public ThObservable {
public:

    /**
     * @brief CCBe33 constructor.
     */
    CCBe33(const StandardModel& SM_i);

    /**
     * @return @CCBe33@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe12
 * @ingroup LeptonFlavour
 * @brief CCBe12.
 */
class CCBe12: public ThObservable {
public:

    /**
     * @brief CCBe12 constructor.
     */
    CCBe12(const StandardModel& SM_i);

    /**
     * @return @CCBe12@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe13
 * @ingroup LeptonFlavour
 * @brief CCBe13.
 */
class CCBe13: public ThObservable {
public:

    /**
     * @brief CCBe13 constructor.
     */
    CCBe13(const StandardModel& SM_i);

    /**
     * @return @CCBe13@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class CCBe23
 * @ingroup LeptonFlavour
 * @brief CCBe23.
 */
class CCBe23: public ThObservable {
public:

    /**
     * @brief CCBe23 constructor.
     */
    CCBe23(const StandardModel& SM_i);

    /**
     * @return @CCBe23@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

#endif	/* CCBBOUNDS_H */
