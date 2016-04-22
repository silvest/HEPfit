/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef SFERMIONDELTAS_H
#define	SFERMIONDELTAS_H

//#include <stdexcept>
#include "ThObservable.h"
#include "SUSY.h"

/**
 * @class deltaLL1_q
 * @ingroup SUSY
 * @brief deltaLL1_q.
 */
class deltaLL1_q: public ThObservable {
public:

    /**
     * @brief deltaLL1_q constructor.
     */
    deltaLL1_q(const StandardModel& SM_i);

    /**
     * @return @f$\delta LL1_q@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaLL2_q
 * @ingroup SUSY
 * @brief deltaLL2_q.
 */
class deltaLL2_q: public ThObservable {
public:

    /**
     * @brief deltaLL2_q constructor.
     */
    deltaLL2_q(const StandardModel& SM_i);

    /**
     * @return @f$\delta LL2_q@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaLL3_q
 * @ingroup SUSY
 * @brief deltaLL3_q.
 */
class deltaLL3_q: public ThObservable {
public:

    /**
     * @brief deltaLL3_q constructor.
     */
    deltaLL3_q(const StandardModel& SM_i);

    /**
     * @return @f$\delta LL3_q@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRR1_u
 * @ingroup SUSY
 * @brief deltaRR1_u.
 */
class deltaRR1_u: public ThObservable {
public:

    /**
     * @brief deltaRR1_u constructor.
     */
    deltaRR1_u(const StandardModel& SM_i);

    /**
     * @return @f$\delta RR1_u@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRR2_u
 * @ingroup SUSY
 * @brief deltaRR2_u.
 */
class deltaRR2_u: public ThObservable {
public:

    /**
     * @brief deltaRR2_u constructor.
     */
    deltaRR2_u(const StandardModel& SM_i);

    /**
     * @return @f$\delta RR2_u@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRR3_u
 * @ingroup SUSY
 * @brief deltaRR3_u.
 */
class deltaRR3_u: public ThObservable {
public:

    /**
     * @brief deltaRR3_u constructor.
     */
    deltaRR3_u(const StandardModel& SM_i);

    /**
     * @return @f$\delta RR3_u@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRR1_d
 * @ingroup SUSY
 * @brief deltaRR1_d.
 */
class deltaRR1_d: public ThObservable {
public:

    /**
     * @brief deltaRR1_d constructor.
     */
    deltaRR1_d(const StandardModel& SM_i);

    /**
     * @return @f$\delta RR1_d@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRR2_d
 * @ingroup SUSY
 * @brief deltaRR2_d.
 */
class deltaRR2_d: public ThObservable {
public:

    /**
     * @brief deltaRR2_d constructor.
     */
    deltaRR2_d(const StandardModel& SM_i);

    /**
     * @return @f$\delta RR2_d@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRR3_d
 * @ingroup SUSY
 * @brief deltaRR3_d.
 */
class deltaRR3_d: public ThObservable {
public:

    /**
     * @brief deltaRR3_d constructor.
     */
    deltaRR3_d(const StandardModel& SM_i);

    /**
     * @return @f$\delta RR3_d@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaLL1_l
 * @ingroup SUSY
 * @brief deltaLL1_l.
 */
class deltaLL1_l: public ThObservable {
public:

    /**
     * @brief deltaLL1_l constructor.
     */
    deltaLL1_l(const StandardModel& SM_i);

    /**
     * @return @f$\delta LL1_l@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaLL2_l
 * @ingroup SUSY
 * @brief deltaLL2_l.
 */
class deltaLL2_l: public ThObservable {
public:

    /**
     * @brief deltaLL2_l constructor.
     */
    deltaLL2_l(const StandardModel& SM_i);

    /**
     * @return @f$\delta LL2_l@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaLL3_l
 * @ingroup SUSY
 * @brief deltaLL3_l.
 */
class deltaLL3_l: public ThObservable {
public:

    /**
     * @brief deltaLL3_l constructor.
     */
    deltaLL3_l(const StandardModel& SM_i);

    /**
     * @return @f$\delta LL3_l@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRR1_e
 * @ingroup SUSY
 * @brief deltaRR1_e.
 */
class deltaRR1_e: public ThObservable {
public:

    /**
     * @brief deltaRR1_e constructor.
     */
    deltaRR1_e(const StandardModel& SM_i);

    /**
     * @return @f$\delta RR1_e@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRR2_e
 * @ingroup SUSY
 * @brief deltaRR2_e.
 */
class deltaRR2_e: public ThObservable {
public:

    /**
     * @brief deltaRR2_e constructor.
     */
    deltaRR2_e(const StandardModel& SM_i);

    /**
     * @return @f$\delta RR2_e@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRR3_e
 * @ingroup SUSY
 * @brief deltaRR3_e.
 */
class deltaRR3_e: public ThObservable {
public:

    /**
     * @brief deltaRR3_e constructor.
     */
    deltaRR3_e(const StandardModel& SM_i);

    /**
     * @return @f$\delta RR3_e@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_12_u
 * @ingroup SUSY
 * @brief deltaRL_12_u.
 */
class deltaRL_12_u: public ThObservable {
public:

    /**
     * @brief deltaRL_12_u constructor.
     */
    deltaRL_12_u(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL12u@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_13_u
 * @ingroup SUSY
 * @brief deltaRL_13_u.
 */
class deltaRL_13_u: public ThObservable {
public:

    /**
     * @brief deltaRL_13_u constructor.
     */
    deltaRL_13_u(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL13u@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_23_u
 * @ingroup SUSY
 * @brief deltaRL_23_u.
 */
class deltaRL_23_u: public ThObservable {
public:

    /**
     * @brief deltaRL_23_u constructor.
     */
    deltaRL_23_u(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL23u@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_12_d
 * @ingroup SUSY
 * @brief deltaRL_12_d.
 */
class deltaRL_12_d: public ThObservable {
public:

    /**
     * @brief deltaRL_12_d constructor.
     */
    deltaRL_12_d(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL12d@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_13_d
 * @ingroup SUSY
 * @brief deltaRL_13_d.
 */
class deltaRL_13_d: public ThObservable {
public:

    /**
     * @brief deltaRL_13_d constructor.
     */
    deltaRL_13_d(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL13d@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_23_d
 * @ingroup SUSY
 * @brief deltaRL_23_d.
 */
class deltaRL_23_d: public ThObservable {
public:

    /**
     * @brief deltaRL_23_d constructor.
     */
    deltaRL_23_d(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL23d@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_12_e
 * @ingroup SUSY
 * @brief deltaRL_12_e.
 */
class deltaRL_12_e: public ThObservable {
public:

    /**
     * @brief deltaRL_12_e constructor.
     */
    deltaRL_12_e(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL12e@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_13_e
 * @ingroup SUSY
 * @brief deltaRL_13_e.
 */
class deltaRL_13_e: public ThObservable {
public:

    /**
     * @brief deltaRL_13_e constructor.
     */
    deltaRL_13_e(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL13e@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_23_e
 * @ingroup SUSY
 * @brief deltaRL_23_e.
 */
class deltaRL_23_e: public ThObservable {
public:

    /**
     * @brief deltaRL_23_e constructor.
     */
    deltaRL_23_e(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL23e@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_21_e
 * @ingroup SUSY
 * @brief deltaRL_21_e.
 */
class deltaRL_21_e: public ThObservable {
public:

    /**
     * @brief deltaRL_21_e constructor.
     */
    deltaRL_21_e(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL21e@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_31_e
 * @ingroup SUSY
 * @brief deltaRL_31_e.
 */
class deltaRL_31_e: public ThObservable {
public:

    /**
     * @brief deltaRL_31_e constructor.
     */
    deltaRL_31_e(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL31e@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

/**
 * @class deltaRL_32_e
 * @ingroup SUSY
 * @brief deltaRL_32_e.
 */
class deltaRL_32_e: public ThObservable {
public:

    /**
     * @brief deltaRL_32_e constructor.
     */
    deltaRL_32_e(const StandardModel& SM_i);

    /**
     * @return @f$\delta_RL32e@f$
     */
    double computeThValue();

private:
    const SUSY& mySUSY;
};

#endif	/* SFERMIONDELTAS_H */
