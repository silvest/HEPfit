/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
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
 * @class m11_2 
 * @ingroup THDM 
 * @brief An observable class for the quadratic Higgs potential coupling @f$m_{11}^2@f$.
 * @details This class is used to compute the quadratic Higgs potential coupling @f$m_{11}^2@f$.
 */
class m11_2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    m11_2(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$m_{11}^2@f$.
     * @return @f$m_{11}^2@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class m22_2 
 * @ingroup THDM 
 * @brief An observable class for the quadratic Higgs potential coupling @f$m_{22}^2@f$.
 * @details This class is used to compute the quadratic Higgs potential coupling @f$m_{22}^2@f$.
 */
class m22_2 : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    m22_2(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$m_{22}^2@f$.
     * @return @f$m_{22}^2@f$
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

/**
 * @class g_hhh
 * @ingroup THDM 
 * @brief An observable class for the triple Higgs coupling @f$g_{hhh}@f$.
 * @details Taken from @cite Baglio:2014nea, equation (5) and multiplied by the SM coupling.
 *          Cross-checked with equation (F1) from @cite Gunion:2002zf
 */
class g_hhh : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    g_hhh(const StandardModel& SM_i);

    /**
     * @brief The triple Higgs coupling @f$g_{hhh}@f$.
     * @return @f$g_{hhh}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class g_hhHh
 * @ingroup THDM 
 * @brief An observable class for the triple Higgs coupling @f$g_{hhH}@f$.
 * @details Taken from @cite Baglio:2014nea, equation (5) and multiplied by the SM coupling.
 *          Cross-checked with equation (F1) from @cite Gunion:2002zf
 */
class g_hhHh : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    g_hhHh(const StandardModel& SM_i);

    /**
     * @brief The triple Higgs coupling @f$g_{hhH}@f$.
     * @return @f$g_{hhH}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class g_hHhHh
 * @ingroup THDM 
 * @brief An observable class for the triple Higgs coupling @f$g_{hHH}@f$.
 * @details Taken from @cite Baglio:2014nea, equation (5) and multiplied by the SM coupling.
 *          Cross-checked with equation (F1) from @cite Gunion:2002zf
 */
class g_hHhHh : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    g_hHhHh(const StandardModel& SM_i);

    /**
     * @brief The triple Higgs coupling @f$g_{hHH}@f$.
     * @return @f$g_{hHH}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class g_HhHhHh
 * @ingroup THDM 
 * @brief An observable class for the triple Higgs coupling @f$g_{HHH}@f$.
 * @details Taken from @cite Baglio:2014nea, equation (5) and multiplied by the SM coupling.
 *          Cross-checked with equation (F1) from @cite Gunion:2002zf
 */
class g_HhHhHh : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    g_HhHhHh(const StandardModel& SM_i);

    /**
     * @brief The triple Higgs coupling @f$g_{HHH}@f$.
     * @return @f$g_{HHH}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class g_hAA
 * @ingroup THDM 
 * @brief An observable class for the triple Higgs coupling @f$g_{hAA}@f$.
 * @details Taken from @cite Baglio:2014nea, equation (5) and multiplied by the SM coupling.
 *          Cross-checked with equation (F1) from @cite Gunion:2002zf
 */
class g_hAA : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    g_hAA(const StandardModel& SM_i);

    /**
     * @brief The triple Higgs coupling @f$g_{hAA}@f$.
     * @return @f$g_{hAA}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class g_HhAA
 * @ingroup THDM 
 * @brief An observable class for the triple Higgs coupling @f$g_{HAA}@f$.
 * @details Taken from @cite Baglio:2014nea, equation (5) and multiplied by the SM coupling.
 *          Cross-checked with equation (F1) from @cite Gunion:2002zf
 */
class g_HhAA : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    g_HhAA(const StandardModel& SM_i);

    /**
     * @brief The triple Higgs coupling @f$g_{HAA}@f$.
     * @return @f$g_{HAA}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class g_hHpHm
 * @ingroup THDM 
 * @brief An observable class for the triple Higgs coupling @f$g_{hH^+H^-}@f$.
 * @details Taken from @cite Baglio:2014nea, equation (5) and multiplied by the SM coupling.
 *          Cross-checked with equation (F1) from @cite Gunion:2002zf
 */
class g_hHpHm : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    g_hHpHm(const StandardModel& SM_i);

    /**
     * @brief The triple Higgs coupling @f$g_{hH^+H^-}@f$.
     * @return @f$g_{hH^+H^-}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class g_HhHpHm
 * @ingroup THDM 
 * @brief An observable class for the triple Higgs coupling @f$g_{HH^+H^-}@f$.
 * @details Taken from @cite Baglio:2014nea, equation (5) and multiplied by the SM coupling.
 *          Cross-checked with equation (F1) from @cite Gunion:2002zf
 */
class g_HhHpHm : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    g_HhHpHm(const StandardModel& SM_i);

    /**
     * @brief The triple Higgs coupling @f$g_{HH^+H^-}@f$.
     * @return @f$g_{HH^+H^-}@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Y1_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Y_1@f$.
 * @details This class is used to compute the quadratic Higgs potential coupling @f$Y_1@f$ in the Higgs basis.
 */
class Y1_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Y1_THDM(const StandardModel& SM_i);

    /**
     * @brief The quadratic coupling @f$Y_1@f$.
     * @return @f$Y_1@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Y2_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Y_2@f$.
 * @details This class is used to compute the quadratic Higgs potential coupling @f$Y_2@f$ in the Higgs basis.
 */
class Y2_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Y2_THDM(const StandardModel& SM_i);

    /**
     * @brief The quadratic coupling @f$Y_2@f$.
     * @return @f$Y_2@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Y3_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Y_3@f$.
 * @details This class is used to compute the quadratic Higgs potential coupling @f$Y_3@f$ in the Higgs basis.
 */
class Y3_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Y3_THDM(const StandardModel& SM_i);

    /**
     * @brief The quadratic coupling @f$Y_3@f$.
     * @return @f$Y_3@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Z1_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Z_1@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$Z_1@f$ in the Higgs basis.
 */
class Z1_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Z1_THDM(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$Z_1@f$.
     * @return @f$Z_1@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Z2_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Z_2@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$Z_2@f$ in the Higgs basis.
 */
class Z2_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Z2_THDM(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$Z_2@f$.
     * @return @f$Z_2@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Z3_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Z_3@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$Z_3@f$ in the Higgs basis.
 */
class Z3_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Z3_THDM(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$Z_3@f$.
     * @return @f$Z_3@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Z4_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Z_4@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$Z_4@f$ in the Higgs basis.
 */
class Z4_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Z4_THDM(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$Z_4@f$.
     * @return @f$Z_4@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Z5_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Z_5@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$Z_5@f$ in the Higgs basis.
 */
class Z5_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Z5_THDM(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$Z_5@f$.
     * @return @f$Z_5@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Z6_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Z_6@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$Z_6@f$ in the Higgs basis.
 */
class Z6_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Z6_THDM(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$Z_6@f$.
     * @return @f$Z_6@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

/**
 * @class Z7_THDM
 * @ingroup THDM 
 * @brief An observable class for the Higgs basis parameter @f$Z_7@f$.
 * @details This class is used to compute the quartic Higgs potential coupling @f$Z_7@f$ in the Higgs basis.
 */
class Z7_THDM : public ThObservable {
public:

    /**
     * @brief Constructor.
     */
    Z7_THDM(const StandardModel& SM_i);

    /**
     * @brief The quartic coupling @f$Z_7@f$.
     * @return @f$Z_7@f$
     */
    double computeThValue();

    const THDM * myTHDM;
};

#endif	/* THDMQUANTITIES_H */
