/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GMQUANTITIES_H
#define	GMQUANTITIES_H

#include <stdexcept>
#include "ThObservable.h"
#include "GeorgiMachacek.h"
#include "GMcache.h"

/**
 * @class tanbetaGM
 * @ingroup GeorgiMachacek 
 * @brief The tangent of beta.
 */
class tanbetaGM: public ThObservable {
public:

    /**
     * @brief tanbetaGM constructor.
     */
    tanbetaGM(const StandardModel& SM_i);

    /**
     * @return @f$\tan \beta@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class m1sqGM
 * @ingroup GeorgiMachacek 
 * @brief The dimension 2 coupling m1sq from the Higgs potential.
 */
class m1sqGM: public ThObservable {
public:

    /**
     * @brief m1sqGM constructor.
     */
    m1sqGM(const StandardModel& SM_i);

    /**
     * @return @f$m_1^2@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class m2sqGM
 * @ingroup GeorgiMachacek 
 * @brief The dimension 2 coupling m2sq from the Higgs potential.
 */
class m2sqGM: public ThObservable {
public:

    /**
     * @brief m2sqGM constructor.
     */
    m2sqGM(const StandardModel& SM_i);

    /**
     * @return @f$m_2^2@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class lambda1GM
 * @ingroup GeorgiMachacek 
 * @brief The quartic coupling lambda1 from the Higgs potential.
 */
class lambda1GM: public ThObservable {
public:

    /**
     * @brief lambda1GM constructor.
     */
    lambda1GM(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_1@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class lambda2GM
 * @ingroup GeorgiMachacek 
 * @brief The quartic coupling lambda2 from the Higgs potential.
 */
class lambda2GM: public ThObservable {
public:

    /**
     * @brief lambda2GM constructor.
     */
    lambda2GM(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_2@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class lambda3GM
 * @ingroup GeorgiMachacek 
 * @brief The quartic coupling lambda3 from the Higgs potential.
 */
class lambda3GM: public ThObservable {
public:

    /**
     * @brief lambda3GM constructor.
     */
    lambda3GM(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_3@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class lambda4GM
 * @ingroup GeorgiMachacek 
 * @brief The quartic coupling lambda4 from the Higgs potential.
 */
class lambda4GM: public ThObservable {
public:

    /**
     * @brief lambda4GM constructor.
     */
    lambda4GM(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_4@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class lambda5GM
 * @ingroup GeorgiMachacek 
 * @brief The quartic coupling lambda5 from the Higgs potential.
 */
class lambda5GM: public ThObservable {
public:

    /**
     * @brief lambda5GM constructor.
     */
    lambda5GM(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_5@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class vPhiGM
 * @ingroup GeorgiMachacek 
 * @brief The vacuum expectation values of the bi-doublet.
 */
class vPhiGM: public ThObservable {
public:

    /**
     * @brief vPhiGM constructor.
     */
    vPhiGM(const StandardModel& SM_i);

    /**
     * @return @f$v_\Phi@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};
//
///**
// * @class vDeltaGM
// * @ingroup GeorgiMachacek 
// * @brief The vacuum expectation values of the bi-triplet.
// */
//class vDeltaGM: public ThObservable {
//public:
//
//    /**
//     * @brief vDeltaGM constructor.
//     */
//    vDeltaGM(const StandardModel& SM_i);
//
//    /**
//     * @return @f$v_\Delta@f$
//     */
//    double computeThValue();
//
//    const GeorgiMachacek& myGM;
//};

/**
 * @class rh_gaga_GM
 * @ingroup GeorgiMachacek 
 * @brief The ratio of the GM partial Higgs decay width to two photons and the corresponding SM decay width.
 */
class rh_gaga_GM: public ThObservable {
public:

    /**
     * @brief rh_gaga_GM constructor.
     */
    rh_gaga_GM(const StandardModel& SM_i);

    /**
     * @return @f$r^{(h)}_{\gamma\gamma}@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class rh_Zga_GM
 * @ingroup GeorgiMachacek 
 * @brief The ratio of the GM partial Higgs decay width into a $Z$ boson and a photon and the corresponding SM decay width.
 */
class rh_Zga_GM: public ThObservable {
public:

    /**
     * @brief rh_Zga_GM constructor.
     */
    rh_Zga_GM(const StandardModel& SM_i);

    /**
     * @return @f$r^{(h)}_{Z\gamma}@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class GMmass_mHh
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmass_mHh: public ThObservable {
public:

    /**
     * @brief GMmass_mHh constructor.
     */
    GMmass_mHh(const StandardModel& SM_i);

    /**
     * @return @f$m_H@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmass_mA
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmass_mA: public ThObservable {
public:

    /**
     * @brief GMmass_mA constructor.
     */
    GMmass_mA(const StandardModel& SM_i);

    /**
     * @return @f$m_A@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmass_mH5
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmass_mH5: public ThObservable {
public:

    /**
     * @brief GMmass_mH5 constructor.
     */
    GMmass_mH5(const StandardModel& SM_i);

    /**
     * @return @f$m_{H_5}@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mHlmmHh
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mHlmmHh: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mHlmmHh constructor.
     */
    GMmassdifference_mHlmmHh(const StandardModel& SM_i);

    /**
     * @return @f$m_h-m_H@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mHhmmHl
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mHhmmHl: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mHhmmHl constructor.
     */
    GMmassdifference_mHhmmHl(const StandardModel& SM_i);

    /**
     * @return @f$m_H-m_h@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mHlmmA
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mHlmmA: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mHlmmA constructor.
     */
    GMmassdifference_mHlmmA(const StandardModel& SM_i);

    /**
     * @return @f$m_h-m_A@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mAmmHl
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mAmmHl: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mAmmHl constructor.
     */
    GMmassdifference_mAmmHl(const StandardModel& SM_i);

    /**
     * @return @f$m_A-m_h@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mHlmmH5
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mHlmmH5: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mHlmmH5 constructor.
     */
    GMmassdifference_mHlmmH5(const StandardModel& SM_i);

    /**
     * @return @f$m_h-m_{H_5}@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mH5mmHl
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mH5mmHl: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mH5mmHl constructor.
     */
    GMmassdifference_mH5mmHl(const StandardModel& SM_i);

    /**
     * @return @f$m_{H_5}-m_h@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mHhmmA
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mHhmmA: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mHhmmA constructor.
     */
    GMmassdifference_mHhmmA(const StandardModel& SM_i);

    /**
     * @return @f$m_H-m_A@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mAmmHh
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mAmmHh: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mAmmHh constructor.
     */
    GMmassdifference_mAmmHh(const StandardModel& SM_i);

    /**
     * @return @f$m_A-m_H@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mHhmmH5
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mHhmmH5: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mHhmmH5 constructor.
     */
    GMmassdifference_mHhmmH5(const StandardModel& SM_i);

    /**
     * @return @f$m_H-m_{H_5}@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mH5mmHh
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mH5mmHh: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mH5mmHh constructor.
     */
    GMmassdifference_mH5mmHh(const StandardModel& SM_i);

    /**
     * @return @f$m_{H_5}-m_H@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mAmmH5
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mAmmH5: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mAmmH5 constructor.
     */
    GMmassdifference_mAmmH5(const StandardModel& SM_i);

    /**
     * @return @f$m_A-m_{H_5}@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMmassdifference_mH5mmA
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMmassdifference_mH5mmA: public ThObservable {
public:

    /**
     * @brief GMmassdifference_mH5mmA constructor.
     */
    GMmassdifference_mH5mmA(const StandardModel& SM_i);

    /**
     * @return @f$m_{H_5}-m_A@f$
     */
    double computeThValue();

    const GeorgiMachacek * myGM;
};

/**
 * @class GMGammah
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMGammah: public ThObservable {
public:

    /**
     * @brief GMGammah constructor.
     */
    GMGammah(const StandardModel& SM_i);

    /**
     * @return @f$\Gamma_h@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class GMGammaH1
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMGammaH1: public ThObservable {
public:

    /**
     * @brief GMGammaH1 constructor.
     */
    GMGammaH1(const StandardModel& SM_i);

    /**
     * @return @f$\Gamma_1@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class GMGammaH3
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMGammaH3: public ThObservable {
public:

    /**
     * @brief GMGammaH3 constructor.
     */
    GMGammaH3(const StandardModel& SM_i);

    /**
     * @return @f$\Gamma_3@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class GMGammaH3p
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMGammaH3p: public ThObservable {
public:

    /**
     * @brief GMGammaH3p constructor.
     */
    GMGammaH3p(const StandardModel& SM_i);

    /**
     * @return @f$\Gamma_{3+}@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class GMGammaH5
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMGammaH5: public ThObservable {
public:

    /**
     * @brief GMGammaH5 constructor.
     */
    GMGammaH5(const StandardModel& SM_i);

    /**
     * @return @f$\Gamma_5@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class GMGammaH5p
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMGammaH5p: public ThObservable {
public:

    /**
     * @brief GMGammaH5p constructor.
     */
    GMGammaH5p(const StandardModel& SM_i);

    /**
     * @return @f$\Gamma_{5+}@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

/**
 * @class GMGammaH5pp
 * @ingroup GeorgiMachacek
 * @brief 
 */
class GMGammaH5pp: public ThObservable {
public:

    /**
     * @brief GMGammaH5pp constructor.
     */
    GMGammaH5pp(const StandardModel& SM_i);

    /**
     * @return @f$\Gamma_{5++}@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

#endif	/* GMQUANTITIES_H */
