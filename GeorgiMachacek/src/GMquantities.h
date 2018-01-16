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

/**
 * @class vDeltaGM
 * @ingroup GeorgiMachacek 
 * @brief The vacuum expectation values of the bi-triplet.
 */
class vDeltaGM: public ThObservable {
public:

    /**
     * @brief vDeltaGM constructor.
     */
    vDeltaGM(const StandardModel& SM_i);

    /**
     * @return @f$v_\Delta@f$
     */
    double computeThValue();

    const GeorgiMachacek& myGM;
};

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

#endif	/* GMQUANTITIES_H */
