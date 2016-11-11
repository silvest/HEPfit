/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMQUANTITIES_H
#define	GENERALTHDMQUANTITIES_H

//#include <stdexcept>
#include "ThObservable.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

/**
 * @class mH1_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The lightest neutral Higgs mass
 */
class mH1_GTHDM: public ThObservable {
public:

    /**
     * @brief mH1_GTHDM constructor.
     */
    mH1_GTHDM(const StandardModel& SM_i);

    /**
     * @return The lightest of the physical neutral Higgs masses
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The medium-weighted neutral Higgs mass
 */
class mH2_GTHDM: public ThObservable {
public:

    /**
     * @brief mH2_GTHDM constructor.
     */
    mH2_GTHDM(const StandardModel& SM_i);

    /**
     * @return The second of the physical neutral Higgs masses
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH3_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The heaviest neutral Higgs mass
 */
class mH3_GTHDM: public ThObservable {
public:

    /**
     * @brief mH3_GTHDM constructor.
     */
    mH3_GTHDM(const StandardModel& SM_i);

    /**
     * @return The heaviest of the physical neutral Higgs masses
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH1sq_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The lightest neutral Higgs mass squared
 */
class mH1sq_GTHDM: public ThObservable {
public:

    /**
     * @brief mH1sq_GTHDM constructor.
     */
    mH1sq_GTHDM(const StandardModel& SM_i);

    /**
     * @return The lightest of the physical neutral Higgs masses squared
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH2sq_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The medium-weighted neutral Higgs mass squared
 */
class mH2sq_GTHDM: public ThObservable {
public:

    /**
     * @brief mH2sq_GTHDM constructor.
     */
    mH2sq_GTHDM(const StandardModel& SM_i);

    /**
     * @return The second of the physical neutral Higgs masses squared
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH3sq_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The heaviest neutral Higgs mass squared
 */
class mH3sq_GTHDM: public ThObservable {
public:

    /**
     * @brief mH3sq_GTHDM constructor.
     */
    mH3sq_GTHDM(const StandardModel& SM_i);

    /**
     * @return The heaviest of the physical neutral Higgs masses squared
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Msq11_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The (1,1) entry of the non-diagonal neutral Higgs mass square matrix
 */
class Msq11_GTHDM: public ThObservable {
public:

    /**
     * @brief Msq11_GTHDM constructor.
     */
    Msq11_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$M_{11}^2@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Msq12_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The (1,2) entry of the non-diagonal neutral Higgs mass square matrix
 */
class Msq12_GTHDM: public ThObservable {
public:

    /**
     * @brief Msq12_GTHDM constructor.
     */
    Msq12_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$M_{12}^2@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Msq13_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The (1,3) entry of the non-diagonal neutral Higgs mass square matrix
 */
class Msq13_GTHDM: public ThObservable {
public:

    /**
     * @brief Msq13_GTHDM constructor.
     */
    Msq13_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$M_{13}^2@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Msq22_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The (2,2) entry of the non-diagonal neutral Higgs mass square matrix
 */
class Msq22_GTHDM: public ThObservable {
public:

    /**
     * @brief Msq22_GTHDM constructor.
     */
    Msq22_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$M_{22}^2@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Msq23_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The (2,3) entry of the non-diagonal neutral Higgs mass square matrix
 */
class Msq23_GTHDM: public ThObservable {
public:

    /**
     * @brief Msq23_GTHDM constructor.
     */
    Msq23_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$M_{23}^2@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Msq33_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The (3,3) entry of the non-diagonal neutral Higgs mass square matrix
 */
class Msq33_GTHDM: public ThObservable {
public:

    /**
     * @brief Msq33_GTHDM constructor.
     */
    Msq33_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$M_{33}^2@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class M2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$M^2@f$
 */
class M2_GTHDM: public ThObservable {
public:

    /**
     * @brief M2_GTHDM constructor.
     */
    M2_GTHDM(const StandardModel& SM_i);

    /**
     * @return M2_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class m11_2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$m_{11}^2@f$
 */
class m11_2_GTHDM: public ThObservable {
public:

    /**
     * @brief m11_2_GTHDM constructor.
     */
    m11_2_GTHDM(const StandardModel& SM_i);

    /**
     * @return m11_2_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class m22_2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$m_{22}^2@f$
 */
class m22_2_GTHDM: public ThObservable {
public:

    /**
     * @brief m22_2_GTHDM constructor.
     */
    m22_2_GTHDM(const StandardModel& SM_i);

    /**
     * @return m22_2_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Imm12_2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$Im(m_{12}^2)@f$
 */
class Imm12_2_GTHDM: public ThObservable {
public:

    /**
     * @brief Imm12_2_GTHDM constructor.
     */
    Imm12_2_GTHDM(const StandardModel& SM_i);

    /**
     * @return Imm12_2_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda1_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\lambda_{1}@f$
 */
class lambda1_GTHDM: public ThObservable {
public:

    /**
     * @brief lambda1_GTHDM constructor.
     */
    lambda1_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambda1_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\lambda_{2}@f$
 */
class lambda2_GTHDM: public ThObservable {
public:

    /**
     * @brief lambda2_GTHDM constructor.
     */
    lambda2_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambda2_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda3_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\lambda_{3}@f$
 */
class lambda3_GTHDM: public ThObservable {
public:

    /**
     * @brief lambda3_GTHDM constructor.
     */
    lambda3_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambda3_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda4_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\lambda_{4}@f$
 */
class lambda4_GTHDM: public ThObservable {
public:

    /**
     * @brief lambda4_GTHDM constructor.
     */
    lambda4_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambda4_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Relambda5_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$Re(\lambda_{5})@f$
 */
class Relambda5_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambda5_GTHDM constructor.
     */
    Relambda5_GTHDM(const StandardModel& SM_i);

    /**
     * @return Relambda5_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class v1_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The vacuum expectation value of @f$\Phi_1@f$ in the generic basis
 */
class v1_GTHDM: public ThObservable {
public:

    /**
     * @brief v1_GTHDM constructor.
     */
    v1_GTHDM(const StandardModel& SM_i);

    /**
     * @return v1
     */
    double computeThValue();

    const GeneralTHDM * myGTHDM;
};

/**
 * @class v2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The vacuum expectation value of @f$\Phi_2@f$ in the generic basis
 */
class v2_GTHDM: public ThObservable {
public:

    /**
     * @brief v2_GTHDM constructor.
     */
    v2_GTHDM(const StandardModel& SM_i);

    /**
     * @return v2
     */
    double computeThValue();

    const GeneralTHDM * myGTHDM;
};

//    const GeneralTHDM * myGeneralTHDM;
//};
//
///**
// * @class 
// * @ingroup GeneralTHDM 
// * @brief parameter of the Higgs potential @f$lambda_{1}@f$
// */
//class GTHDM_lambda1: public ThObservable {
//public:
//
//    /**
//     * @brief lambda1 constructor.
//     */
//    GTHDM_lambda1(const StandardModel& SM_i);
//
//    /**
//     * @return lambda1
//     */
//    double computeThValue();
//
//    const GeneralTHDM * myGeneralTHDM;
//};
//
///**
// * @class 
// * @ingroup GeneralTHDM 
// * @brief parameter of the Higgs potential @f$lambda_{2}@f$
// */
//class GTHDM_lambda2: public ThObservable {
//public:
//
//    /**
//     * @brief lambda2 constructor.
//     */
//    GTHDM_lambda2(const StandardModel& SM_i);
//
//    /**
//     * @return lambda2
//     */
//    double computeThValue();
//
//    const GeneralTHDM * myGeneralTHDM;
//};
//
///**
// * @class 
// * @ingroup GeneralTHDM 
// * @brief parameter of the Higgs potential @f$lambda_{3}@f$
// */
//class GTHDM_lambda3: public ThObservable {
//public:
//
//    /**
//     * @brief lambda3 constructor.
//     */
//    GTHDM_lambda3(const StandardModel& SM_i);
//
//    /**
//     * @return lambda3
//     */
//    double computeThValue();
//
//    const GeneralTHDM * myGeneralTHDM;
//};
//
///**
// * @class 
// * @ingroup GeneralTHDM 
// * @brief parameter of the Higgs potential @f$lambda_{4}@f$
// */
//class GTHDM_lambda4: public ThObservable {
//public:
//
//    /**
//     * @brief lambda4 constructor.
//     */
//    GTHDM_lambda4(const StandardModel& SM_i);
//
//    /**
//     * @return lambda4
//     */
//    double computeThValue();
//
//    const GeneralTHDM * myGeneralTHDM;
//};
//
///**
// * @class 
// * @ingroup GeneralTHDM 
// * @brief parameter @f$M^2@f$ of @cite Kanemura:2015ska
// */
//class GTHDM_M_2: public ThObservable {
//public:
//
//    /**
//     * @brief M_2 constructor.
//     */
//    GTHDM_M_2(const StandardModel& SM_i);
//
//    /**
//     * @return M_2
//     */
//    double computeThValue();
//
//    const GeneralTHDM * myGeneralTHDM;
//};
//
///**
// * @class 
// * @ingroup GeneralTHDM 
// * @brief mass squared of neutral Higgs
// */
//class GTHDM_mH2_2: public ThObservable {
//public:
//
//    /**
//     * @brief mH2_2 constructor.
//     */
//    GTHDM_mH2_2(const StandardModel& SM_i);
//
//    /**
//     * @return mH2_2
//     */
//    double computeThValue();
//
//    const GeneralTHDM * myGeneralTHDM;
//};
//
///**
// * @class 
// * @ingroup GeneralTHDM 
// * @brief mass squared of neutral Higgs
// */
//class GTHDM_mH3_2: public ThObservable {
//public:
//
//    /**
//     * @brief mH3_2 constructor.
//     */
//    GTHDM_mH3_2(const StandardModel& SM_i);
//
//    /**
//     * @return mH3_2
//     */
//    double computeThValue();
//
//    const GeneralTHDM * myGeneralTHDM;
//};

#endif	/* GENERALTHDMQUANTITIES_H */
