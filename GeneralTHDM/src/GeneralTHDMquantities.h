/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMQUANTITIES_H
#define	GENERALTHDMQUANTITIES_H

#include <stdexcept>
#include "ThObservable.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"


/**
 * @class mH1_2
 * @ingroup GeneralTHDM 
 * @brief The mass of the SM Higgs
 */
class mH1_2: public ThObservable {
public:

    /**
     * @brief mH1_GTHDM constructor.
     */
    mH1_2(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH2_2
 * @ingroup GeneralTHDM 
 * @brief The mass of the SM Higgs
 */
class mH2_2: public ThObservable {
public:

    /**
     * @brief mH2_GTHDM constructor.
     */
    mH2_2(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class mH3_2
 * @ingroup GeneralTHDM 
 * @brief The mass of the SM Higgs
 */
class mH3_2: public ThObservable {
public:

    /**
     * @brief mH3_GTHDM constructor.
     */
    mH3_2(const StandardModel& SM_i);

    /**
     * @return The value of mH3
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class mH1_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass of the SM Higgs
 */
class mH1_GTHDM: public ThObservable {
public:

    /**
     * @brief mH1_GTHDM constructor.
     */
    mH1_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass of the second neutral Higgs
 */
class mH2_GTHDM: public ThObservable {
public:

    /**
     * @brief mH2_GTHDM constructor.
     */
    mH2_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of mH2
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH3_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass of the third neutral Higgs
 */
class mH3_GTHDM: public ThObservable {
public:

    /**
     * @brief mH3_GTHDM constructor.
     */
    mH3_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of mH3
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mHlight_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass of the lightest neutral Higgs
 */
class mHlight_GTHDM: public ThObservable {
public:

    /**
     * @brief mHlight_GTHDM constructor.
     */
    mHlight_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of mHlight
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mHmedium_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass of the medium-weighted neutral Higgs
 */
class mHmedium_GTHDM: public ThObservable {
public:

    /**
     * @brief mHmedium_GTHDM constructor.
     */
    mHmedium_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of mHmedium
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mHheavy_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass of the heaviest neutral Higgs
 */
class mHheavy_GTHDM: public ThObservable {
public:

    /**
     * @brief mHheavy_GTHDM constructor.
     */
    mHheavy_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of mHheavy
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mHp_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass of the charged Higgs
 */
class mHp_GTHDM: public ThObservable {
public:

    /**
     * @brief mHp_GTHDM constructor.
     */
    mHp_GTHDM(const StandardModel& SM_i);

    /**
     * @return The mass of the charged Higgs
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class mH3mmH2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass difference between the third and the second neutral Higgs
 */
class mH3mmH2_GTHDM: public ThObservable {
public:

    /**
     * @brief mH3mmH2_GTHDM constructor.
     */
    mH3mmH2_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of the mass of H3 minus the mass of H2
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH3mmH1_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass difference between the third neutral Higgs and the SM Higgs
 */
class mH3mmH1_GTHDM: public ThObservable {
public:

    /**
     * @brief mH3mmH1_GTHDM constructor.
     */
    mH3mmH1_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of the mass of H3 minus the mass of H1
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH3mmHp_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass difference between the third neutral Higgs and the charged Higgs
 */
class mH3mmHp_GTHDM: public ThObservable {
public:

    /**
     * @brief mH3mmHp_GTHDM constructor.
     */
    mH3mmHp_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of the mass of H3 minus the mass of the charged Higgs
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH2mmHp_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass difference between the second neutral Higgs and the charged Higgs
 */
class mH2mmHp_GTHDM: public ThObservable {
public:

    /**
     * @brief mH2mmHp_GTHDM constructor.
     */
    mH2mmHp_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of the mass of H2 minus the mass of the charged Higgs
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH2mmH1_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass difference between the second neutral Higgs and the SM Higgs
 */
class mH2mmH1_GTHDM: public ThObservable {
public:

    /**
     * @brief mH2mmH1_GTHDM constructor.
     */
    mH2mmH1_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of the mass of H2 minus the mass of H1
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mHpmmH1_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The mass difference between the charged Higgs and the SM Higgs
 */
class mHpmmH1_GTHDM: public ThObservable {
public:

    /**
     * @brief mHpmmH1_GTHDM constructor.
     */
    mHpmmH1_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of the mass of the charged Higgs minus the mass of H1
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class mH1sq_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The SM Higgs mass squared
 */
class mH1sq_GTHDM: public ThObservable {
public:

    /**
     * @brief mH1sq_GTHDM constructor.
     */
    mH1sq_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of mH1 squared
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH2sq_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The second neutral Higgs mass squared
 */
class mH2sq_GTHDM: public ThObservable {
public:

    /**
     * @brief mH2sq_GTHDM constructor.
     */
    mH2sq_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of mH2 squared
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class mH3sq_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The third neutral Higgs mass squared
 */
class mH3sq_GTHDM: public ThObservable {
public:

    /**
     * @brief mH3sq_GTHDM constructor.
     */
    mH3sq_GTHDM(const StandardModel& SM_i);

    /**
     * @return The value of mH3 squared
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
 * @class Rem12_2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$Re(m_{12}^2)@f$
 */
class Rem12_2_GTHDM: public ThObservable {
public:

    /**
     * @brief Rem12_2_GTHDM constructor.
     */
    Rem12_2_GTHDM(const StandardModel& SM_i);

    /**
     * @return Rem12_2_GTHDM
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
 * @class R11_GTHDM
 * @ingroup GeneralTHDM 
 * @brief Element (1,1) of the ortogonal matrix determining the mass e/states
 */
class R11_GTHDM: public ThObservable {
public:

    /**
     * @brief R11_GTHDM constructor.
     */
    R11_GTHDM(const StandardModel& SM_i);

    /**
     * @return R11_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class R12_GTHDM
 * @ingroup GeneralTHDM 
 * @brief Element (1,2) of the ortogonal matrix determining the mass e/states
 */
class R12_GTHDM: public ThObservable {
public:

    /**
     * @brief R12_GTHDM constructor.
     */
    R12_GTHDM(const StandardModel& SM_i);

    /**
     * @return R12_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class R13_GTHDM
 * @ingroup GeneralTHDM 
 * @brief Element (1,3) of the ortogonal matrix determining the mass e/states
 */
class R13_GTHDM: public ThObservable {
public:

    /**
     * @brief R13_GTHDM constructor.
     */
    R13_GTHDM(const StandardModel& SM_i);

    /**
     * @return R13_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class R11_GTHDM
 * @ingroup GeneralTHDM 
 * @brief Element (2,1) of the ortogonal matrix determining the mass e/states
 */
class R21_GTHDM: public ThObservable {
public:

    /**
     * @brief R21_GTHDM constructor.
     */
    R21_GTHDM(const StandardModel& SM_i);

    /**
     * @return R21_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class R22_GTHDM
 * @ingroup GeneralTHDM 
 * @brief Element (2,2) of the ortogonal matrix determining the mass e/states
 */
class R22_GTHDM: public ThObservable {
public:

    /**
     * @brief R22_GTHDM constructor.
     */
    R22_GTHDM(const StandardModel& SM_i);

    /**
     * @return R22_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class R23_GTHDM
 * @ingroup GeneralTHDM 
 * @brief Element (2,3) of the ortogonal matrix determining the mass e/states
 */
class R23_GTHDM: public ThObservable {
public:

    /**
     * @brief R23_GTHDM constructor.
     */
    R23_GTHDM(const StandardModel& SM_i);

    /**
     * @return R23_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class R31_GTHDM
 * @ingroup GeneralTHDM 
 * @brief Element (3,1) of the ortogonal matrix determining the mass e/states
 */
class R31_GTHDM: public ThObservable {
public:

    /**
     * @brief R31_GTHDM constructor.
     */
    R31_GTHDM(const StandardModel& SM_i);

    /**
     * @return R31_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class R32_GTHDM
 * @ingroup GeneralTHDM 
 * @brief Element (3,2) of the ortogonal matrix determining the mass e/states
 */
class R32_GTHDM: public ThObservable {
public:

    /**
     * @brief R32_GTHDM constructor.
     */
    R32_GTHDM(const StandardModel& SM_i);

    /**
     * @return R32_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class R33_GTHDM
 * @ingroup GeneralTHDM 
 * @brief Element (3,3) of the ortogonal matrix determining the mass e/states
 */
class R33_GTHDM: public ThObservable {
public:

    /**
     * @brief R13_GTHDM constructor.
     */
    R33_GTHDM(const StandardModel& SM_i);

    /**
     * @return R33_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class maa_2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @maa_2_GTHDM@f,$
 */
class maa_2_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambda5_GTHDM constructor.
     */
    maa_2_GTHDM(const StandardModel& SM_i);

    /**
     * @return maa_2_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class mbb_2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @mbb_2_GTHDM@f,$
 */
class mbb_2_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambda5_GTHDM constructor.
     */
    mbb_2_GTHDM(const StandardModel& SM_i);

    /**
     * @return mbb_2_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Remab_2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @Remab_2_GTHDM@f,$
 */
class Remab_2_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambda5_GTHDM constructor.
     */
    Remab_2_GTHDM(const StandardModel& SM_i);

    /**
     * @return Remab_2_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Immab_2_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @Immab_2_GTHDM@f,$
 */
class Immab_2_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambda5_GTHDM constructor.
     */
    Immab_2_GTHDM(const StandardModel& SM_i);

    /**
     * @return Immab_2_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambdaa_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @lambdaa_GTHDM@f,$
 */
class lambdaa_GTHDM: public ThObservable {
public:

    /**
     * @brief lambdaa_GTHDM constructor.
     */
    lambdaa_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambdaa_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class lambdab_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @lambdab_GTHDM@f,$
 */
class lambdab_GTHDM: public ThObservable {
public:

    /**
     * @brief lambdab_GTHDM constructor.
     */
    lambdab_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambdab_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class lambdac_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @lambdac_GTHDM@f,$
 */
class lambdac_GTHDM: public ThObservable {
public:

    /**
     * @brief lambdac_GTHDM constructor.
     */
    lambdac_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambdac_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class lambdad_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @lambdad_GTHDM@f,$
 */
class lambdad_GTHDM: public ThObservable {
public:

    /**
     * @brief lambdad_GTHDM constructor.
     */
    lambdad_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambdad_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Relambdae_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @Relambdae_GTHDM@f,$
 */
class Relambdae_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambdae_GTHDM constructor.
     */
    Relambdae_GTHDM(const StandardModel& SM_i);

    /**
     * @return Relambdae_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Relambdaf_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @Relambdaf_GTHDM@f,$
 */
class Relambdaf_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambdaf_GTHDM constructor.
     */
    Relambdaf_GTHDM(const StandardModel& SM_i);

    /**
     * @return Relambdaf_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class Relambdag_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @Relambdag_GTHDM@f,$
 */
class Relambdag_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambdag_GTHDM constructor.
     */
    Relambdag_GTHDM(const StandardModel& SM_i);

    /**
     * @return Relambdag_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Imlambdae_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @Imlambdae_GTHDM@f,$
 */
class Imlambdae_GTHDM: public ThObservable {
public:

    /**
     * @brief Imlambdae_GTHDM constructor.
     */
    Imlambdae_GTHDM(const StandardModel& SM_i);

    /**
     * @return Imlambdae_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Imlambdaf_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @Imlambdaf_GTHDM@f,$
 */
class Imlambdaf_GTHDM: public ThObservable {
public:

    /**
     * @brief Imlambdaf_GTHDM constructor.
     */
    Imlambdaf_GTHDM(const StandardModel& SM_i);

    /**
     * @return Imlambdaf_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Imlambdag_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential (in the Higgs basis) @Imlambdag_GTHDM@f,$
 */
class Imlambdag_GTHDM: public ThObservable {
public:

    /**
     * @brief Imlambdag_GTHDM constructor.
     */
    Imlambdag_GTHDM(const StandardModel& SM_i);

    /**
     * @return Imlambdag_GTHDM
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

/**
 * @class Resigmau
 * @ingroup GeneralTHDM 
 * @brief Real part of the up-type Yukawa coupling parameter sigma
 */
class Resigmau: public ThObservable {
public:

    /**
     * @brief Resigmau constructor.
     */
    Resigmau(const StandardModel& SM_i);

    /**
     * @return Resigmau
     */
    double computeThValue();

    const GeneralTHDM * myGTHDM;
};

#endif	/* GENERALTHDMQUANTITIES_H */
