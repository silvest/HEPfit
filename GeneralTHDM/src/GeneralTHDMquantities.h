/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMQUANTITIES_H
#define	GENERALTHDMQUANTITIES_H

#include "ThObservable.h"


class GeneralTHDM;
class GeneralTHDMcache;

/**
 * @class tanbeta_GTHDM
 * @ingroup GeneralTHDM 
 * @brief The tangent of beta.
 */
class tanbeta_GTHDM: public ThObservable {
public:

    /**
     * @brief tanbeta constructor.
     */
    tanbeta_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\tan \beta@f$
     */
    double computeThValue();

    const GeneralTHDM * myGTHDM;
};


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
//    const GeneralTHDM& myGTHDM;
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
//    const GeneralTHDM& myGTHDM;
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
//    const GeneralTHDM& myGTHDM;
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
//    const GeneralTHDM& myGTHDM;
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
//    const GeneralTHDM& myGTHDM;
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
//    const GeneralTHDM& myGTHDM;
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
//    const GeneralTHDM& myGTHDM;
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
 * @class lambda1H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\lambda_{1}@f$ in the Higgs basis
 */
class lambda1H_GTHDM: public ThObservable {
public:

    /**
     * @brief lambda1H_GTHDM constructor.
     */
    lambda1H_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambda1H_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda2H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\lambda_{2}@f$ in the Higgs basis
 */
class lambda2H_GTHDM: public ThObservable {
public:

    /**
     * @brief lambda2H_GTHDM constructor.
     */
    lambda2H_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambda2H_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda3H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\lambda_{3}@f$ in the Higgs basis
 */
class lambda3H_GTHDM: public ThObservable {
public:

    /**
     * @brief lambda3H_GTHDM constructor.
     */
    lambda3H_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambda3H_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda4H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\lambda_{4}@f$ in the Higgs basis
 */
class lambda4H_GTHDM: public ThObservable {
public:

    /**
     * @brief lambda4H_GTHDM constructor.
     */
    lambda4H_GTHDM(const StandardModel& SM_i);

    /**
     * @return lambda4H_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Relambda5H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\Relambda_{5}@f$ in the Higgs basis
 */
class Relambda5H_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambda5H_GTHDM constructor.
     */
    Relambda5H_GTHDM(const StandardModel& SM_i);

    /**
     * @return Relambda5H_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Imlambda5H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\Imlambda_{5}@f$ in the Higgs basis
 */
class Imlambda5H_GTHDM: public ThObservable {
public:

    /**
     * @brief Imlambda5H_GTHDM constructor.
     */
    Imlambda5H_GTHDM(const StandardModel& SM_i);

    /**
     * @return Imlambda5H_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class Relambda6H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\Relambda_{6}@f$ in the Higgs basis
 */
class Relambda6H_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambda6H_GTHDM constructor.
     */
    Relambda6H_GTHDM(const StandardModel& SM_i);

    /**
     * @return Relambda6H_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Imlambda6H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\Imlambda_{6}@f$ in the Higgs basis
 */
class Imlambda6H_GTHDM: public ThObservable {
public:

    /**
     * @brief Imlambda6H_GTHDM constructor.
     */
    Imlambda6H_GTHDM(const StandardModel& SM_i);

    /**
     * @return Imlambda6H_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Relambda7H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\Relambda_{7}@f$ in the Higgs basis
 */
class Relambda7H_GTHDM: public ThObservable {
public:

    /**
     * @brief Relambda7H_GTHDM constructor.
     */
    Relambda7H_GTHDM(const StandardModel& SM_i);

    /**
     * @return Relambda7H_GTHDM
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Imlambda7H_GTHDM
 * @ingroup GeneralTHDM 
 * @brief parameter of the Higgs potential @f$\Imlambda_{7}@f$ in the Higgs basis
 */
class Imlambda7H_GTHDM: public ThObservable {
public:

    /**
     * @brief Imlambda7H_GTHDM constructor.
     */
    Imlambda7H_GTHDM(const StandardModel& SM_i);

    /**
     * @return Imlambda7H_GTHDM
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


/**
 * @class cosalpha1_GTHDM
 * @ingroup GeneralTHDM 
 * @brief cosalpha1_GTHDM
 */
class cosalpha1_GTHDM: public ThObservable {
public:

    /**
     * @brief cosalpha1 constructor.
     */
    cosalpha1_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$cosalpha1_GTHDM@f$
     */
    double computeThValue();

   const GeneralTHDM * myGTHDM;
};

/**
 * @class m1_2
 * @ingroup GeneralTHDM 
 * @brief 
 */
class m1_2: public ThObservable {
public:

    /**
     * @brief 
     */
    m1_2(const StandardModel& SM_i);

    /**
     * @return m1_2
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class m2_2
 * @ingroup GeneralTHDM 
 * @brief 
 */
class m2_2: public ThObservable {
public:

    /**
     * @brief 
     */
    m2_2(const StandardModel& SM_i);

    /**
     * @return m2_2
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};


/**
 * @class m3_2
 * @ingroup GeneralTHDM 
 * @brief 
 */
class m3_2: public ThObservable {
public:

    /**
     * @brief 
     */
    m3_2(const StandardModel& SM_i);

    /**
     * @return m3_2
     */
    double computeThValue();
private:
    const GeneralTHDM& myGTHDM;
};



/**
 * @class Q_stGTHDM
 * @ingroup GeneralTHDM
 * @brief Q_stTHDMW.
 */
class Q_stGTHDM: public ThObservable {
public:

    /**
     * @brief Q_stGTHDM constructor.
     */
    Q_stGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$Q_{\text{stability}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class DeltaQ_GTHDM
 * @ingroup GeneralTHDM
 * @brief DeltaQ_GTHDM.
 */
class DeltaQ_GTHDM: public ThObservable {
public:

    /**
     * @brief DeltaQ_GTHDM constructor.
     */
    DeltaQ_GTHDM(const StandardModel& SM_i);

    /**
     * @return @f$Q_{\text{GTHDM}}-Q_{\text{stability}}@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class g1atQGTHDM
 * @ingroup GeneralTHDM
 * @brief g1atQGTHDM.
 */
class g1atQGTHDM: public ThObservable {
public:

    /**
     * @brief g1atQGTHDM constructor.
     */
    g1atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$g_1(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class g2atQGTHDM
 * @ingroup GeneralTHDM
 * @brief g2atQGTHDM.
 */
class g2atQGTHDM: public ThObservable {
public:

    /**
     * @brief g2atQGTHDM constructor.
     */
    g2atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$g_2(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class g3atQGTHDM
 * @ingroup GeneralTHDM
 * @brief g3atQGTHDM.
 */
class g3atQGTHDM: public ThObservable {
public:

    /**
     * @brief g3atQGTHDM constructor.
     */
    g3atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$g_3(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class etaU1atQGTHDM
 * @ingroup GeneralTHDM
 * @brief etaU1atQGTHDM.
 */
class etaU1atQGTHDM: public ThObservable {
public:

    /**
     * @brief etaU1atQGTHDM constructor.
     */
    etaU1atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$Y_{t1}(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class etaU2atQGTHDM
 * @ingroup GeneralTHDM
 * @brief etaU2atQGTHDM.
 */
class etaU2atQGTHDM: public ThObservable {
public:

    /**
     * @brief etaU2atQGTHDM constructor.
     */
    etaU2atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$Y_{t2}(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class etaD1atQGTHDM
 * @ingroup GeneralTHDM
 * @brief etaD1atQGTHDM.
 */
class etaD1atQGTHDM: public ThObservable {
public:

    /**
     * @brief etaD1atQGTHDM constructor.
     */
    etaD1atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$Y_{b1}(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class etaD2atQGTHDM
 * @ingroup GeneralTHDM
 * @brief etaD2atQGTHDM.
 */
class etaD2atQGTHDM: public ThObservable {
public:

    /**
     * @brief etaD2atQGTHDM constructor.
     */
    etaD2atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$Y_{b2}(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class etaL1atQGTHDM
 * @ingroup GeneralTHDM
 * @brief etaL1atQGTHDM.
 */
class etaL1atQGTHDM: public ThObservable {
public:

    /**
     * @brief etaL1atQGTHDM constructor.
     */
    etaL1atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$Y_{\tau 1}(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class etaL2atQGTHDM
 * @ingroup GeneralTHDM
 * @brief etaL2atQGTHDM.
 */
class etaL2atQGTHDM: public ThObservable {
public:

    /**
     * @brief etaL2atQGTHDM constructor.
     */
    etaL2atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$Y_{\tau 2}(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda1atQGTHDM
 * @ingroup GeneralTHDM
 * @brief lambda1atQGTHDM.
 */
class lambda1atQGTHDM: public ThObservable {
public:

    /**
     * @brief lambda1atQGTHDM constructor.
     */
    lambda1atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_1(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda2atQGTHDM
 * @ingroup GeneralTHDM
 * @brief lambda2atQGTHDM.
 */
class lambda2atQGTHDM: public ThObservable {
public:

    /**
     * @brief lambda2atQGTHDM constructor.
     */
    lambda2atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_2(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda3atQGTHDM
 * @ingroup GeneralTHDM
 * @brief lambda3atQGTHDM.
 */
class lambda3atQGTHDM: public ThObservable {
public:

    /**
     * @brief lambda3atQGTHDM constructor.
     */
    lambda3atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_3(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class lambda4atQGTHDM
 * @ingroup GeneralTHDM
 * @brief lambda4atQGTHDM.
 */
class lambda4atQGTHDM: public ThObservable {
public:

    /**
     * @brief lambda4atQGTHDM constructor.
     */
    lambda4atQGTHDM(const StandardModel& SM_i);

    /**
     * @return @f$\lambda_4(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Relambda5atQGTHDM
 * @ingroup GeneralTHDM
 * @brief Relambda5atQGTHDM.
 */
class Relambda5atQGTHDM: public ThObservable {
public:

    /**
     * @brief Relambda5atQGTHDM constructor.
     */
    Relambda5atQGTHDM(const StandardModel& SM_i);

    /**
     * @return Re @f$\lambda_5(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Relambda6atQGTHDM
 * @ingroup GeneralTHDM
 * @brief Relambda6atQGTHDM.
 */
class Relambda6atQGTHDM: public ThObservable {
public:

    /**
     * @brief Relambda6atQGTHDM constructor.
     */
    Relambda6atQGTHDM(const StandardModel& SM_i);

    /**
     * @return Re @f$\lambda_6(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

/**
 * @class Relambda7atQGTHDM
 * @ingroup GeneralTHDM
 * @brief Relambda7atQGTHDM.
 */
class Relambda7atQGTHDM: public ThObservable {
public:

    /**
     * @brief Relambda7atQGTHDM constructor.
     */
    Relambda7atQGTHDM(const StandardModel& SM_i);

    /**
     * @return Re @f$\lambda_7(Q)@f$
     */
    double computeThValue();

private:
    const GeneralTHDM& myGTHDM;
};

#endif	/* GENERALTHDMQUANTITIES_H */
