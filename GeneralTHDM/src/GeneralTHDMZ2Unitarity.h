/* 
 * Copyright (C) 2025 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMZ2UNITARITY_H
#define GENERALTHDMZ2UNITARITY_H

#include "ThObservable.h"
#include "GeneralTHDMZ2Runner.h"
#include <gslpp.h>

/**
 * @class GeneralTHDMZ2Unitarity
 * @ingroup GeneralTHDM
 * @brief Base class for NLO perturbative unitarity conditions to the THDM with a Z2 symmetry.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the NLO unitarity constraints in terms
 * of the quartic parameters from the scalar potential, as defined in Cacchio:2016qyh
 */

/**
 * @class unitarity_Z2
 * @ingroup GeneralTHDM
 * @brief An observable class for NLO perturbative unitarity conditions to the THDM with a Z2 symmetry.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the necessary running variables to be used in the th. observables.
 */
class unitarity_Z2 : public ThObservable {
public:
    /**
     * @brief unitarity_Z2 constructor.
     */
   unitarity_Z2(const StandardModel& SM_i);

protected:
    GeneralTHDMZ2Runner myGTHDM;

    gslpp::matrix<double> myZ2_at_Q;

    void computeZ2_at_Q();
};

/************************************/
/* Eigenvalues of the even 00 block */
/************************************/

/**
 * @class unitarity00eveP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the even 00 block.
 */
class unitarity00eveP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity00eveP_Z2 constructor.
     */
    unitarity00eveP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{00+}@f$
     */
    double computeThValue();
};

/**
 * @class unitarity00eveM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the even 00 block.
 */
class unitarity00eveM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity00eveM_Z2 constructor.
     */
    unitarity00eveM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{00-}@f$
     */
    double computeThValue();
};


/***********************************/
/* Eigenvalues of the odd 00 block */
/***********************************/

/**
 * @class unitarity00oddP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the odd 00 block.
 */
class unitarity00oddP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity00oddP_Z2 constructor.
     */
    unitarity00oddP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{00+}@f$
     */
    double computeThValue();
};

/**
 * @class unitarity00oddM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the odd 00 block.
 */
class unitarity00oddM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity00oddM_Z2 constructor.
     */
    unitarity00oddM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{00-}@f$
     */
    double computeThValue();
};


/************************************/
/* Eigenvalues of the even 01 block */
/************************************/

/**
 * @class unitarity01eveP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the even 01 block.
 */
class unitarity01eveP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity01eveP_Z2 constructor.
     */
    unitarity01eveP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{01+}@f$
     */
    double computeThValue();
};

/**
 * @class unitarity01eveM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the even 01 block.
 */
class unitarity01eveM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity01eveM_Z2 constructor.
     */
    unitarity01eveM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{01-}@f$
     */
    double computeThValue();
};


/***********************************/
/* Eigenvalues of the odd 01 block */
/***********************************/

/**
 * @class unitarity01oddP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the odd 01 block.
 */
class unitarity01oddP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity01oddP_Z2 constructor.
     */
    unitarity01oddP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{01+}@f$
     */
    double computeThValue();
};

/**
 * @class unitarity01oddM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the odd 01 block.
 */
class unitarity01oddM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity01oddM_Z2 constructor.
     */
    unitarity01oddM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{01-}@f$
     */
    double computeThValue();
};


/**********************************/
/* Eigenvalue of the odd 10 block */
/**********************************/

/**
 * @class unitarity10odd_Z2
 * @ingroup GeneralTHDM
 * @brief Eigenvalue of the odd 10 block.
 */
class unitarity10odd_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity10odd_Z2 constructor.
     */
    unitarity10odd_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{10}@f$
     */
    double computeThValue();
};


/************************************/
/* Eigenvalues of the even 11 block */
/************************************/

/**
 * @class unitarity11eveP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the even 11 block.
 */
class unitarity11eveP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity11eveP_Z2 constructor.
     */
    unitarity11eveP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{11+}@f$
     */
    double computeThValue();
};

/**
 * @class unitarity11eveM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the even 11 block.
 */
class unitarity11eveM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity11eveM_Z2 constructor.
     */
    unitarity11eveM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{11-}@f$
     */
    double computeThValue();
};


/**********************************/
/* Eigenvalue of the odd 11 block */
/**********************************/

/**
 * @class unitarity11odd_Z2
 * @ingroup GeneralTHDM
 * @brief Eigenvalue of the odd 11 block.
 */
class unitarity11odd_Z2: public unitarity_Z2 {
public:

    /**
     * @brief unitarity11odd_Z2 constructor.
     */
    unitarity11odd_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{11}@f$
     */
    double computeThValue();
};


/***********************************/
/* R1p ratios of the even 00 block */
/***********************************/

/**
 * @class R1p00eveP_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the plus-sign eigenvalue of the even 00 block.
 */

//class R1p00eveP_Z2: public unitarity_Z2 {
//public:
//
//    /**
//     * @brief unitarity00eveP_Z2 constructor.
//     */
//    R1p00eveP_Z2(const StandardModel& SM_i);
//
//    /**
//     * @return @f$a^{\text{even}}_{00+}@f$
//     */
//    double computeThValue();
//    
//    private:
//    const GeneralTHDMZ2 * myGTHDMZ2;
//};

class R1p00eveP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p00eveP_Z2 constructor.
     */
    R1p00eveP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{even,\, NLO}}_{00+}}{a^{\text{even,\, LO}}_{00+}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};

/**
 * @class R1p00eveM_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the minus-sign eigenvalue of the even 00 block.
 */
class R1p00eveM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p00eveM_Z2 constructor.
     */
    R1p00eveM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{even,\, NLO}}_{00-}}{a^{\text{even,\, LO}}_{00-}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};


/**********************************/
/* R1p ratios of the odd 00 block */
/**********************************/

/**
 * @class R1p00oddP_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the plus-sign eigenvalue of the odd 00 block.
 */
class R1p00oddP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p00oddP_Z2 constructor.
     */
    R1p00oddP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{odd,\, NLO}}_{00+}}{a^{\text{odd,\, LO}}_{00+}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};

/**
 * @class R1p00oddM_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the minus-sign eigenvalue of the odd 00 block.
 */
class R1p00oddM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p00oddM_Z2 constructor.
     */
    R1p00oddM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{odd,\, NLO}}_{00-}}{a^{\text{odd,\, LO}}_{00-}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};


/***********************************/
/* R1p ratios of the even 01 block */
/***********************************/

/**
 * @class R1p01eveP_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the plus-sign eigenvalue of the even 01 block.
 */
class R1p01eveP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p01eveP_Z2 constructor.
     */
    R1p01eveP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{even,\, NLO}}_{01+}}{a^{\text{even,\, LO}}_{01+}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};

/**
 * @class R1p01eveM_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the minus-sign eigenvalue of the even 01 block.
 */
class R1p01eveM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p01eveM_Z2 constructor.
     */
    R1p01eveM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{even,\, NLO}}_{01-}}{a^{\text{even,\, LO}}_{01-}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};


/**********************************/
/* R1p ratios of the odd 01 block */
/**********************************/

/**
 * @class R1p01oddP_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the plus-sign eigenvalue of the odd 01 block.
 */
class R1p01oddP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p01oddP_Z2 constructor.
     */
    R1p01oddP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{odd,\, NLO}}_{01+}}{a^{\text{odd,\, LO}}_{01+}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};

/**
 * @class R1p01oddM_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the minus-sign eigenvalue of the odd 01 block.
 */
class R1p01oddM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p01oddM_Z2 constructor.
     */
    R1p01oddM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{odd,\, NLO}}_{01-}}{a^{\text{odd,\, LO}}_{01-}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};


/*********************************/
/* R1p ratio of the odd 10 block */
/*********************************/

/**
 * @class R1p10odd_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the eigenvalue of the odd 10 block.
 */
class R1p10odd_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p10odd_Z2 constructor.
     */
    R1p10odd_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{odd,\, NLO}}_{10}}{a^{\text{odd,\, LO}}_{10}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};


/***********************************/
/* R1p ratios of the even 11 block */
/***********************************/

/**
 * @class R1p11eveP_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the plus-sign eigenvalue of the even 11 block.
 */
class R1p11eveP_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p11eveP_Z2 constructor.
     */
    R1p11eveP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{even,\, NLO}}_{11+}}{a^{\text{even,\, LO}}_{11+}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};

/**
 * @class R1p11eveM_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the minus-sign eigenvalue of the even 11 block.
 */
class R1p11eveM_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p11eveM_Z2 constructor.
     */
    R1p11eveM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{even,\, NLO}}_{11-}}{a^{\text{even,\, LO}}_{11-}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};


/*********************************/
/* R1p ratio of the odd 11 block */
/*********************************/

/**
 * @class R1p11odd_Z2
 * @ingroup GeneralTHDM
 * @brief Perturbative unitarity ratio for the eigenvalue of the odd 11 block.
 */
class R1p11odd_Z2: public unitarity_Z2 {
public:

    /**
     * @brief R1p11odd_Z2 constructor.
     */
    R1p11odd_Z2(const StandardModel& SM_i);

    /**
     * @return @f$R_1^{\prime} = \vert\frac{a^{\text{odd,\, NLO}}_{11}}{a^{\text{odd,\, LO}}_{11}}\vert@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};


///////////////////////////////////////////////////////////////////////////
// An observable for tanb

/**
 * @class tanbetaZ2
 * @ingroup GeneralTHDM
 * @brief The tangent of beta.
 */
class tanbetaZ2: public ThObservable {
public:

    /**
     * @brief tanbetaZ2 constructor.
     */
    tanbetaZ2(const StandardModel& SM_i);

    /**
     * @return @f$ \tan(\beta) @f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2 * myGTHDMZ2;
};

#endif /* GENERALTHDMZ2UNITARITY_H */
