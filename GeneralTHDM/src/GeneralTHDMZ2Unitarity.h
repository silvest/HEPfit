/* 
 * Copyright (C) 2025 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMZ2UNITARITY_H
#define GENERALTHDMZ2UNITARITY_H

#include "ThObservable.h"

class GeneralTHDMZ2;
class GeneralTHDMcache;

/**
 * @class GeneralTHDMZ2Unitarity
 * @ingroup GeneralTHDM
 * @brief Base class for NLO perturbative unitarity conditions to the THDM with a Z2 symmetry.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the NLO unitarity constraints in terms
 * of the quartic parameters from the scalar potential, as defined in Cacchio:2016qyh
 */

/************************************/
/* Eigenvalues of the even 00 block */
/************************************/

/**
 * @class unitarity00eveP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the even 00 block.
 */
class unitarity00eveP_Z2: public ThObservable {
public:

    /**
     * @brief unitarity00eveP_Z2 constructor.
     */
    unitarity00eveP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{00+}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};

/**
 * @class unitarity00eveM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the even 00 block.
 */
class unitarity00eveM_Z2: public ThObservable {
public:

    /**
     * @brief unitarity00eveM_Z2 constructor.
     */
    unitarity00eveM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{00-}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};


/***********************************/
/* Eigenvalues of the odd 00 block */
/***********************************/

/**
 * @class unitarity00oddP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the odd 00 block.
 */
class unitarity00oddP_Z2: public ThObservable {
public:

    /**
     * @brief unitarity00oddP_Z2 constructor.
     */
    unitarity00oddP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{00+}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};

/**
 * @class unitarity00oddM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the odd 00 block.
 */
class unitarity00oddM_Z2: public ThObservable {
public:

    /**
     * @brief unitarity00oddM_Z2 constructor.
     */
    unitarity00oddM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{00-}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};


/************************************/
/* Eigenvalues of the even 01 block */
/************************************/

/**
 * @class unitarity01eveP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the even 01 block.
 */
class unitarity01eveP_Z2: public ThObservable {
public:

    /**
     * @brief unitarity01eveP_Z2 constructor.
     */
    unitarity01eveP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{01+}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};

/**
 * @class unitarity01eveM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the even 01 block.
 */
class unitarity01eveM_Z2: public ThObservable {
public:

    /**
     * @brief unitarity01eveM_Z2 constructor.
     */
    unitarity01eveM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{01-}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};


/***********************************/
/* Eigenvalues of the odd 01 block */
/***********************************/

/**
 * @class unitarity01oddP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the odd 01 block.
 */
class unitarity01oddP_Z2: public ThObservable {
public:

    /**
     * @brief unitarity01oddP_Z2 constructor.
     */
    unitarity01oddP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{01+}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};

/**
 * @class unitarity01oddM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the odd 01 block.
 */
class unitarity01oddM_Z2: public ThObservable {
public:

    /**
     * @brief unitarity01oddM_Z2 constructor.
     */
    unitarity01oddM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{01-}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};


/**********************************/
/* Eigenvalue of the odd 10 block */
/**********************************/

/**
 * @class unitarity10odd_Z2
 * @ingroup GeneralTHDM
 * @brief Eigenvalue of the odd 10 block.
 */
class unitarity10odd_Z2: public ThObservable {
public:

    /**
     * @brief unitarity10odd_Z2 constructor.
     */
    unitarity10odd_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{10}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};


/************************************/
/* Eigenvalues of the even 11 block */
/************************************/

/**
 * @class unitarity11eveP_Z2
 * @ingroup GeneralTHDM
 * @brief Plus-sign eigenvalue of the even 11 block.
 */
class unitarity11eveP_Z2: public ThObservable {
public:

    /**
     * @brief unitarity11eveP_Z2 constructor.
     */
    unitarity11eveP_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{11+}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};

/**
 * @class unitarity11eveM_Z2
 * @ingroup GeneralTHDM
 * @brief Minus-sign eigenvalue of the even 11 block.
 */
class unitarity11eveM_Z2: public ThObservable {
public:

    /**
     * @brief unitarity11eveM_Z2 constructor.
     */
    unitarity11eveM_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{even}}_{11-}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};


/**********************************/
/* Eigenvalue of the odd 11 block */
/**********************************/

/**
 * @class unitarity11odd_Z2
 * @ingroup GeneralTHDM
 * @brief Eigenvalue of the odd 11 block.
 */
class unitarity11odd_Z2: public ThObservable {
public:

    /**
     * @brief unitarity11odd_Z2 constructor.
     */
    unitarity11odd_Z2(const StandardModel& SM_i);

    /**
     * @return @f$a^{\text{odd}}_{11}@f$
     */
    double computeThValue();

private:
    const GeneralTHDMZ2& myGTHDM;
};

#endif /* GENERALTHDMZ2UNITARITY_H */
