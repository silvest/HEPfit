/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWUNITARITY_H
#define	THDMWUNITARITY_H

#include "ThObservable.h"

class THDMW;
class THDMWcache;

/**
 * @class THDMWunitarityLO
 * @ingroup THDMW
 * @brief An observable class for the requirement of perturbative unitarity at leading order.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to require unitarity for all the tree level 
 * scalar-scalar scattering amplitudes.
 * The eigenvalues of the S-matrix can be found in ?.
 * They should be smaller than ? in magnitude to preserve the unitarity of the S-matrix.
 */
class THDMWunitarityLO: public ThObservable {
public:

    /**
     * @brief THDMWunitarityLO constructor.
     */
    THDMWunitarityLO(const StandardModel& SM_i, unsigned int index_i);

    /**
     * @brief Destructor.
     */
    virtual ~THDMWunitarityLO();

    /**
     * @return Unitarity eigenvalues
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
    unsigned int index;
};

/**
 * @class THDMWunitarityNLO
 * @ingroup THDMW
 * @brief An observable class for the requirement of perturbative unitarity at next-to-leading order.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to require NLO unitarity for all the tree level 
 * scalar-scalar scattering amplitudes.
 * The eigenvalues of the S-matrix can be found in ?.
 * They should be smaller than ? in magnitude to preserve the unitarity of the S-matrix.
 */
class THDMWunitarityNLO: public ThObservable {
public:

    /**
     * @brief THDMWunitarityNLO constructor.
     */
    THDMWunitarityNLO(const StandardModel& SM_i, unsigned int index_i);

    /**
     * @brief Destructor.
     */
    virtual ~THDMWunitarityNLO();

    /**
     * @return NLO unitarity conditions
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
    unsigned int index;
};

/**
 * @class THDMWunitarityNLOp
 * @ingroup THDMW
 * @brief Another observable class for the requirement of perturbative unitarity at next-to-leading order.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class THDMWunitarityNLOp: public ThObservable {
public:

    /**
     * @brief THDMWunitarityNLOp constructor.
     */
    THDMWunitarityNLOp(const StandardModel& SM_i, unsigned int index_i);

    /**
     * @brief Destructor.
     */
    virtual ~THDMWunitarityNLOp();

    /**
     * @return Real part of the NLO unitarity eigenvalues
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
    unsigned int index;
};

/**
 * @class THDMWunitarityRp
 * @ingroup THDMW
 * @brief An observable class for the requirement of perturbativity of the unitarity criteria.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class THDMWunitarityRp: public ThObservable {
public:

    /**
     * @brief THDMWunitarityRp constructor.
     */
    THDMWunitarityRp(const StandardModel& SM_i, unsigned int index_i);

    /**
     * @brief Destructor.
     */
    virtual ~THDMWunitarityRp();

    /**
     * @return Rp ratios
     */
    double computeThValue();
private:
    const THDMW& myTHDMW;
    unsigned int index;
};

#endif	/* THDMWUNITARITY_H */
