/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMWUNITARITY_H
#define	THDMWUNITARITY_H

#include "ThObservable.h"
#include "THDMW.h"
#include <gslpp.h>

/**
 * @class THDMWunitarity
 * @ingroup THDMW
 * @brief An observable class for the requirement of perturbative unitarity.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to require the unitarity for all the tree level 
 * scalar-scalar scattering amplitudes.
 * The eigenvalues of the S-matrix can be found in @cite ?.
 * They should be smaller than ? in magnitude to preserve the unitarity of the S-matrix.
 */
class THDMWunitarity {
public:
    /**
     * @brief Constructor.
     * @param[in] ?
     */
   THDMWunitarity(const StandardModel& SM_i);

    /**
     * @brief Destructor.
     */
    virtual ~THDMWunitarity();

    /**
     * @brief Computes the eigenvalues of the S matrix
     */
    bool CalcSeigen(gslpp::matrix<gslpp::complex>& Seigenvectors_i, gslpp::vector<double>& Seigenvalues_i);

    /**
     * @brief Assigns to a vector the eigenvalues of the S-matrix
     */
    gslpp::vector<double> getSeigenvalues();
    
    
private:
    const THDMW * myTHDMW;
    
    gslpp::matrix<gslpp::complex> Smatrix; 
    gslpp::matrix<gslpp::complex> Seigenvectors;
    gslpp::vector<double> Seigenvalues;
};


class THDMWunitarity1: public ThObservable {
public:

    /**
     * @brief THDMWunitarity1 constructor.
     */
    THDMWunitarity1(const StandardModel& SM_i);

    /**
     * @return
     */
    double computeThValue();

private:
    THDMWunitarity myTHDMWunitarity;
};

#endif	/* THDMWUNITARITY_H */
