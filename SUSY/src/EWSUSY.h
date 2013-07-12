/*
 * Copyright (C) 2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSUSY_H
#define	EWSUSY_H

#include <gslpp.h>
#include <PVfunctions.h>
#include "SUSY.h"

using namespace gslpp;

/**
 * @class EWSUSY
 * @ingroup SUSY
 * @brief A class for SUSY contributions to the EW precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class EWSUSY {
public:

    /**
     * @brief A constructor.
     * @param[in] SUSY_in A reference to a SUSY object. 
     */
    EWSUSY(const SUSY& SUSY_in);

    /**
     * @brief Sets parameters in Rosiek's notation. 
     */
    void SetRosiekParameters();

    /**
     * @brief Fermionic contribuiton to the transverse part of a gauge-boson self-energy.
     * @param[in] mu The renormalization scale. 
     * @param[in] p2
     * @param[in] mi
     * @param[in] mj
     * @param[in] cV_aij
     * @param[in] cV_bji
     * @param[in] cA_aij
     * @param[in] cA_bji
     * @return
     */
    complex FA(const double mu, const double p2, const double mi, const double mj,
               const double cV_aij, const double cV_bji,
               const double cA_aij, const double cA_bji) const;

    complex PiT_Z(const double mu, const double p2, const double Mw) const;


private:

    /* An object of PVdunctions class */
    const PVfunctions PV;

    /* A reference to the SUSY object passed to the constructor. */
    const SUSY& mySUSY;

    /* Yukawa couplings in Rosiek's notation */
    matrix<complex> Yu, Yd, Yl;

    /* Trilinear couplings in Rosiek's notation */
    matrix<complex> Au, Ad, Al;

    /* rotation matrices in Rosiek's notation */
    matrix<complex> Zm, Zp, Zn, ZU, ZD, ZL;

};

#endif	/* EWSUSY_H */

