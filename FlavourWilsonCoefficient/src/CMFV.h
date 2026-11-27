/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CMFV_H
#define CMFV_H

#include "StandardModel.h"
#include "CMFVMatching.h"

/**
 * @addtogroup FlavourWilsonCoefficent
 * @brief A module for EFT models for flavour observaqbles.
 * @details This module has implementations of several EFTs for flavour
 * analyses like the Weak Effective Theory at low energy and generalized
 * EFT models for @f$\Delta F = 1@f$ and @f$\Delta F = 2@f$.
 * @{
 */

/**
 * @class CMFV
 * @ingroup FlavourWilsonCoefficent
 * @brief %Model for CMFV contributions to flavour.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details New physics enters through real, flavour-universal shifts of the
 * gauge-independent master functions:
 * @li Ftt: shift of the top-top box function @f$S@f$ (all @f$\Delta F=2@f$ sectors),
 *   implemented by shifting @f$S_0(x_t,x_t)@f$ in the matching;
 * @li CMFV_dX, CMFV_dY, CMFV_dZ: shifts of @f$X@f$, @f$Y@f$, @f$Z@f$ ,
 *   added at leading order with the SM normalisation to @f$s\to d\nu\bar\nu@f$,
 *   @f$b\to q\nu\bar\nu@f$, @f$B_q\to\mu\mu@f$ and @f$b\to s\ell\ell@f$, with
 *   @f$\Delta C_{10}=-\Delta Y/s_W^2@f$ and @f$\Delta C_9=\Delta Y/s_W^2-4\Delta Z@f$;
 * @li CMFV_dC7: shift of @f$C_7^{\rm eff}@f$ at the scale CMFV_muNP .
 *
 * Optional parameters CMFV_r_bsg, CMFV_r_bsll, CMFV_r_Bmumu (default 1): multiplicative theory-uncertainty factors,
 * applied in the configuration as "CMFV_r_bsg*BR_bsgamma" etc. (see ScaledObs).
 *
 * Model flag CMFV_ZPenguinDominance (default false): @f$Z^0@f$-penguin dominance,
 * @f$\Delta X=\Delta Y=\Delta Z=\Delta C@f$, with @f$\Delta C@f$ given by CMFV_dX; the parameters CMFV_dY and
 * CMFV_dZ are then ignored.
 */

class CMFV : public StandardModel {
public:

    static const int NCMFVvars = 6;

    static const std::string CMFVvars[NCMFVvars];

    /**
     * @brief CMFV constructor
     */
    CMFV();

    /**
     * @brief A method to check if all the mandatory parameters for %CMFV
     * have been provided in model initialization.
     * @param[in] DPars a map of the parameters that are being updated in the Monte Carlo run
     * (including parameters that are varied and those that are held constant)
     * @return a boolean that is true if the execution is successful
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

    /**
     * @brief A get method to access the member reference of type %CMFVMatching
     * @return a reference to a %CMFVMatching object
     */
    virtual CMFVMatching& getMatching() const
    {
        return CMFVM.getObj();
    }

    /**
     *
     * @return Ftt
     */
    double getFtt() const
    {
        return Ftt;
    }

    /** @return the shift of the master function X */
    double getdX() const { return dX; }
    /** @return the shift of the master function Y (CMFV_dX with the flag CMFV_ZPenguinDominance) */
    double getdY() const { return FlagZPenguinDominance ? dX : dY; }
    /** @return the shift of the master function Z (CMFV_dX with the flag CMFV_ZPenguinDominance) */
    double getdZ() const { return FlagZPenguinDominance ? dX : dZ; }
    /** @return true with the flag CMFV_ZPenguinDominance: dX = dY = dZ = dC, CMFV_dY and CMFV_dZ ignored */
    bool getFlagZPenguinDominance() const { return FlagZPenguinDominance; }
    /** @return the shift of C7eff at the scale muNP */
    double getdC7() const { return dC7; }
    /** @return the scale of the low-energy shifts dC7, dC9, dC10 */
    double getmuNP() const { return muNP; }


    /**
     * @brief Model flags: CMFV_ZPenguinDominance (default false), dX = dY = dZ = dC given by CMFV_dX.
     */
    virtual bool setFlag(const std::string name, const bool value);

protected:

    virtual void setParameter(const std::string, const double&);
    mutable Matching<CMFVMatching,CMFV> CMFVM;

private:
    double Ftt;
    double dX, dY, dZ, dC7, muNP; ///< 
    bool FlagZPenguinDominance = false; ///< Z0-penguin dominance: dY and dZ follow dX
    double r_bsg = 1., r_bsll = 1., r_Bmumu = 1.; ///< theory-uncertainty factors (optional parameters)

};

/**
 * @}
 */

#endif /* CMFV_H */

