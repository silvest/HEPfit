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
    /** @return the shift of the master function Y */
    double getdY() const { return dY; }
    /** @return the shift of the master function Z */
    double getdZ() const { return dZ; }
    /** @return the shift of C7eff at the scale muNP */
    double getdC7() const { return dC7; }
    /** @return the scale of the low-energy shifts dC7, dC9, dC10 */
    double getmuNP() const { return muNP; }


protected:

    virtual void setParameter(const std::string, const double&);
    mutable Matching<CMFVMatching,CMFV> CMFVM;

private:
    double Ftt;
    double dX, dY, dZ, dC7, muNP; ///< 

};

/**
 * @}
 */

#endif /* CMFV_H */

