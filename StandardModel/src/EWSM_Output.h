/*
 * Copyright (C) 2013 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSM_OUTPUT_H
#define	EWSM_OUTPUT_H

#include "StandardModel.h"

/**
 * @class EWSM_Output
 * @ingroup StandardModel
 * @brief A class for testing SM radiative corrections to the %EW precision
 * obsrvables. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class EWSM_Output : public StandardModel {
public:

    /**
     * @brief Constructor.
     * @param[in] EWSM_in a reference to an object of type EWSM
     */
    EWSM_Output(const StandardModel& SM_in);


    ////////////////////////////////////////////////////////////////////////

    void outputEachDeltaR(const double Mw_i) const;

    void outputEachDeltaRhoZ_l(const StandardModel::lepton l, const double Mw_i) const;

    void outputEachDeltaRhoZ_q(const QCD::quark q, const double Mw_i) const;

    void outputEachDeltaRhoZ(const double f_AlphaToGF,
            const double DeltaRho[StandardModel::orders_EW_size],
            const double deltaRho_rem[StandardModel::orders_EW_size],
            const double DeltaRbar_rem, const bool bool_Zbb,
            const double taub[StandardModel::orders_EW_size],
            const double ZbbSubtract) const;

    void outputEachDeltaKappaZ_l(const StandardModel::lepton l, const double Mw_i) const;

    void outputEachDeltaKappaZ_q(const QCD::quark q, const double Mw_i) const;

    void outputEachDeltaKappaZ(const double f_AlphaToGF,
            const double cW2overSW2,
            const double DeltaRho[StandardModel::orders_EW_size],
            const double deltaKappa_rem[StandardModel::orders_EW_size],
            const double DeltaRbar_rem, const bool bool_Zbb,
            const double taub[StandardModel::orders_EW_size],
            const double ZbbSubtract,
            const double Zgamma_EW2) const;


    ////////////////////////////////////////////////////////////////////////
private:

};

#endif	/* EWSM_OUTPUT_H */

