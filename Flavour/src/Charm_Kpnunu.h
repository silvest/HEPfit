/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHARM_KPNUNU_H
#define	CHARM_KPNUNU_H

#include "StandardModel.h"
#include <sstream>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_clausen.h>

/**
 * 
 * @class: Charm_Kpnunu
 * @details: class to compute the charm contribution to the process \f$  K^{0} \rightarrow \pi^{0}  
 * \nu \bar{\nu}\f$ at the NNLO in QCD corrections, according to hep-ph/0603079.
 */
class Charm_Kpnunu {

public:
    
    /**
     * 
     * @brief constructor
     */
    Charm_Kpnunu(const StandardModel& model_i);
    
    /**
     * 
     * @brief destructor
     */
    ~Charm_Kpnunu();
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return Wilson coefficients related to the Z-penguin contribution, given
     * at the renormalization scale \f$ \mu_{W} \f$
     */
    gslpp::vector<double> Cp(orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order
     * @param nf, number of flavours
     * @return LO, NLO and NNLO RG evolution matrix for the Z-penguin contribution
     */
    gslpp::matrix<double> RGevolP(int nf, orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order
     * @param nf, number of flavours
     * @return non trivial threshold matching at NNLO level for the Wilson coefficients 
     * related to the Z-penguin contribution
     */
    gslpp::vector<double> ThresholdCp(orders order);    
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return Wilson coefficients related to the Z-penguin contribution evolved
     * down to the renormalization scale \f$ mu_{c} \f$ 
     */
    gslpp::vector<double> C_p(orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return coefficient recasting the total Z-penguin contribution to BR of the process
     */
    double C_P(orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return Wilson coefficients related to the EW box contribution, given
     * at the renormalization scale \f$ \mu_{W} \f$
     */
    gslpp::vector<double> Cb(orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order
     * @param nf, number of flavours
     * @return LO, NLO and NNLO RG evolution matrix for the EW box contribution
     */
    gslpp::matrix<double> RGevolB(int nf, orders order);
    /**
     * 
     * @param order, QCD perturbation theory order
     * @return non trivial threshold matching at NNLO level for the Wilson coefficients 
     * related to the EW box contribution
     */
    gslpp::vector<double> ThresholdCb(orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return Wilson coefficients related to the EW box contribution evolved
     * down to the renormalization scale \f$ \mu_{c} \f$ 
     */
    gslpp::vector<double> C_b(orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return coefficient recasting the EW box contribution related to the light leptons
     *  to the BR of the process
     */
    double C_Be(orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return coefficient recasting the EW box contribution related to the tau lepton
     *  to the BR of the process
     */
    double C_Bt(orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return the phenomenological function P_C which contains the appropriate C_P, C_Be and C_Bt 
     *  linear combination appearing explicitly in the final BR formula of the process 
     */
    double P_C(orders order);
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return P_C + isospin correction + peculiar contribution of the top quark coming from the 
     * loop function X_t for this process (in respect with the one present also in 
     * \f$  K^{0} \rightarrow \pi^{0}  \nu \bar{\nu} \f$)
     */
    double C_TOT(orders order, orders_ew order_ew);
    
private:
    const StandardModel& model;
    const StandardModelMatching& modelmatching;
    gslpp::vector<double> cp, dcp, c_p, cpmuW0, cpmuW1, cpmuW2, cb, dcb, c_b, cbmuW0,
                   cbmuW1, cbmuW2;
    gslpp::matrix<double> U4p, U5p, J5p1, J4p1, J5p2, J4p2, dc_p, 
                   U4b, U5b, J5b1, J4b1, J5b2, J4b2, dc_b;
    double etab, etacb, etac, mc, kc ,xi1c, xi2c, xc, CP, CBe, CBt;
};

#endif	/* CHARM_KPNUNU_H */

