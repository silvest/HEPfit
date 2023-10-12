/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHARM_KPNUNU_H
#define CHARM_KPNUNU_H

class StandardModel;
#include "OrderScheme.h"
#include "gslpp.h"
#include "WilsonCoefficient.h"
#include <vector>
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
     * @return coefficient recasting the total Z-penguin contribution to BR of the process
     */
    double C_P(orders order);

    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return coefficient recasting the QED Z-penguin contribution to BR of the process
     */
    double C_P_qed(orders_qed order_qed);

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
     * @brief This method solves the RGE equation: d/dlog(mu) C = gamma_transp C 
     * from mu_w to mu_c relevant for the box and penguin contributions in Kpnunu decay.
     * References: 0805.4119v1 and hep-ph/0603079v3
     * @param order, order at which final Wilson coefficients are returned
     * @param order_qed, qed order at which final Wilson coefficients are returned
     * @param contribution, 0 for Box (dim=2) and 1 for Penguin (dim=3)
     * @return The Wilson coefficient array evolved to \mu_c at specified order. 
     * For box is: ( 4*CW^2 , CnuB ) and for penguin is: ( 4*Cp*CA , 4*Cm*CA , CnuP ) for qed and for qcd the same with CW=1 and CA=1.
     */
    gslpp::vector<double> C(orders order,orders_qed order_qed, int contribution);
   
    /**
     * @brief This method returns the J matrices useful for the RGE of the charm penguin/box diagrams.
     * @param order, order at which final J is returned
     * @param contribution, 0 for Box (dim=2) and 1 for Penguin (dim=3)
     * @return J matrices useful for the RGE of the charm penguin/box diagrams
     */
    gslpp::matrix<double> RGevol_J(orders order, int nf, int contribution);
    
    /**
     * @brief This method returns the R matrices useful for the RGE of the charm penguin/box diagrams in qed.
     * @param order_qed, qed order at which final R is returned
     * @param contribution, 0 for Box (dim=2) and 1 for Penguin (dim=3)
     * @return R at specified order, where R has to be considered multiplied like this: (alpha/4pi)/alphaS(muC)*R at O(1/alphaS) or LO_QED here and (alpha/4pi)*R at O(1) or NLO_QED11 here 
     */
    gslpp::matrix<double> RGevol_R(orders_qed order_qed, int nf, int contribution);


    /**
     * @brief 0805.4119v1 and hep-ph/0603079v3. To retrieve QCD CW at MuW just put NO_QED to order_qed and the qcd_order chosen to "order". To retrieve QED CW at MuW just put LO to order and the qed_order chosen to "qed_order"
     * @param order_qed, QED order at which initial conditions are returned 
     * @param order, QCD order at which initial conditions are returned 
     * @param contribution, 0 for Box (dim=2) and 1 for Penguin (dim=3)
     * @return The Wilson coefficient array  at \f$ \mu_{W} \f$ . For box is: ( 4*CW^2 , CnuB ) and for penguin is: ( 4*Cp*CA , 4*Cm*CA , CnuP ) for qed and for qcd the same with CW=1 and CA=1.
     */
    gslpp::vector<double> CWin_muw(orders order,orders_qed order_qed, int contribution);


    /**
     * 
     * @brief Reference: 0805.4119v1 and hep-ph/0603079v3. To retrieve QCD ADM just put NO_QED to order_qed and the qcd_order chosen to "order". To retrieve QED ADM just put LO to order and the qed_order chosen to "qed_order"
     * @param qed_order, QED perturbation theory order
     * @param order, QCD perturbation theory order
     * @param nf, Number of total active flavours
     * @param contribution, 0 for Box (dim=2) and 1 for Penguin (dim=3)
     * @return The TRANSPOSED Anomalous Dimension Matrix with nf flavours at specified order in QED for kpnunu decay.
     */
    gslpp::matrix<double> ADM(orders order,orders_qed order_qed , double nf, int contribution);


    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return coefficient recasting the EW box contribution related to the light leptons
     *  to the BR of the process
     */
    double C_Be(orders order);

    /**
     * 
     * @param order_qed, QED perturbation theory order 
     * @return coefficient recasting the QED box contribution related to the light leptons
     *  to the BR of the process
     */
    double C_Be_qed(orders_qed order_qed);

    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return coefficient recasting the EW box contribution related to the tau lepton
     *  to the BR of the process
     */
    double C_Bt(orders order);

    /**
     * 
     * @param order_qed, QED perturbation theory order 
     * @return coefficient recasting the QED box contribution related to the tau lepton
     *  to the BR of the process
     */
    double C_Bt_qed(orders_qed order_qed);

    /**
     * 
     * @param order, QCD perturbation theory order 
     * @return the phenomenological function P_C which contains the appropriate C_P, C_Be and C_Bt 
     *  linear combination appearing explicitly in the final BR formula of the process 
     */
    std::vector<WilsonCoefficient>& EVOCkpnn();
    
    
    /**
     * 
     * @param order, QCD perturbation theory order 
     * @param order_qed, QED perturbation theory order
     * @return P_C + isospin correction + peculiar contribution of the top quark coming from the 
     * loop function X_t for this process (in respect with the one present also in 
     * \f$  K^{0} \rightarrow \pi^{0}  \nu \bar{\nu} \f$)
     */
    double C_TOT(orders order, orders_qed order_qed);

private:
    const StandardModel& model;
    WilsonCoefficient evoCkpnn;
    gslpp::vector<double>  dcp,dcb;
    double mc_mc,etab,etacb,etac,kc,xc_mc_qed,L;
    double xi1c,xi2c,xice,xices;
    //box
    gslpp::vector<double> CW0b,CW1b,CW2b,CWeb,CWesb;
    gslpp::matrix<double> U0_4b,U0_5b,J1_4b,J2_4b,J1_5b,J2_5b,R0_4b,R0_5b,R1_4b,R1_5b;
    //penguin
    gslpp::vector<double> CW0p,CW1p,CW2p,CWep,CWesp;
    gslpp::matrix<double> U0_4p,U0_5p,J1_4p,J2_4p,J1_5p,J2_5p,R0_4p,R0_5p,R1_4p,R1_5p;
    
protected:
    std::vector<WilsonCoefficient> vevoCkpnn;
    
    
};

#endif /* CHARM_KPNUNU_H */

