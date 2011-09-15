/* 
 * File:   EWSM.h
 * Author: mishima
 */

#ifndef EWSM_H
#define	EWSM_H

#include <StandardModel.h>
#include "EWSMcommon.h"
#include "OneLoopEW.h"
#include "TwoLoopQCD.h"
#include "ThreeLoopQCD.h"
#include "TwoLoopEW.h"
#include "ThreeLoopEW2QCD.h"
#include "ThreeLoopEW.h"
#include "ApproximateFormulae.h"

using namespace gslpp;


class EWSM {
public:
    
    /**
     * @brief the order of radiative corrections
     * 
     * The number of elements is set in "orders_EW_size".
     */
    enum orders_EW {EW1=0, EW1QCD1, EW1QCD2, EW2, EW2QCD1, EW3, orders_EW_size};    


    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @brief EWSM constructor
     * @param[in] SM_i reference to a StandardModel object
     */
    EWSM(const StandardModel& SM_i);

    /**
     * @brief EWSM copy constructor
     * @param[in] orig reference to an EWSM object
     */
    EWSM(const EWSM& orig);
    
    /**
     * @brief EWSM destructor
     */
    virtual ~EWSM();

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @return a pointer to the EWSMcommon object in EWSM class
     */
    EWSMcommon* getEWSMC() const {
        return EWSMC;
    }

    /**
     * @brief leptonic contribution to alpha
     * @param[in] order the order of the contribution
     * @return Delta alpha_{lept} at each order
     */
    double getDeltaAlpha_l(const orders_EW order) const {
        return DeltaAlpha_l[order];
    }

    /**
     * @brief top-quark contribution to alpha
     * @param[in] order the order of the contribution
     * @return Delta alpha_{top} at each order
     */
    double getDeltaAlpha_t(const orders_EW order) const {
        return DeltaAlpha_t[order];
    }    
    
    /**
     * @brief leading contribution to Delta r
     * @param[in] order the order of the contribution
     * @return Delta rho at each order
     */
    double getDeltaRho(const orders_EW order) const {
        return DeltaRho[order];
    }

    /**
     * @brief remainder contribution to Delta r
     * @param[in] order the order of the contribution
     * @return Delta r_{rem} at each order
     */
    double getDeltaR_rem(const orders_EW order) const {
        return DeltaR_rem[order];
    }
    
    /**
     * @brief remainder contribution for rho_Z^f and kappa_Z^f
     * @return Delta rbar_{rem} at each order
     */
    double getDeltaRbar_rem() const {
        return DeltaRbar_rem;
    }

    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of a lepton 
     * @param[in] order the order of the contribution
     * @return delta rho_{rem}^l at each order
     */        
    complex getDeltaRho_rem_l(const StandardModel::lepton l,
                              const orders_EW order) const {
        return deltaRho_rem_l[l][order];
    }

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of a quark
     * @param[in] order the order of the contribution 
     * @return delta rho_{rem}^q at each order
     */
    complex getDeltaRho_rem_q(const StandardModel::quark q,
                              const orders_EW order) const {
        return deltaRho_rem_q[q][order];
    }
    
    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of a lepton 
     * @param[in] order the order of the contribution
     * @return delta kappa_{rem}^l at each order
     */
    complex getDeltaKappa_rem_l(const StandardModel::lepton l,
                                const orders_EW order) const {
        return deltaKappa_rem_l[l][order];
    }

    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of a quark
     * @param[in] order the order of the contribution
     * @return delta kappa_{rem}^q at each order
     */
    complex getDeltaKappa_rem_q(const StandardModel::quark q,
                                const orders_EW order) const {
        return deltaKappa_rem_q[q][order];
    }
    

    //////////////////////////////////////////////////////////////////////// 

    /**
     * @brief computes radiative corrections to alpha
     * @param[in] flag_order_i
     */
    void ComputeDeltaAlpha(const bool flag_order_i[orders_EW_size]);    

    /**
     * @brief computes radiative corrections to Charged-Current interactions
     * @param[in] Mw_i
     * @param[in] flag_order
     */
    void ComputeCC(double Mw_i, const bool flag_order[orders_EW_size]);  
    
    /**
     * @brief computes radiative corrections to Neutral-Current interactions 
     * @param[in] Mw_i
     * @param[in] flag_order
     */
    void ComputeNC(double Mw_i, const bool flag_order[orders_EW_size]);
     
    
    ////////////////////////////////////////////////////////////////////////     
    
protected:
    const StandardModel& SM;

    
private:
    EWSMcommon* EWSMC;
    OneLoopEW* myOneLoopEW;
    TwoLoopQCD* myTwoLoopQCD;
    ThreeLoopQCD* myThreeLoopQCD;
    TwoLoopEW* myTwoLoopEW;
    ThreeLoopEW2QCD* myThreeLoopEW2QCD;
    ThreeLoopEW* myThreeLoopEW; 
    
    double DeltaAlpha_l[orders_EW_size];
    double DeltaAlpha_t[orders_EW_size];
    double DeltaRho[orders_EW_size];
    double DeltaR_rem[orders_EW_size];
    complex deltaRho_rem_l[6][orders_EW_size], deltaRho_rem_q[6][orders_EW_size];
    complex deltaKappa_rem_l[6][orders_EW_size], deltaKappa_rem_q[6][orders_EW_size];
    double DeltaRbar_rem;

    
    ////////////////////////////////////////////////////////////////////////
    

//    /**
//     * @return the charm-quak mass at Mz, mc(Mz)
//     */
//    double mcMz() const;
//
//    /**
//     * @return the bottom-quak mass at Mz, mb(Mz)
//     */
//    double mbMz() const;
//



};

#endif	/* EWSM_H */

