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

    /**
     * @brief schemes for the resummations in Mw, rho_Z^f and kappa_Z^f
     * 
     * APPROXIMATEFORMULA for the use of approximate formulae
     */
    enum schemes_EW {NORESUM=0, OMSI, INTERMEDIATE, OMSII, APPROXIMATEFORMULA}; 

    
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
     * @return resummation scheme for Mw
     */
    schemes_EW getSchemeMw() const {
        return schemeMw;
    }

    /**
     * @return resummation scheme for rho_Z^f
     */
    schemes_EW getSchemeRhoZ() const {
        return schemeRhoZ;
    }

    /**
     * @return resummation scheme for kappa_Z^f
     */
    schemes_EW getSchemeKappaZ() const {
        return schemeKappaZ;
    }

    /**
     * @param[in] order
     * @return 
     */
    bool getFlag_order(orders_EW order) const {
        return flag_order[order];
    }    
    
    
    ////////////////////////////////////////////////////////////////////////     
    
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
    
    /**
     * @brief electromagnetic coupling alpha at Mz
     * @return alpha(Mz)
     */
    double getAlphaMz() const {
        return alphaMz;
    }   

    /**
     * @return the total radiative corrections to alpha at Mz
     */
    double getDeltaAlpha() const {
        return DeltaAlpha;
    }

    /**
     * @return the leptonic and 5-flavour hadronic corrections to alpha at Mz
     */
    double getDeltaAlpha_l5q() const {
        return DeltaAlpha_l5q;
    }    
    
    /**
     * @return the W boson mass at tree level
     */
    double getMw_tree() const {
        return Mw_tree;
    }   
    
    /**
     * @brief the W-boson mass
     * @return M_W
     */
    double getMw() const {
        return Mw;
    }    
    
    /**
     * @brief effective coupling rho_Z^l
     * @param[in] l name of a lepton 
     * @return rho_Z^l for lepton "l" 
     */
    complex getRhoZ_l(const StandardModel::lepton l) const {
        return rhoZ_l[l];
    }

    /**
     * @brief effective coupling rho_Z^q
     * @param[in] q name of a quark
     * @return rho_Z^q for quark "q" 
     */
    complex getRhoZ_q(const StandardModel::quark q) const {
        return rhoZ_q[q];
    }
    
    /**
     * @brief effective coupling kappa_Z^l
     * @param[in] l name of a lepton 
     * @return kappa_Z^l for lepton "l" 
     */
    complex getKappaZ_l(const StandardModel::lepton l) const {
        return kappaZ_l[l];
    }
    
    /**
     * @brief effective coupling kappa_Z^q
     * @param[in] q name of a quark
     * @return kappa_Z^q for quark "q" 
     */
    complex getKappaZ_q(const StandardModel::quark q) const {
        return kappaZ_q[q];
    }
    

    //////////////////////////////////////////////////////////////////////// 

    /**
     * @brief 
     * @param[in] schemeMw_i resummation scheme for Mw
     * @param[in] schemeRhoZ_i resummation scheme for rho_Z^f
     * @param[in] schemeKappaZ_i resummation scheme for kappa_Z^f
     * @param[in] flag_order_i
     */
    void setFlags(const schemes_EW schemeMw_i, 
                  const schemes_EW schemeRhoZ_i, const schemes_EW schemeKappaZ_i,
                  const bool flag_order_i[orders_EW_size]);
    
    /**
     * @brief computes alpha(Mz), the W-boson mass, and rho_Z^f and kappa_Z^f
     * 
     * EWSM::setFlags() has to be called in advance.
     */
    void Compute();
    
    
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
    ApproximateFormulae* myApproximateFormulae;
    
    /* flags */
    schemes_EW schemeMw, schemeRhoZ, schemeKappaZ;
    bool flag_order[orders_EW_size];
    
    /* The last elements store the total contribution.  */
    double DeltaAlpha_l[orders_EW_size+1];
    double DeltaAlpha_t[orders_EW_size+1];
    double DeltaRho[orders_EW_size+1];
    double DeltaR_rem[orders_EW_size+1];
    complex deltaRho_rem_l[6][orders_EW_size+1], deltaRho_rem_q[6][orders_EW_size+1];
    complex deltaKappa_rem_l[6][orders_EW_size+1], deltaKappa_rem_q[6][orders_EW_size+1];
    double DeltaRbar_rem;
    
    double alphaMz, DeltaAlpha, DeltaAlpha_l5q;    
    double Mw_tree, Mw;
    complex rhoZ_l[6], rhoZ_q[6];
    complex kappaZ_l[6], kappaZ_q[6];
    
    
   //////////////////////////////////////////////////////////////////////// 

    /**
     * @brief computes alpha(Mz)
     */
    void ComputeAlphaMz();    

    /**
     * @brief computes the W-boson mass
     * 
     * EWSM::AlphaMz() has to be called in advance. 
     */
    void ComputeMw();  
    
    /**
     * @brief computes rho_Z^f
     * 
     * EWSM::ComputeAlphaMz() and EWSM::CoumuteMw() have to be called in advance. 
     */
    void ComputeRhoZ();
    
     /**
     * @brief computes kappa_Z^f
     * 
     * EWSM::ComputeAlphaMz() and EWSM::CoumuteMw() have to be called in advance. 
     */
    void ComputeKappaZ();
    
    /**
     * @return resummed Mw
     */
    double resumMw();
    
    /**
     * @return resummed rho_Z^f
     */
    complex resumRhoZ(const complex deltaRho_rem[orders_EW_size+1]);
    
    /**
     * @return resummed kappa_Z^f
     */
    complex resumKappaZ(const complex deltaKappa_rem[orders_EW_size+1]);
    
    
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

