/* 
 * File:   EWSM.h
 * Author: mishima
 */

#ifndef EWSM_H
#define	EWSM_H

#include <gslpp.h>
#include "StandardModel.h"
#include "EWModel.h"
#include "EWSMcache.h"
#include "EWSMOneLoopEW.h"
#include "EWSMTwoLoopQCD.h"
#include "EWSMThreeLoopQCD.h"
#include "EWSMTwoLoopEW.h"
#include "EWSMThreeLoopEW2QCD.h"
#include "EWSMThreeLoopEW.h"
#include "EWSMApproximateFormulae.h"
using namespace gslpp;


class EWSM : public EWModel {
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

    
    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return a reference to the EWSMcache object
     */
    const EWSMcache* getMyCache() const {
        return myCache;
    }
    
    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @return the leptonic corrections to alpha at Mz
     */
    double DeltaAlphaLepton() const;    

    /**
     * @return the sum of the leptonic and hadronic corrections to alpha at Mz
     */
    double DeltaAlphaL5q() const;
    
    /**
     * @return the total (leptonic+hadronic+top) corrections to alpha at Mz
     */
    double DeltaAlpha() const;
    
    /**
     * @brief electromagnetic coupling alpha at Mz
     * @return alpha(Mz)
     */
    double alphaMz() const;

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return the W boson mass
     */
    virtual double Mw() const;

    /**
     * @return Mw^2/Mz^2
     */
    virtual double cW2() const;
    
    /**
     * @return 1-Mw^2/Mz^2
     */
    virtual double sW2() const;
    
    /**
     * @brief effective coupling rho_Z^l
     * @param[in] l lepton
     * @return rho_Z^l
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;
    
    /**
     * @brief effective coupling rho_Z^q
     * @param[in] q quark
     * @return rho_Z^q
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;
    
    /**
     * @brief the ratio of the effective couplings for neutral-current interactions
     * @param[in] l lepton
     * @return g_V^l/g_A^l
     */
    virtual complex gZl_over_gAl(const StandardModel::lepton l) const;

    /**
     * @brief the ratio of the effective couplings for neutral-current interactions
     * @param[in] q quark
     * @return g_V^q/g_A^q
     */
    virtual complex gZq_over_gAq(const StandardModel::quark q) const;
    
    /**
     * @return the total width of the W boson
     */
    virtual double GammaW() const;    
    
    /**
     * @return NP contribution to oblique parameter S
     */
    virtual double obliqueS() const {
        return 0.0;
    };
        
    /**
     * @return NP contribution to oblique parameter T
     */
    virtual double obliqueT() const {
        return 0.0;
    };
    
    /**
     * @return NP contribution to oblique parameter U
     */
    virtual double obliqueU() const {
        return 0.0;
    };
    
    
    ////////////////////////////////////////////////////////////////////////     
protected:

    const StandardModel& SM;

            
    ////////////////////////////////////////////////////////////////////////           
        
    /**
     * @return the W boson mass in the SM
     */
    double Mw_SM() const;
    
    /** 
     * @brief computes Delta r from Mw()
     * @return Delta r in the SM
     */
    double DeltaR_SM() const;
    
    /**
     * @return Mw^2/Mz^2 in the SM
     */
    double cW2_SM() const;
    
    /**
     * @return 1-Mw^2/Mz^2 in the SM
     */
    double sW2_SM() const;
    
    /**
     * @brief SM contribution to effective coupling rho_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @return rho_Z^f in the SM
     */
    template<typename T> complex rhoZ_f_SM(const T f) const;

    /**
     * @brief SM contribution to effective coupling kappa_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @return kappa_Z^f in the SM
     */
    template<typename T> complex kappaZ_f_SM(const T f) const;
    
    /**
     * @param[in] li name of a neutrino
     * @param[in] lj name of a charged lepton
     * @return rho_ij^W for Gamma_W in the SM
     * @attention Fermion masses are neglected. 
     */    
    double rho_GammaW_l_SM(const StandardModel::lepton li, 
                           const StandardModel::lepton lj) const;
    
    /**
     * @param[in] qi name of a up-type quark
     * @param[in] qj name of a down-type quark
     * @return rho_ij^W for Gamma_W in the SM
     * @attention Fermion masses are neglected. 
     */    
    double rho_GammaW_q_SM(const StandardModel::quark qi, 
                           const StandardModel::quark qj) const;
    
    /**
     * @param[in] li name of a neutrino
     * @param[in] lj name of a charged lepton
     * @return the partial width of W^+ decay into an l_i\bar{l_j} pair in the SM
     * @attention Mixings in the lepton sector are neglected. 
     * @attention Fermion masses are neglected. 
     */
    double GammaW_l_SM(const StandardModel::lepton li, 
                       const StandardModel::lepton lj) const;

    /**
     * @param[in] qi name of a up-type quark
     * @param[in] qj name of a down-type quark
     * @return the partial width of W^+ decay into an q_i\bar{q_j} pair in the SM
     * @attention Fermion masses are neglected. 
     */
    double GammaW_q_SM(const StandardModel::quark qi, 
                       const StandardModel::quark qj) const;    
    
    /**
     * @return the total width of the W boson in the SM
     */
    double GammaW_SM() const;   
    

    ////////////////////////////////////////////////////////////////////////     
    // The functions below are used in NP models with S, T and U parameters. 
    
    /**
     * @param[in] S oblique parameter
     * @param[in] T oblique parameter
     * @param[in] U oblique parameter
     * @return the W-boson mass in a NP model from oblique parameters
     */
    double Mw_NP_fromSTU(const double S, const double T, const double U);
    
    /**
     * @brief effective coupling rho_Z^f
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @param[in] S oblique parameter
     * @param[in] T oblique parameter
     * @return rho_Z^f in a NP model from oblique parameters
     */
    template<typename Type> 
    complex rhoZ_f_NP_fromSTU(const Type f, const double T);

    /**
     * @brief the ratio of the effective couplings for neutral-current interactions
     * @param[in] f StandardModel::quark or StandardModel::lepton 
     * @param[in] S oblique parameter
     * @param[in] T oblique parameter
     * @return g_V^f/g_A^f in a NP model from oblique parameters
     */
    template<typename Type> 
    complex gZf_over_gAf_NP_fromSTU(const Type f,
                                    const double S, const double T);
    
    /**
     * @param[in] S oblique parameter
     * @param[in] T oblique parameter
     * @param[in] U oblique parameter
     * @return the total width of the W boson in a NP model from oblique parameters
     */    
    double GammaW_NP_fromSTU(const double S, const double T, const double U);    

    
    ////////////////////////////////////////////////////////////////////////         
private:
    
    bool flag_order[orders_EW_size]; 
    schemes_EW schemeMw, schemeRhoZ, schemeKappaZ;
    
    EWSMcache* myCache;
    EWSMOneLoopEW* myOneLoopEW;
    EWSMTwoLoopQCD* myTwoLoopQCD;
    EWSMThreeLoopQCD* myThreeLoopQCD;
    EWSMTwoLoopEW* myTwoLoopEW;
    EWSMThreeLoopEW2QCD* myThreeLoopEW2QCD;
    EWSMThreeLoopEW* myThreeLoopEW; 
    EWSMApproximateFormulae* myApproximateFormulae;
        
    // accuracy in the iterative calculation of Mw
    static const double Mw_error;
    
    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief computes Delta rho
     * @param[in] Mw_i the W boson mass
     * @param[out] DeltaRho[]
     */
    void ComputeDeltaRho(const double Mw_i, double DeltaRho[orders_EW_size]) const;  
     
    /**
     * @brief computes Delta r_rem
     * @param[in] Mw_i the W boson mass
     * @param[out] DeltaR_rem
     */
    void ComputeDeltaR_rem(const double Mw_i, double DeltaR_rem[orders_EW_size]) const;  
    
    /**
     * @param[in] DeltaRho
     * @param[in] DeltaR_rem
     * @return resummed Mw
     */
    double resumMw(const double DeltaRho[orders_EW_size],
                   const double DeltaR_rem[orders_EW_size]) const;
    
    /**
     * @param[in] DeltaRho
     * @param[in] deltaRho_rem 
     * @param[in] DeltaRbar_rem
     * @return resummed Re[rho_Z^f]
     */
    double resumRhoZ(const double DeltaRho[orders_EW_size],
                     const double deltaRho_rem[orders_EW_size], 
                     const double DeltaRbar_rem) const;
    
    /**
     * @param[in] DeltaRho
     * @param[in] deltaKappa_rem 
     * @param[in] DeltaRbar_rem
     * @return resummed Re[kappa_Z^f]
     */
    double resumKappaZ(const double DeltaRho[orders_EW_size],
                       const double deltaKappa_rem[orders_EW_size],
                       const double DeltaRbar_rem) const;    


};

#endif	/* EWSM_H */

