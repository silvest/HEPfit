/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSM_H
#define	EWSM_H

#include <cstring>
#include <gslpp.h>
#include "StandardModel.h"
#include "EWSMcache.h"
#include "EWSMOneLoopEW.h"
#include "EWSMTwoLoopQCD.h"
#include "EWSMThreeLoopQCD.h"
#include "EWSMTwoLoopEW.h"
#include "EWSMThreeLoopEW2QCD.h"
#include "EWSMThreeLoopEW.h"
#include "EWSMApproximateFormulae.h"
#include "EWSMOneLoopEW_HV.h"
#include "EWSMTwoFermionsLEP2.h"

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
    
    // The number of the parameters relevant to EW observables
    static const int NumSMParams = 21;
        
    
    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @brief EWSM constructor
     * @param[in] SM_i reference to a StandardModel object
     */
    EWSM(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return a pointer to the EWSMTwoFermionsLEP2 object
     */
    EWSMTwoFermionsLEP2* getMyTwoFermionsLEP2() const {
        return myTwoFermionsLEP2;
    }
    
    /**
     * @brief set a flag for M_W
     * @param[in] schemeMw NORESUM, OMSI, INTERMEDIATE, OMSII or APPROXIMATEFORMULA
     */
    void setSchemeMw(schemes_EW schemeMw) {
        this->schemeMw = schemeMw;
    }

    /**
     * @brief set a flag for rho_Z
     * @param[in] schemeRhoZ NORESUM, OMSI, INTERMEDIATE, OMSII or APPROXIMATEFORMULA
     * @attention NORESUM is not available, since reducible two-loop EW corrections have not been implemented yet. 
     */
    void setSchemeRhoZ(schemes_EW schemeRhoZ) {
        this->schemeRhoZ = schemeRhoZ;
    }
    
    /**
     * @brief set a flag for kappa_Z
     * @param[in] schemeKappaZ NORESUM, OMSI, INTERMEDIATE, OMSII or APPROXIMATEFORMULA
     */
    void setSchemeKappaZ(schemes_EW schemeKappaZ) {
        this->schemeKappaZ = schemeKappaZ;
    }

    
    //////////////////////////////////////////////////////////////////////// 

    bool checkSMparams(double Params_cache[]) const;
    
    
    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @param[in] s invariant mass squared 
     * @return the leptonic corrections to alpha
     */
    double DeltaAlphaLepton(const double s) const;    

    /**
     * @return the sum of the leptonic and hadronic corrections to alpha at Mz
     */
    double DeltaAlphaL5q() const;

    /**
     * @param[in] s invariant mass squared 
     * @return the top-quark corrections to alpha
     */
    double DeltaAlphaTop(const double s) const;    

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
     * @brief SM contribution to effective coupling rho_Z^l
     * @param[in] l name of lepton
     * @return rho_Z^l in the SM
     */
    complex rhoZ_l_SM(const StandardModel::lepton l) const;

    /**
     * @brief SM contribution to effective coupling rho_Z^q
     * @param[in] q name of quark
     * @return rho_Z^q in the SM
     */
    complex rhoZ_q_SM(const StandardModel::quark q) const;    
    
    /**
     * @brief SM contribution to effective coupling kappa_Z^l
     * @param[in] l name of lepton
     * @return kappa_Z^l in the SM
     */
    complex kappaZ_l_SM(const StandardModel::lepton l) const;

    /**
     * @brief SM contribution to effective coupling kappa_Z^q
     * @param[in] q name of quark
     * @return kappa_Z^q in the SM
     */
    complex kappaZ_q_SM(const StandardModel::quark q) const;    

    /**
     * @brief SM contribution to g_V^l
     * @param[in] l lepton
     * @return g_V^l
     */
   complex gVl_SM(const StandardModel::lepton l) const;

    /**
     * @brief SM contribution to g_V^q
     * @param[in] q quark
     * @return g_V^q
     */
    complex gVq_SM(const StandardModel::quark q) const;

    /**
     * @brief SM contribution to g_A^l
     * @param[in] l lepton
     * @return g_A^l
     */
    complex gAl_SM(const StandardModel::lepton l) const;

    /**
     * @brief SM contribution to g_A^q
     * @param[in] q quark
     * @return g_A^q
     */
    complex gAq_SM(const StandardModel::quark q) const;    
    
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

    /**
     * @brief R_b^0 with the complete fermionic EW two-loop corrections
     * @param[in] DeltaAlphaL5q_i the sum of the leptonic and hadronic corrections to alpha at Mz
     * @return R_b^0 in the SM, obtained from an approximate formula
     * @attention This function will be used if boolR0bApproximate is true. 
     */
    double R0_bottom_SM() const;    
    
    
    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return the top-quark corrections to the Z-b-bbar vertex
     */
    double taub() const;    
    

    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @param[in] l lepton
     * @return flavor-dependent correction to rho_Z^l with respect to that for the charged leptons
     */
    complex rhoZ_l_SM_FlavorDep(const StandardModel::lepton l) const;

    /**
     * @param[in] q quark
     * @return flavor-dependent correction to rho_Z^q with respect to that for the charged leptons
     */
    complex rhoZ_q_SM_FlavorDep(StandardModel::quark q) const;

    /**
     * @param[in] l lepton
     * @return flavor-dependent correction to kappa_Z^l with respect to that for the charged leptons
     */
    complex kappaZ_l_SM_FlavorDep(const StandardModel::lepton l) const;

    /**
     * @param[in] q quark
     * @return flavor-dependent correction to kappa_Z^q with respect to that for the charged leptons
     */
    complex kappaZ_q_SM_FlavorDep(StandardModel::quark q) const;
    
    
    ////////////////////////////////////////////////////////////////////////     
protected:

    const StandardModel& SM;


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
    
    EWSMTwoFermionsLEP2* myTwoFermionsLEP2;
        
    // accuracy in the iterative calculation of Mw
    static const double Mw_error;
    
    
    ////////////////////////////////////////////////////////////////////////     
    //caches
    
    bool bUseCacheEWSM; // true for caching
    
    mutable double DeltaAlphaLepton_params_cache[NumSMParams];
    mutable double DeltaAlphaLepton_cache;
    
    mutable double DeltaAlpha_params_cache[NumSMParams];
    mutable double DeltaAlpha_cache;
    
    mutable double Mw_params_cache[NumSMParams];
    mutable double Mw_cache;

    mutable double rhoZ_l_params_cache[6][NumSMParams];
    mutable complex rhoZ_l_cache[6];

    mutable double rhoZ_q_params_cache[6][NumSMParams];
    mutable complex rhoZ_q_cache[6];

    mutable double kappaZ_l_params_cache[6][NumSMParams];
    mutable complex kappaZ_l_cache[6];

    mutable double kappaZ_q_params_cache[6][NumSMParams];
    mutable complex kappaZ_q_cache[6];
    
    mutable double GammaW_params_cache[NumSMParams];
    mutable double GammaW_cache;
    
    mutable double R0b_params_cache[NumSMParams];
    mutable double R0b_cache;
        
    
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
     * @param[in] Mw_i the W boson mass
     * @param[in] DeltaRho
     * @param[in] DeltaR_rem
     * @return resummed Mw
     */
    double resumMw(const double Mw_i, const double DeltaRho[orders_EW_size],
                   const double DeltaR_rem[orders_EW_size]) const;
    
    /**
     * @param[in] DeltaRho
     * @param[in] deltaRho_rem 
     * @param[in] DeltaRbar_rem
     * @param[in] bool_Zbb true for f=b
     * @return resummed Re[rho_Z^f]
     */
    double resumRhoZ(const double DeltaRho[orders_EW_size],
                     const double deltaRho_rem[orders_EW_size], 
                     const double DeltaRbar_rem, const bool bool_Zbb) const;
    
    /**
     * @param[in] DeltaRho
     * @param[in] deltaKappa_rem 
     * @param[in] DeltaRbar_rem
     * @param[in] bool_Zbb true for f=b
     * @return resummed Re[kappa_Z^f]
     */
    double resumKappaZ(const double DeltaRho[orders_EW_size],
                       const double deltaKappa_rem[orders_EW_size],
                       const double DeltaRbar_rem, const bool bool_Zbb) const;

    
};

#endif	/* EWSM_H */

