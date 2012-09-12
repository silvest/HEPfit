/* 
 * File:   EWSM.h
 * Author: mishima
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
#include "EWSMOneLoopLEP2.h"

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
     * @param[in] bDebug_i boolean value for debugging (true for debugging)
     */
    EWSM(const StandardModel& SM_i, const bool bDebug_i=false);

    
    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return a reference to the EWSMcache object
     */
    const EWSMcache* getMyCache() const {
        return myCache;
    }
    
    /**
     * @return boolean variable: if true, the approximate formula for R_b^0 is employed. 
     */
    bool isBoolR0bApproximate() const {
        return boolR0bApproximate;
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
    // LEP-II observables

    double dsigmaLEP2_l(const StandardModel::lepton l,const double s,const double Mw_i,
                        const double cos_theta,const double W,const double X,const double Y,
                        const double GammaZ) const;
    
    double dsigmaLEP2_q(const StandardModel::quark q,const double s,const double Mw_i,
                        const double cos_theta,const double W,const double X,const double Y,
                        const double GammaZ) const;

    
    ////////////////////////////////////////////////////////////////////////     
protected:

    const StandardModel& SM;


    ////////////////////////////////////////////////////////////////////////         
private:
    bool flag_order[orders_EW_size]; 
    schemes_EW schemeMw, schemeRhoZ, schemeKappaZ;
    bool boolR0bApproximate;
    
    EWSMcache* myCache;
    EWSMOneLoopEW* myOneLoopEW;
    EWSMTwoLoopQCD* myTwoLoopQCD;
    EWSMThreeLoopQCD* myThreeLoopQCD;
    EWSMTwoLoopEW* myTwoLoopEW;
    EWSMThreeLoopEW2QCD* myThreeLoopEW2QCD;
    EWSMThreeLoopEW* myThreeLoopEW; 
    EWSMApproximateFormulae* myApproximateFormulae;
    
    EWSMOneLoopLEP2* myLEP2;
        
    // accuracy in the iterative calculation of Mw
    static const double Mw_error;
    
    
    ////////////////////////////////////////////////////////////////////////     
    //caches
    
    bool bUseCacheEWSM; // true for caching
    
    // The number of the parameters relevant to EW observables
    static const int NumSMParams = 21;
        
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
    
    bool checkSMparams(double Params_cache[]) const {
        // 11 parameters in QCD:
        // "AlsMz","Mz","mup","mdown","mcharm","mstrange", "mtop","mbottom",
        // "mut","mub","muc"
        // 10 parameters in StandardModel
        // "GF", "ale", "dAle5Mz", "mHl", 
        // "mneutrino_1", "mneutrino_2", "mneutrino_3", "melectron", "mmu", "mtau"
        double SMparams[NumSMParams] = { 
            SM.getAlsMz(), SM.getMz(), SM.getGF(), SM.getAle(), SM.getDAle5Mz(),
            SM.getMHl(), SM.getMtpole(), 
            SM.getLeptons(SM.NEUTRINO_1).getMass(), 
            SM.getLeptons(SM.NEUTRINO_2).getMass(),
            SM.getLeptons(SM.NEUTRINO_3).getMass(),
            SM.getLeptons(SM.ELECTRON).getMass(),
            SM.getLeptons(SM.MU).getMass(),
            SM.getLeptons(SM.TAU).getMass(),
            SM.getQuarks(SM.UP).getMass(),
            SM.getQuarks(SM.DOWN).getMass(),
            SM.getQuarks(SM.CHARM).getMass(),
            SM.getQuarks(SM.STRANGE).getMass(),
            SM.getQuarks(SM.BOTTOM).getMass(),
            SM.getMut(), SM.getMub(), SM.getMuc()
        };
        
        // check updated parameters
        bool bCache = true;
        for(int i=0; i<NumSMParams; ++i) {
             if (Params_cache[i] != SMparams[i]) { 
                 Params_cache[i] = SMparams[i];                 
                 bCache &= false;
             }
        }
        
        return bCache;
    }
    
    
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

