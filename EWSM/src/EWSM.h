/* 　
 * File:   EWSM.h
 * Author: mishima
 */

#ifndef EWSM_H
#define	EWSM_H

#include <StandardModel.h>
#include "EWModel.h"
#include "EWSMcommon.h"
#include "OneLoopEW.h"
#include "TwoLoopQCD.h"
#include "ThreeLoopQCD.h"
#include "TwoLoopEW.h"
#include "ThreeLoopEW2QCD.h"
#include "ThreeLoopEW.h"
#include "ApproximateFormulae.h"

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
     * @return the leptonic corrections to alpha at Mz
     */
    double DeltaAlphaLepton() const;    

    /**
     * @return the sum of the leptonic and hadronic corrections to alpha at Mz
     */
    double DeltaAlphaL5q() const;
    
    /**
     * @return the total radiative corrections to alpha at Mz
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
    virtual double Mw();

    /**
     * @return Mw^2/Mz^2
     */
    virtual double cW2();
    
    /**
     * @return 1-Mw^2/Mz^2
     */
    virtual double sW2();
    
    /**
     * @brief effective coupling rho_Z^l
     * @param[in] l name of a lepton 
     * @return rho_Z^l
     */
    virtual complex rhoZ_l(const StandardModel::lepton l);

    /**
     * @brief 
     * @param[in] q name of a quark
     * @return rho_Z^q
     */
    virtual complex rhoZ_q(const StandardModel::quark q);
    
    /**
     * @brief the ratio of the effective couplings for neutral-current interactions
     * @param[in] l name of a lepton 
     * @return g_V^l/g_A^l
     */
    virtual complex gZl_over_gAl(const StandardModel::lepton l);
    
    /**
     * @brief the ratio of the effective couplings for neutral-current interactions
     * @param[in] q name of a quark
     * @return g_V^q/g_A^q
     */
    virtual complex gZq_over_gAq(const StandardModel::quark q);

    /**
     * @return the total width of the W boson
     */
    virtual double GammaW();    
    
    
    ////////////////////////////////////////////////////////////////////////     
protected:
    
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
     * @param[in] l name of a lepton 
     * @return rho_Z^l in the SM
     */
    complex rhoZ_l_SM(const StandardModel::lepton l) const;

    /**
     * @brief SM contribution to effective coupling rho_Z^q
     * @param[in] q name of a quark
     * @return rho_Z^q in the SM
     */
    complex rhoZ_q_SM(const StandardModel::quark q) const;
    
    /**
     * @brief SM contribution to effective coupling kappa_Z^l
     * @param[in] l name of a lepton 
     * @return kappa_Z^l in the SM
     */
    complex kappaZ_l_SM(const StandardModel::lepton l) const;
    
    /**
     * @brief SM contribution to effective coupling kappa_Z^q
     * @param[in] q name of a quark
     * @return kappa_Z^q in the SM
     */
    complex kappaZ_q_SM(const StandardModel::quark q) const;
    
    /**
     * @return rho_ij^W for Gamma(W->l nu) in the SM
     * @attention Fermion masses are neglected. 
     */
    double rho_GammaW_l_SM() const;
    
    /**
     * @return rho_ij^W for Gamma(W->q qbar) in the SM
     * @attention Fermion masses are neglected. 
     */
    double rho_GammaW_q_SM() const;
    
    /**
     * @param[in] li name of neutrinos
     * @param[in] lj name of charged leptons
     * @return the partial width of W^+ decay into an l_i\bar{l_j} pair in the SM
     * @attention Mixings in the lepton sector are neglected. 
     * @attention Fermion masses are neglected. 
     */
    double GammaW_l_SM(const StandardModel::lepton li, 
                       const StandardModel::lepton lj) const;
        
    /**
     * @param[in] qi name of up-type quark
     * @param[in] qj name of down-type quark
     * @return the partial width of W^+ decay into a q_i\bar{q_j} pair in the SM 
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
     * @brief effective coupling rho_Z^l
     * @param[in] l name of a lepton 
     * @param[in] S oblique parameter
     * @param[in] T oblique parameter
     * @return rho_Z^l in a NP model from oblique parameters
     */
    complex rhoZ_l_NP_fromSTU(const StandardModel::lepton l, const double T);

    /**
     * @brief effective coupling rho_Z^q
     * @param[in] q name of a quark
     * @param[in] S oblique parameter
     * @param[in] T oblique parameter
     * @return rho_Z^q in a NP model from oblique parameters
     */
    complex rhoZ_q_NP_fromSTU(const StandardModel::quark q, const double T);
    
    /**
     * @brief the ratio of the effective couplings for neutral-current interactions
     * @param[in] l name of a lepton 
     * @param[in] S oblique parameter
     * @param[in] T oblique parameter
     * @return g_V^l/g_A^l in a NP model from oblique parameters
     */
    complex gZl_over_gAl_NP_fromSTU(const StandardModel::lepton l,
                                    const double S, const double T);
    
    /**
     * @brief the ratio of the effective couplings for neutral-current interactions
     * @param[in] q name of a quark
     * @param[in] S oblique parameter
     * @param[in] T oblique parameter
     * @return g_V^q/g_A^q in a NP model from oblique parameters
     */
    complex gZq_over_gAq_NP_fromSTU(const StandardModel::quark q, 
                                    const double S, const double T);
    
    /**
     * @param[in] S oblique parameter
     * @param[in] T oblique parameter
     * @param[in] U oblique parameter
     * @return the total width of the W boson in a NP model from oblique parameters
     */    
    double GammaW_NP_fromSTU(const double S, const double T, const double U);    
    
    
    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return a pointer to the EWSMcommon object in EWSM class
     */
    EWSMcommon* getEWSMC() const {
        return EWSMC;
    }

    
    ////////////////////////////////////////////////////////////////////////     
protected:

    const StandardModel& SM;
    
    
    ////////////////////////////////////////////////////////////////////////         
private:
    
    bool flag_order[orders_EW_size]; 
    schemes_EW schemeMw, schemeRhoZ, schemeKappaZ;
    
    EWSMcommon* EWSMC;
    OneLoopEW* myOneLoopEW;
    TwoLoopQCD* myTwoLoopQCD;
    ThreeLoopQCD* myThreeLoopQCD;
    TwoLoopEW* myTwoLoopEW;
    ThreeLoopEW2QCD* myThreeLoopEW2QCD;
    ThreeLoopEW* myThreeLoopEW; 
    ApproximateFormulae* myApproximateFormulae;
        
    // accuracy in the iterative calculation of Mw
    static const double Mw_error = 0.00001; /* 0.01 MeV */ 
    


    
    
    
    
    // Cache
    static const int CacheSize = 5;
    mutable double Mw_cache[CacheSize];
    void CacheShift(double cache[][CacheSize], int n) const {
        int i,j;
        for(i=CacheSize-1;i>0;i--)
            for(j=0;j<n;j++)
                cache[j][i] = cache[j][i-1];
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

