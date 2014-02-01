/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSM_H
#define	EWSM_H

#include <cstring>
#include <stdexcept>
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

class EWSM_Output; // forward reference to EWSM_Output class

/**
 * @class EWSM
 * @ingroup StandardModel
 * @brief A class for radiative corrections to the EW precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details
 *
 *
 * the on-shell scheme:
 * @cite Bardin:1980fe
 * @cite Bardin:1981sv
 *
 *
 */
class EWSM {
public:

    /**
     * @brief Friend class of EWSM class.
     */
    friend class EWSM_Output;

    // accuracy in the iterative calculation of Mw
    static const double Mw_error;
    
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
    enum schemes_EW {NORESUM=0, OMSI, INTERMEDIATE, OMSII, APPROXIMATEFORMULA,
                     FIXED, schemes_EW_size};

    std::string SchemeString(const schemes_EW scheme) const
    {
        if (scheme == NORESUM) return "NORESUM";
        else if (scheme == OMSI) return "OMSI";
        else if (scheme == INTERMEDIATE) return "INTERMEDIATE";
        else if (scheme == OMSII) return "OMSII";
        else if (scheme == APPROXIMATEFORMULA) return "APPROXIMATEFORMULA";
        else if (scheme == FIXED) return "FIXED";
        else
            throw std::runtime_error("Undefined scheme in EWSM::SchemeString");
    }

    // The number of the parameters relevant to EW observables
    static const int NumSMParams = 24;
        
    
    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @brief Constructor. 
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EWSM(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief
     * @return
     */
    EWSMcache* getMyCache() const
    {
        return myCache;
    }

    /**
     * @brief
     * @return
     */
    EWSMOneLoopEW* getMyOneLoopEW() const
    {
        return myOneLoopEW;
    }

    /**
     * @brief
     * @return
     */
    EWSMApproximateFormulae* getMyApproximateFormulae() const
    {
        return myApproximateFormulae;
    }

    /**
     * @return a pointer to the EWSMTwoFermionsLEP2 object
     */
    EWSMTwoFermionsLEP2* getMyTwoFermionsLEP2() const 
    {
        return myTwoFermionsLEP2;
    }

    schemes_EW getSchemeMw() const
    {
        return schemeMw;
    }

    /**
     * @brief set a flag for M_W
     * @param[in] schemeMw NORESUM, OMSI, INTERMEDIATE, OMSII or APPROXIMATEFORMULA
     */
    void setSchemeMw(schemes_EW schemeMw) 
    {
        this->schemeMw = schemeMw;
        std::cout << "Mw: " << SchemeString(schemeMw) << std::endl;
    }

    schemes_EW getSchemeRhoZ() const
    {
        return schemeRhoZ;
    }

    /**
     * @brief set a flag for rho_Z
     * @param[in] schemeRhoZ NORESUM, OMSI, INTERMEDIATE, OMSII or APPROXIMATEFORMULA
     * @attention NORESUM is not available, since reducible two-loop EW corrections have not been implemented yet. 
     */
    void setSchemeRhoZ(schemes_EW schemeRhoZ) 
    {
        this->schemeRhoZ = schemeRhoZ;
        std::cout << "rhoZf: " << SchemeString(schemeRhoZ) << std::endl;
    }

    schemes_EW getSchemeKappaZ() const
    {
        return schemeKappaZ;
    }
    
    /**
     * @brief set a flag for kappa_Z
     * @param[in] schemeKappaZ NORESUM, OMSI, INTERMEDIATE, OMSII or APPROXIMATEFORMULA
     */
    void setSchemeKappaZ(schemes_EW schemeKappaZ)
    {
        this->schemeKappaZ = schemeKappaZ;
        std::cout << "kappaZf: " << SchemeString(schemeKappaZ) << std::endl;
    }

    
    //////////////////////////////////////////////////////////////////////// 

    bool checkSMparams(double Params_cache[], const bool bUpdate=true) const;

    bool checkScheme(schemes_EW& scheme_cache, const schemes_EW scheme_current,
                     const bool bUpdate=true) const;

    
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
     * @return The W boson mass without weak corrections.
     */
    double Mw0() const;

    /**
     * @return @f$sin^2\theta_W@f$ without weak corrections.
     */
    double s02() const;

    /**
     * @return @f$\cos^2\theta_W@f$ without weak corrections.
     */
    double c02() const;


    ////////////////////////////////////////////////////////////////////////     
        
    /**
     * @return the W boson mass in the SM
     */
    double Mw_SM() const;
    
    /** 
     * @brief Delta r derived from Mw().
     * @return Delta r in the SM. 
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

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief The Z boson mass in the complex-pole scheme.
     * @return The Z boson mass in the complex-pole scheme.
     */
    double Mzbar() const;

    /**
     * @brief Converts the W-boson mass from the experimental scheme to the complex-pole one.
     * @param[in] Mw The W boson mass in the experimental scheme.
     * @return The W boson mass in the complex-pole scheme.
     */
    double MwbarFromMw(const double Mw) const;

    /**
     * @brief Converts the W-boson mass from the complex-pole scheme to the experimental one.
     * @param[in] Mwbar The W boson mass in the complex-pole scheme.
     * @return The W boson mass in the experimental scheme.
     */
    double MwFromMwbar(const double Mwbar) const;

    /**
     * @brief The SM prediction for Delta r derived from Mwbar in the complex-pole scheme.
     * @return The SM prediction for Delta r in the complex-pole scheme.
     */
    double DeltaRbar_SM() const;

    
    ////////////////////////////////////////////////////////////////////////
    // SM contribution to the effective couplings

    /**
     * @brief SM contribution to the effective coupling @f$\rho_Z^l@f$.
     * @param[in] l Name of a lepton.
     * @return @f$\rho_Z^l@f$ in the SM.
     */
    complex rhoZ_l_SM(const StandardModel::lepton l) const;

    /**
     * @brief SM contribution to the effective coupling @f$\rho_Z^q@f$.
     * @param[in] q Name of a quark.
     * @return @f$\rho_Z^q@f$ in the SM.
     */
    complex rhoZ_q_SM(const StandardModel::quark q) const;    
    
    /**
     * @brief SM contribution to the effective coupling @f$\kappa_Z^l@f$.
     * @param[in] l Name of a lepton.
     * @return @f$\kappa_Z^l@f$ in the SM.
     */
    complex kappaZ_l_SM(const StandardModel::lepton l) const;

    /**
     * @brief SM contribution to the effective coupling @f$\kappa_Z^q@f$.
     * @param[in] q Name of a quark.
     * @return @f$\kappa_Z^q@f$ in the SM.
     */
    complex kappaZ_q_SM(const StandardModel::quark q) const;    

    /**
     * @brief SM contribution to the effective coupling @f$g_V^l@f$.
     * @param[in] l Name of a lepton.
     * @return @f$g_V^l@f$ in the SM.
     */
    complex gVl_SM(const StandardModel::lepton l) const;

    /**
     * @brief SM contribution to the effective coupling @f$g_V^q@f$.
     * @param[in] q Name of a quark.
     * @return @f$g_V^q@f$ in the SM.
     */
    complex gVq_SM(const StandardModel::quark q) const;

    /**
     * @brief SM contribution to the effective coupling @f$g_A^l@f$.
     * @param[in] l Name of a lepton.
     * @return @f$g_A^l@f$ in the SM.
     */
    complex gAl_SM(const StandardModel::lepton l) const;

    /**
     * @brief SM contribution to the effective coupling @f$g_A^q@f$.
     * @param[in] q Name of a quark.
     * @return @f$g_A^q@f$ in the SM.
     */
    complex gAq_SM(const StandardModel::quark q) const;
    

    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] l Name of a lepton.
     * @return The effective coupling @f$\rho_Z^l@f$.
     */
    virtual complex rhoZ_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q Name of a quark.
     * @return The effective coupling neutral-current interactions @f$\rho_Z^q@f$.
     */
    virtual complex rhoZ_q(const StandardModel::quark q) const;

    /**
     * @param[in] l Name of a lepton.
     * @return The effective coupling neutral-current interactions @f$\kappa_Z^l@f$.
     */
    virtual complex kappaZ_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q Name of a quark.
     * @return The effective coupling neutral-current interactions @f$\kappa_Z^q@f$.
     */
    virtual complex kappaZ_q(const StandardModel::quark q) const;

    /**
     * @param[in] l Name of a lepton.
     * @return The effective vector coupling for neutral-current interactions @f$g_V^l@f$.
     */
    virtual complex gVl(const StandardModel::lepton l) const;

    /**
     * @param[in] q Name of a quark.
     * @return The effective vector coupling for neutral-current interactions @f$g_V^q@f$.
     */
    virtual complex gVq(const StandardModel::quark q) const;

    /**
     * @param[in] l Name of a lepton.
     * @return The effective axial-vector coupling for neutral-current interactions @f$g_A^l@f$.
     */
    virtual complex gAl(const StandardModel::lepton l) const;

    /**
     * @param[in] q Name of a quark.
     * @return The effective axial-vector coupling for neutral-current interactions @f$g_A^q@f$.
     */
    virtual complex gAq(const StandardModel::quark q) const;

    
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

    /**
     * @return SM contribution to @f$\epsilon_1@f$.
     */
    double epsilon1_SM() const;

    /**
     * @return SM contribution to @f$\epsilon_2@f$.
     */
    double epsilon2_SM() const;

    /**
     * @return SM contribution to @f$\epsilon_3@f$.
     */
    double epsilon3_SM() const;

    /**
     * @return SM contribution to @f$\epsilon_b@f$.
     */
    double epsilonb_SM() const;
    

    ////////////////////////////////////////////////////////////////////////     
    // The W-boson decay width

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
     * @return the partial width of W^+ decay into an @f$l_i\bar{l_j}@f$ pair in the SM
     * @attention Mixings in the lepton sector are neglected. 
     * @attention Fermion masses are neglected. 
     */
    double GammaW_l_SM(const StandardModel::lepton li, 
                       const StandardModel::lepton lj) const;

    /**
     * @param[in] qi name of a up-type quark
     * @param[in] qj name of a down-type quark
     * @return the partial width of W^+ decay into an @f$q_i\bar{q_j}@f$ pair in the SM
     * @attention Fermion masses are neglected. 
     */
    double GammaW_q_SM(const StandardModel::quark qi, 
                       const StandardModel::quark qj) const;    
    
    /**
     * @return the total width of the W boson in the SM
     */
    double GammaW_SM() const;   


    ////////////////////////////////////////////////////////////////////////
protected:
    const StandardModel& SM;

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

    mutable schemes_EW schemeMw_cache, schemeRhoZ_cache, schemeKappaZ_cache;


    ////////////////////////////////////////////////////////////////////////

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

