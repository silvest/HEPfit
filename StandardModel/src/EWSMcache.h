/* 
 * File:   EWSMcache.h
 * Author: mishima
 */

#ifndef EWSMCACHE_H
#define	EWSMCACHE_H

#include <cmath>
#include <PVfunctions.h>
#include <Polylogarithms.h>
#include <ClausenFunctions.h>
#include "StandardModel.h"
using namespace gslpp;


class EWSMcache {

public:

    /**
     * @brief EWSMcache constructor
     * @param[in] SM_i reference to a StandardModel object
     */
    EWSMcache(const StandardModel& SM_i);

    static const int CacheSize = 5;

    int CacheCheck(const double cache[][CacheSize], 
                   const int NumPar, const double params[]) const;
    int CacheCheck(const complex cache[][CacheSize], 
                   const int NumPar, const double params[]) const;
    void CacheShift(double cache[][CacheSize], const int NumPar, 
                    const double params[], const double newResult) const;
    void CacheShift(complex cache[][CacheSize], const int NumPar, 
                    const double params[], const complex newResult) const;    
    
    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return a reference to the StandardModel object in the current class
     */
    //const StandardModel& getSM() const {
    //    return SM;
    //}

    /**
     * @return an object of ClausenFunctions class
     */
    const ClausenFunctions GetClausen() const {
        return Clausen;
    }

    
    ////////////////////////////////////////////////////////////////////////         
    
    /**
     * @return the zeta function zeta(2)
     */
    double GetZeta2() const {
        return zeta2;
    }

    /**
     * @return the zeta function zeta(3)
     */
    double GetZeta3() const {
        return zeta3;
    }
    
    /**
     * @return the zeta function zeta(4)
     */
    double GetZeta4() const {
        return zeta4;
    }        
    
    /**
     * @return the zeta function zeta(5)
     */
    double GetZeta5() const {
        return zeta5;
    }
    
    double GetS2() const {
        return S2;
    }
 
    double GetD3() const {
        return D3;
    }   
    
    double GetB4() const {
        return B4;
    }    

    /**
     * @return the logarithmic function log(2)
     */
    double GetLog2() const {
        return log2;
    }    

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @param[in] l name of lepton
     * @return mass of lepton
     */
    double ml(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of quark
     * @return mass of quark
     */
    double mq(const StandardModel::quark q) const;    
    
    /**
     * @return the top-quark mass
     */
    double Mt() const;
    
    /**
     * @return alpha_s(M_z^2)
     */
    double alsMz() const;
    
    /**
     * @return G_F
     */
    double GF() const;

    /**
     * @return electromagnetic coupling at q^2=0
     */
    double ale() const;

    /**
     * @return five-flavour hadronic correction to alpha at Mz^2
     */
    double dAle5Mz() const;
    
    /**
     * @return the Z-boson mass
     */
    double Mz() const;

    /**
     * @return the higgs mass
     */
    double mh() const;
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return the W-boson mass
     */
    double Mw(const double Mw_i) const;
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return c_W^2
     */
    double cW2(const double Mw_i) const;    
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return s_W^2
     */
    double sW2(const double Mw_i) const;    
        
    
    /**
     * @param[in] l name of lepton
     * @return electric charge of a lepton "l"
     */
    double Ql(const StandardModel::lepton l) const;    

    /**
     * @param[in] name of quark
     * @return electric charge of a quark "q"
     */
    double Qq(const StandardModel::quark q) const;  
    
    /**
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return the tree-level vector coupling for Z->l lbar
     * @attention depends on sW2
     */
    double vl(const StandardModel::lepton l, const double Mw_i) const;

    /**
     * @param[in] f name of quark
     * @param[in] Mw_i the W-boson mass
     * @return the tree-level vector coupling for Z->q qbar
     * @attention depends on sW2
     */
    double vq(const StandardModel::quark q, const double Mw_i) const;
    
    /**
     * @param[in] l name of lepton
     * @return the tree-level axial-vector coupling for Z->l lbar
     */
    double al(const StandardModel::lepton l) const;

    /**
     * @param[in] name of quark
     * @return the tree-level axial-vector coupling for Z->q qbar
     */
    double aq(const StandardModel::quark q) const;    
    
    /**
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return |v_l+a_l| 
     * @attention depends on sW2
     */
    double sigmal(const StandardModel::lepton l, const double Mw_i) const;

    /**
     * @param[in] name of quark
     * @param[in] Mw_i the W-boson mass
     * @return |v_q+a_q| 
     * @attention depends on sW2
     */
    double sigmaq(const StandardModel::quark q, const double Mw_i) const;    
    
    /**
     * @param[in] l name of lepton
     * @param[in] Mw_i the W-boson mass
     * @return v_l-a_l 
     * @attention depends on sW2
     */    
    double deltal(const StandardModel::lepton l, const double Mw_i) const;    

    /**
     * @param[in] name of quark
     * @param[in] Mw_i the W-boson mass
     * @return v_q-a_q 
     * @attention depends on sW2
     */    
    double deltaq(const StandardModel::quark q, const double Mw_i) const;   
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return the conversion factor from alpha to GF
     */
    double f_AlphaToGF(const double Mw_i) const;
   
    /**
     * @return X_t with G_F
     */
    double Xt_GF() const;
    
    /**
     * @param[in] Mw_i the W-boson mass
     * @return X_t with alpha(0)
     */
    double Xt_alpha(const double Mw_i) const;
    
    /**
     * @return alpha_s(M_t^2)
     */
    double alsMt() const;

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Phi_QCD2() const;

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double gamma_QCD2() const;

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double h_QCD2() const;
    
    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double logV1primeAndA1prime() const;
    
    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Cl2_2Phi() const;

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Cl2_4Phi() const;

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Cl3_2Phi() const;

    /**
     * @brief for two-loop QCD corrections
     * @return 
     */
    double Cl3_4Phi() const;

    
    ////////////////////////////////////////////////////////////////////////  
    
    double logMZtoME() const;
    double logMZtoMMU() const;
    double logMZtoMTAU() const;
    double logMZtoMTOP() const;
    double logMTOPtoMH() const;
    double log_cW2(const double Mw_i) const;
    
    double Li2_MW2toMTOP2(const double Mw_i) const;
    double Li3_MW2toMTOP2(const double Mw_i) const;
    double Li3_for_F1(const double Mw_i) const;
    
    double A0_Mz_Mw(const double Mw_i) const;
    double A0_Mz_mh() const;
    double A0_Mw_Mz(const double Mw_i) const;
    double A0_Mw_mh(const double Mw_i) const;
    double A0_Mz_Mz() const;    
    double A0_Mw_Mw(const double Mw_i) const;    
    
    complex B0_Mz_Mw2_mh_Mw(const double Mw_i) const;    
    complex B0_Mz_0_mh_Mw(const double Mw_i) const;    
    complex B0_Mw_Mz2_Mt_Mt(const double Mw_i) const;    
    complex B0_Mz_Mz2_Mw_Mw(const double Mw_i) const;    
    complex B0_Mz_Mz2_mh_Mz() const;
    complex B0_Mz_Mw2_Mz_Mw(const double Mw_i) const;    
    complex B0_Mz_Mw2_0_Mw(const double Mw_i) const;    
    complex B0_Mz_0_Mz_Mw(const double Mw_i) const;    
    complex B0_Mz_0_0_Mw(const double Mw_i) const;    
    complex B0_Mw_Mz2_Mw_Mw(const double Mw_i) const;    
    complex B0_Mw_Mw2_Mz_Mw(const double Mw_i) const;    
    complex B0_Mw_Mw2_mh_Mw(const double Mw_i) const;    
    complex B0_Mw_Mw2_0_Mw(const double Mw_i) const;    
    complex B0_Mz_Mz2_ml_ml(const StandardModel::lepton l) const;
    complex B0_Mz_Mz2_mq_mq(const StandardModel::quark q) const;
    
    complex B0p_Mz_0_mh_Mw(const double Mw_i) const;
    complex B0p_Mz_Mz2_mh_Mz() const;
    complex B0p_Mz_0_Mz_Mw(const double Mw_i) const;
    complex B0p_Mz_Mz2_Mw_Mw(const double Mw_i) const;
    complex B0p_Mw_Mw2_Mz_Mw(const double Mw_i) const;
    complex B0p_Mw_Mw2_mh_Mw(const double Mw_i) const;
    complex B0p_Mw_Mw2_0_Mw(const double Mw_i) const;
    complex B0p_Mz_Mz2_ml_ml(const StandardModel::lepton l) const;
    complex B0p_Mz_Mz2_mq_mq(const StandardModel::quark q) const;
    
    complex B1_Mz_0_ml_mlprime(const int gen) const;
    complex B1_Mz_0_mq_mqprime(const int gen) const;
    complex B1_Mz_0_mlprime_ml(const int gen) const;
    complex B1_Mz_0_mqprime_mq(const int gen) const;
    complex B1_Mz_Mw2_ml_mlprime(const int gen, const double Mw_i) const;
    complex B1_Mz_Mw2_mq_mqprime(const int gen, const double Mw_i) const;
    complex B1_Mz_Mw2_mlprime_ml(const int gen, const double Mw_i) const;
    complex B1_Mz_Mw2_mqprime_mq(const int gen, const double Mw_i) const;
    
    complex B1p_Mw_Mw2_ml_mlprime(const int gen, const double Mw_i) const;
    complex B1p_Mw_Mw2_mq_mqprime(const int gen, const double Mw_i) const;
    complex B1p_Mw_Mw2_mlprime_ml(const int gen, const double Mw_i) const;
    complex B1p_Mw_Mw2_mqprime_mq(const int gen, const double Mw_i) const;
    
    complex Bf_Mz_Mz2_ml_ml(const StandardModel::lepton l) const;
    complex Bf_Mz_Mz2_mq_mq(const StandardModel::quark q) const;
    complex Bf_Mz_0_ml_ml(const StandardModel::lepton l) const;
    complex Bf_Mz_0_mq_mq(const StandardModel::quark q) const;
    complex Bf_Mz_Mw2_mlprime_ml(const int gen, const double Mw_i) const;
    complex Bf_Mz_Mw2_mqprime_mq(const int gen, const double Mw_i) const;
    complex Bf_Mz_0_mlprime_ml(const int gen) const;
    complex Bf_Mz_0_mqprime_mq(const int gen) const;
    complex Bf_Mw_Mw2_mlprime_ml(const int gen, const double Mw_i) const;
    complex Bf_Mw_Mw2_mqprime_mq(const int gen, const double Mw_i) const;
    
    complex Bfp_Mz_Mz2_ml_ml(const StandardModel::lepton l) const;
    complex Bfp_Mz_Mz2_mq_mq(const StandardModel::quark q) const;
    complex Bfp_Mw_Mw2_mlprime_ml(const int gen, const double Mw_i) const;
    complex Bfp_Mw_Mw2_mqprime_mq(const int gen, const double Mw_i) const;
    
    complex C0_Mz2_Mw_Mt_Mw(const double Mw_i) const;
    complex C0_Mz2_Mt_Mw_Mt(const double Mw_i) const;
    complex C0_Mz2_0_Mw_0(const double Mw_i) const;
    complex C0_Mz2_Mw_0_Mw(const double Mw_i) const;
    complex C0_Mw2_Mw_0_Mz(const double Mw_i) const;
    complex C0_Mw2_0_Mz_0(const double Mw_i) const;
    complex C0_Mz2_0_Mz_0() const;
    
    
    //////////////////////////////////////////////////////////////////////// 

private:

    const StandardModel& SM;
    const PVfunctions PV;
    const ClausenFunctions Clausen;
    const Polylogarithms PolyLog;
    
    // Constants 
    double zeta2, zeta3, zeta4, zeta5;
    double S2, D3, B4;    
    double log2;

    
    ////////////////////////////////////////////////////////////////////////     
    // Caches 
    
    /* Logarithms */
    mutable double logMZtoME_cache[3][CacheSize];
    mutable double logMZtoMMU_cache[3][CacheSize];
    mutable double logMZtoMTAU_cache[3][CacheSize];
    mutable double logMZtoMTOP_cache[3][CacheSize]; 
    mutable double logMTOPtoMH_cache[3][CacheSize]; 
    mutable double log_cW2_cache[3][CacheSize]; 
    
    /* Dilogarithm and Trilogarithm for two-loop QCD corrections */
    mutable double Li2_MW2toMTOP2_cache[3][CacheSize]; 
    mutable double Li3_MW2toMTOP2_cache[3][CacheSize]; 
    mutable double Li3_for_F1_cache[3][CacheSize]; 
    
    /* One-loop functions */
    mutable double A0_Mz_Mw_cache[3][CacheSize]; 
    mutable double A0_Mz_mh_cache[3][CacheSize]; 
    mutable double A0_Mw_Mz_cache[3][CacheSize]; 
    mutable double A0_Mw_mh_cache[3][CacheSize]; 
    mutable double A0_Mz_Mz_cache[2][CacheSize]; 
    mutable double A0_Mw_Mw_cache[2][CacheSize]; 

    mutable complex B0_Mz_Mw2_mh_Mw_cache[4][CacheSize];
    mutable complex B0_Mz_0_mh_Mw_cache[4][CacheSize];
    mutable complex B0_Mw_Mz2_Mt_Mt_cache[4][CacheSize];
    mutable complex B0_Mz_Mz2_Mw_Mw_cache[3][CacheSize];
    mutable complex B0_Mz_Mz2_mh_Mz_cache[3][CacheSize];    
    mutable complex B0_Mz_Mw2_Mz_Mw_cache[3][CacheSize];
    mutable complex B0_Mz_Mw2_0_Mw_cache[3][CacheSize];
    mutable complex B0_Mz_0_Mz_Mw_cache[3][CacheSize];
    mutable complex B0_Mz_0_0_Mw_cache[3][CacheSize];  
    mutable complex B0_Mw_Mz2_Mw_Mw_cache[3][CacheSize];
    mutable complex B0_Mw_Mw2_Mz_Mw_cache[3][CacheSize];
    mutable complex B0_Mw_Mw2_mh_Mw_cache[3][CacheSize];
    mutable complex B0_Mw_Mw2_0_Mw_cache[2][CacheSize];    
    mutable complex B0_Mz_Mz2_ml_ml_cache[6][3][CacheSize];
    mutable complex B0_Mz_Mz2_mq_mq_cache[6][3][CacheSize];
    
    mutable complex B0p_Mz_0_mh_Mw_cache[4][CacheSize]; 
    mutable complex B0p_Mz_Mz2_mh_Mz_cache[3][CacheSize]; 
    mutable complex B0p_Mz_0_Mz_Mw_cache[3][CacheSize]; 
    mutable complex B0p_Mz_Mz2_Mw_Mw_cache[3][CacheSize]; 
    mutable complex B0p_Mw_Mw2_Mz_Mw_cache[3][CacheSize]; 
    mutable complex B0p_Mw_Mw2_mh_Mw_cache[3][CacheSize]; 
    mutable complex B0p_Mw_Mw2_0_Mw_cache[2][CacheSize]; 
    mutable complex B0p_Mz_Mz2_ml_ml_cache[6][3][CacheSize];
    mutable complex B0p_Mz_Mz2_mq_mq_cache[6][3][CacheSize];
    
    mutable complex B1_Mz_0_ml_mlprime_cache[3][4][CacheSize];
    mutable complex B1_Mz_0_mq_mqprime_cache[3][4][CacheSize];
    mutable complex B1_Mz_0_mlprime_ml_cache[3][4][CacheSize];
    mutable complex B1_Mz_0_mqprime_mq_cache[3][4][CacheSize];
    mutable complex B1_Mz_Mw2_ml_mlprime_cache[3][5][CacheSize];
    mutable complex B1_Mz_Mw2_mq_mqprime_cache[3][5][CacheSize];
    mutable complex B1_Mz_Mw2_mlprime_ml_cache[3][5][CacheSize];
    mutable complex B1_Mz_Mw2_mqprime_mq_cache[3][5][CacheSize];
    
    mutable complex B1p_Mw_Mw2_ml_mlprime_cache[3][4][CacheSize];
    mutable complex B1p_Mw_Mw2_mq_mqprime_cache[3][4][CacheSize];
    mutable complex B1p_Mw_Mw2_mlprime_ml_cache[3][4][CacheSize];
    mutable complex B1p_Mw_Mw2_mqprime_mq_cache[3][4][CacheSize];
        
    mutable complex Bf_Mz_Mz2_ml_ml_cache[6][3][CacheSize];
    mutable complex Bf_Mz_Mz2_mq_mq_cache[6][3][CacheSize];
    mutable complex Bf_Mz_0_ml_ml_cache[6][3][CacheSize];
    mutable complex Bf_Mz_0_mq_mq_cache[6][3][CacheSize];
    mutable complex Bf_Mz_Mw2_mlprime_ml_cache[3][5][CacheSize];
    mutable complex Bf_Mz_Mw2_mqprime_mq_cache[3][5][CacheSize];
    mutable complex Bf_Mz_0_mlprime_ml_cache[3][4][CacheSize];
    mutable complex Bf_Mz_0_mqprime_mq_cache[3][4][CacheSize];
    mutable complex Bf_Mw_Mw2_mlprime_ml_cache[3][4][CacheSize];
    mutable complex Bf_Mw_Mw2_mqprime_mq_cache[3][4][CacheSize];

    mutable complex Bfp_Mz_Mz2_ml_ml_cache[6][3][CacheSize];
    mutable complex Bfp_Mz_Mz2_mq_mq_cache[6][3][CacheSize];    
    mutable complex Bfp_Mw_Mw2_mlprime_ml_cache[3][4][CacheSize];
    mutable complex Bfp_Mw_Mw2_mqprime_mq_cache[3][4][CacheSize];
    
    mutable complex C0_Mz2_Mw_Mt_Mw_cache[4][CacheSize]; 
    mutable complex C0_Mz2_Mt_Mw_Mt_cache[4][CacheSize]; 
    mutable complex C0_Mz2_0_Mw_0_cache[3][CacheSize]; 
    mutable complex C0_Mz2_Mw_0_Mw_cache[3][CacheSize]; 
    mutable complex C0_Mw2_Mw_0_Mz_cache[3][CacheSize]; 
    mutable complex C0_Mw2_0_Mz_0_cache[3][CacheSize]; 
    mutable complex C0_Mz2_0_Mz_0_cache[2][CacheSize]; 
    
};

#endif	/* EWSMCACHE_H */

