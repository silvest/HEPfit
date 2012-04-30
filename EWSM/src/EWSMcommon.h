/* 
 * File:   EWSMcommon.h
 * Author: mishima
 */

#ifndef EWSMCOMMON_H
#define	EWSMCOMMON_H

#include <StandardModel.h>
#include <PVfunctions.h>
#include <Polylogarithms.h>
#include <ClausenFunctions.h>


class EWSMcommon {

public:

    /**
     * @brief EWSMcommon constructor
     * @param[in] SM_i reference to a StandardModel object
     */
    EWSMcommon(const StandardModel& SM_i);

    /**
     * @brief EWSMcommon copy constructor
     * @param[in] orig reference to an EWSMcommon object
     */
    //EWSMcommon(const EWSMcommon& orig);

    /**
     * @brief EWSMcommon destructor
     */
    virtual ~EWSMcommon();

    
    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @brief computes common variables which are independent of the W boson mass
     * @attention called from the constructor of the current class
     */
    void SetConstants();

    /**
     * @param[in] Mw_i the W boson mass
     * @brief computes Mw-dependent variables for computations of Mw
     */
    void ComputeForCC(const double Mw_i);

    /**
     * @param[in] Mw_i the W boson mass
     * @brief computes Mw-dependent variables for computations of rho_Z^f and kappa_Z^f
     */    
    void ComputeForNC(const double Mw_i);    
    
    /**
     * @param[in] Mw_i the W boson mass
     * @brief computes Mw-dependent variables for computations of rho^W_{ij}
     */       
    void ComputeForRhoWij(const double Mw_i);
    
    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @param[in] l name of a lepton
     * @return electric charge of lepton "l"
     */
    double Qf(const StandardModel::lepton l) const;
    
    /**
     * @param[in] q name of a quark
     * @return electric charge of quark "q"
     */
    double Qf(const StandardModel::quark q) const;    
    
    /**
     * @param[in] l name of a lepton
     * @return the tree-level vector coupling for Z->l+lbar
     * @attention depends on sW2
     */
    double vf(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the tree-level vector coupling for Z->q+qbar
     * @attention depends on sW2
     */    
    double vf(const StandardModel::quark q) const;    
    
    /**
     * @param[in] l name of a lepton
     * @return the tree-level axial-vector coupling for Z->l+lbar
     */
    double af(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return the tree-level axial-vector coupling for Z->q+qbar
     */    
    double af(const StandardModel::quark q) const;  

    /**
     * @param[in] l name of a lepton
     * @return |v_f+a_f| for f=l
     * @attention depends on sW2
     */
    double sigmaf(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return |v_f+a_f| for f=q
     * @attention depends on sW2
     */
    double sigmaf(const StandardModel::quark q) const;    
    
    /**
     * @param[in] l name of a lepton
     * @return v_f-a_f for f=l
     * @attention depends on sW2
     */    
    double deltaf(const StandardModel::lepton l) const;    
    
    /**
     * @param[in] q name of a quark
     * @return v_f-a_f for f=q
     * @attention depends on sW2
     */    
    double deltaf(const StandardModel::quark q) const; 
    
    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @return a reference to the StandardModel object in EWSMcommon class
     */
    const StandardModel& GetSM() const {
        return SM;
    }
    
    /**
     * @return the W boson mass
     */
    double GetMw() const {
        return Mw;
    }

    /**
     * @return c_W^2
     */
    double GetCW2() const {
        return cW2;
    }

    /**
     * @return s_W^2
     */
    double GetSW2() const {
        return sW2;
    }
    
    /**
     * @return the conversion factor from alpha to GF
     */
    double GetF_AlphaToGF() const {
        return f_AlphaToGF;
    }
   
    double GetXt_GF() const {
        return Xt_GF;
    }    
    
    double GetXt_alpha() const {
        return Xt_alpha;
    }  
    
    double GetAlsMt() const {
        return AlsMt;
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

    /**
     * @return the logarithmic function log(2)
     */
    double GetLog2() const {
        return log2;
    }

    /**
     * @return the log of Mz/me
     */
    double GetLogMZtoME() const {
        return logMZtoME;
    }

    /**
     * @return the log of Mz/mmu
     */
    double GetLogMZtoMMU() const {
        return logMZtoMMU;
    }

    /**
     * @return the log of Mz/mtau
     */
    double GetLogMZtoMTAU() const {
        return logMZtoMTAU;
    }
    
    /**
     * @return the log of Mz/mtop
     */
    double GetLogMZtoMTOP() const {
        return logMZtoMTOP;
    }

    double GetLogMTOPtoMH() const {
        return logMTOPtoMH;
    }

    double GetLog_cW2() const {
        return log_cW2;
    }
    
    double GetLi2_MW2toMTOP2() const {
        return Li2_MW2toMTOP2;
    }

    double GetLi3_MW2toMTOP2() const {
        return Li3_MW2toMTOP2;
    }

    double GetLi3_for_F1() const {
        return Li3_for_F1;
    }
    
    double GetPhi_QCD2() const {
        return Phi_QCD2;
    }

    double GetGamma_QCD2() const {
        return gamma_QCD2;
    }

    double GetH_QCD2() const {
        return h_QCD2;
    }
    
    double GetLogV1primeAndA1prime() const {
        return logV1primeAndA1prime;
    }
    
    double GetCl2_2Phi() const {
        return Cl2_2Phi;
    }

    double GetCl2_4Phi() const {
        return Cl2_4Phi;
    }

    double GetCl3_2Phi() const {
        return Cl3_2Phi;
    }

    double GetCl3_4Phi() const {
        return Cl3_4Phi;
    }
    
    double GetA0_Mz_Mw() const {
        return A0_Mz_Mw;
    }

    double GetA0_Mz_Mz() const {
        return A0_Mz_Mz;
    }

    double GetA0_Mz_mh() const {
        return A0_Mz_mh;
    }

    complex GetB0_Mz_Mw2_0_Mw() const {
        return B0_Mz_Mw2_0_Mw;
    }

    complex GetB0_Mz_Mw2_Mz_Mw() const {
        return B0_Mz_Mw2_Mz_Mw;
    }

    complex GetB0_Mz_Mw2_mh_Mw() const {
        return B0_Mz_Mw2_mh_Mw;
    }
    complex GetB0_Mz_0_0_Mw() const {
        return B0_Mz_0_0_Mw;
    }

    complex GetB0_Mz_0_Mz_Mw() const {
        return B0_Mz_0_Mz_Mw;
    }

    complex GetB0_Mz_0_mh_Mw() const {
        return B0_Mz_0_mh_Mw;
    }
    
    complex GetB0_Mz_Mz2_Mw_Mw() const {
        return B0_Mz_Mz2_Mw_Mw;
    }
    
    complex GetB0_Mz_Mz2_mh_Mz() const {
        return B0_Mz_Mz2_mh_Mz;
    }

    complex GetB0_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
        return B0_Mz_Mz2_ml_ml[l];
    }

    complex GetB0_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
        return B0_Mz_Mz2_mq_mq[q];
    }

    complex GetBf_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
        return Bf_Mz_Mz2_ml_ml[l];
    }

    complex GetBf_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
        return Bf_Mz_Mz2_mq_mq[q];
    }

    complex GetBf_Mz_0_ml_ml(const StandardModel::lepton l) const {
        return Bf_Mz_0_ml_ml[l];
    }

    complex GetBf_Mz_0_mq_mq(const StandardModel::quark q) const {
        return Bf_Mz_0_mq_mq[q];
    }    
    
    complex GetB1_Mz_Mw2_ml_mlprime(const int gen) const {
        return B1_Mz_Mw2_ml_mlprime[gen];
    }

    complex GetB1_Mz_Mw2_mq_mqprime(const int gen) const {
        return B1_Mz_Mw2_mq_mqprime[gen];
    }
    
    complex GetB1_Mz_Mw2_mlprime_ml(const int gen) const {
        return B1_Mz_Mw2_mlprime_ml[gen];
    }

    complex GetB1_Mz_Mw2_mqprime_mq(const int gen) const {
        return B1_Mz_Mw2_mqprime_mq[gen];
    }

    complex GetBf_Mz_Mw2_mlprime_ml(const int gen) const {
        return Bf_Mz_Mw2_mlprime_ml[gen];
    }

    complex GetBf_Mz_Mw2_mqprime_mq(const int gen) const {
        return Bf_Mz_Mw2_mqprime_mq[gen];
    }
    
    complex GetB1_Mz_0_ml_mlprime(const int gen) const {
        return B1_Mz_0_ml_mlprime[gen];
    }

    complex GetB1_Mz_0_mq_mqprime(const int gen) const {
        return B1_Mz_0_mq_mqprime[gen];
    }
    
    complex GetB1_Mz_0_mlprime_ml(const int gen) const {
        return B1_Mz_0_mlprime_ml[gen];
    }

    complex GetB1_Mz_0_mqprime_mq(const int gen) const {
        return B1_Mz_0_mqprime_mq[gen];
    }

    complex GetBf_Mz_0_mlprime_ml(const int gen) const {
        return Bf_Mz_0_mlprime_ml[gen];
    }

    complex GetBf_Mz_0_mqprime_mq(const int gen) const {
        return Bf_Mz_0_mqprime_mq[gen];
    }
        
    complex GetB0p_Mz_0_Mz_Mw() const {
        return B0p_Mz_0_Mz_Mw;
    }

    complex GetB0p_Mz_0_mh_Mw() const {
        return B0p_Mz_0_mh_Mw;
    }
    
    complex GetB0p_Mz_Mz2_Mw_Mw() const {
        return B0p_Mz_Mz2_Mw_Mw;
    }

    complex GetB0p_Mz_Mz2_mh_Mz() const {
        return B0p_Mz_Mz2_mh_Mz;
    }
    
    complex GetB0p_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
        return B0p_Mz_Mz2_ml_ml[l];
    }

    complex GetB0p_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
        return B0p_Mz_Mz2_mq_mq[q];
    }

    complex GetBfp_Mz_Mz2_ml_ml(const StandardModel::lepton l) const {
        return Bfp_Mz_Mz2_ml_ml[l];
    }

    complex GetBfp_Mz_Mz2_mq_mq(const StandardModel::quark q) const {
        return Bfp_Mz_Mz2_mq_mq[q];
    }

    complex GetB0_Mw_Mz2_Mt_Mt() const {
        return B0_Mw_Mz2_Mt_Mt;
    }

    complex GetB0_Mw_Mz2_Mw_Mw() const {
        return B0_Mw_Mz2_Mw_Mw;
    }

    complex GetC0_Mz2_0_Mw_0() const {
        return C0_Mz2_0_Mw_0;
    }

    complex GetC0_Mz2_0_Mz_0() const {
        return C0_Mz2_0_Mz_0;
    }

    complex GetC0_Mz2_Mt_Mw_Mt() const {
        return C0_Mz2_Mt_Mw_Mt;
    }

    complex GetC0_Mz2_Mw_0_Mw() const {
        return C0_Mz2_Mw_0_Mw;
    }

    complex GetC0_Mz2_Mw_Mt_Mw() const {
        return C0_Mz2_Mw_Mt_Mw;
    }
    
    double GetA0_Mw_Mw() const {
        return A0_Mw_Mw;
    }

    double GetA0_Mw_Mz() const {
        return A0_Mw_Mz;
    }

    double GetA0_Mw_mh() const {
        return A0_Mw_mh;
    }

    complex GetB0_Mw_Mw2_0_Mw() const {
        return B0_Mw_Mw2_0_Mw;
    }

    complex GetB0_Mw_Mw2_Mz_Mw() const {
        return B0_Mw_Mw2_Mz_Mw;
    }

    complex GetB0_Mw_Mw2_mh_Mw() const {
        return B0_Mw_Mw2_mh_Mw;
    }

    complex GetB0p_Mw_Mw2_0_Mw() const {
        return B0p_Mw_Mw2_0_Mw;
    }

    complex GetB0p_Mw_Mw2_Mz_Mw() const {
        return B0p_Mw_Mw2_Mz_Mw;
    }

    complex GetB0p_Mw_Mw2_mh_Mw() const {
        return B0p_Mw_Mw2_mh_Mw;
    }

    complex GetB1p_Mw_Mw2_ml_mlprime(const int gen) const {
        return B1p_Mw_Mw2_ml_mlprime[gen];
    }

    complex GetB1p_Mw_Mw2_mlprime_ml(const int gen) const {
        return B1p_Mw_Mw2_mlprime_ml[gen];
    }

    complex GetB1p_Mw_Mw2_mq_mqprime(const int gen) const {
        return B1p_Mw_Mw2_mq_mqprime[gen];
    }

    complex GetB1p_Mw_Mw2_mqprime_mq(const int gen) const {
        return B1p_Mw_Mw2_mqprime_mq[gen];
    }

    complex GetBf_Mw_Mw2_mlprime_ml(const int gen) const {
        return Bf_Mw_Mw2_mlprime_ml[gen];
    }

    complex GetBf_Mw_Mw2_mqprime_mq(const int gen) const {
        return Bf_Mw_Mw2_mqprime_mq[gen];
    }

    complex GetBfp_Mw_Mw2_mlprime_ml(const int gen) const {
        return Bfp_Mw_Mw2_mlprime_ml[gen];
    }

    complex GetBfp_Mw_Mw2_mqprime_mq(const int gen) const {
        return Bfp_Mw_Mw2_mqprime_mq[gen];
    }
    
    complex GetC0_Mw2_0_Mz_0() const {
        return C0_Mw2_0_Mz_0;
    }

    complex GetC0_Mw2_Mw_0_Mz() const {
        return C0_Mw2_Mw_0_Mz;
    }

    
    //////////////////////////////////////////////////////////////////////// 

protected:
    const StandardModel& SM;
    PVfunctions PV;
    ClausenFunctions Clausen;
    Polylogarithms PolyLog;
    
    
    // Cache
    static const int CacheSize = 5;
    mutable double cache[CacheSize];
    
    
    
    double Mw;
    double cW2, sW2;    

    /* local variables, copies of the corresponding ones in StandardModel */
    double Mz, mh, Mt, ml[6], mq[6]; 

    double f_AlphaToGF;    
    double Xt_GF; /* X_t with G_F */
    double Xt_alpha; /* X_t with alpha(0) */
    double AlsMt; /* alpha_s(M_t) */
    
    double S2, D3, B4;
    
    double zeta2, zeta3, zeta4, zeta5;

    /* Logarithms */
    double log2;
    double logMZtoME, logMZtoMMU, logMZtoMTAU, logMZtoMTOP; 
    double logMTOPtoMH;
    double log_cW2;
    
    /* Dilogarithm and Trilogarithm */
    double Li2_MW2toMTOP2;
    double Li3_MW2toMTOP2, Li3_for_F1;
    
    /* Logarithms, Clausen functions, etc for two-loop QCD corrections */
    double Phi_QCD2, gamma_QCD2, h_QCD2;
    double logV1primeAndA1prime;
    double Cl3_2Phi, Cl3_4Phi, Cl2_2Phi, Cl2_4Phi;
    
    /* One-loop functions in self-energies */
    double A0_Mz_Mw;
    double A0_Mz_Mz;    
    double A0_Mz_mh;
    complex B0_Mz_Mw2_Mz_Mw;
    complex B0_Mz_Mw2_0_Mw;
    complex B0_Mz_Mw2_mh_Mw;
    complex B0_Mz_Mz2_Mw_Mw;
    complex B0_Mz_Mz2_mh_Mz;    
    complex B0_Mz_0_Mz_Mw;
    complex B0_Mz_0_0_Mw;
    complex B0_Mz_0_mh_Mw;
    complex B0_Mz_Mz2_ml_ml[6];
    complex B0_Mz_Mz2_mq_mq[6];
    complex B1_Mz_Mw2_ml_mlprime[3];
    complex B1_Mz_Mw2_mq_mqprime[3];
    complex B1_Mz_Mw2_mlprime_ml[3];
    complex B1_Mz_Mw2_mqprime_mq[3];
    complex B1_Mz_0_ml_mlprime[3];
    complex B1_Mz_0_mq_mqprime[3];
    complex B1_Mz_0_mlprime_ml[3];
    complex B1_Mz_0_mqprime_mq[3];
    complex Bf_Mz_Mw2_mlprime_ml[3];
    complex Bf_Mz_Mw2_mqprime_mq[3];
    complex Bf_Mz_0_mlprime_ml[3];
    complex Bf_Mz_0_mqprime_mq[3];    
    complex Bf_Mz_Mz2_ml_ml[6];
    complex Bf_Mz_Mz2_mq_mq[6]; 
    complex Bf_Mz_0_ml_ml[6];
    complex Bf_Mz_0_mq_mq[6];    
    complex B0p_Mz_0_Mz_Mw;
    complex B0p_Mz_0_mh_Mw; 
    complex B0p_Mz_Mz2_Mw_Mw;
    complex B0p_Mz_Mz2_mh_Mz; 
    complex B0p_Mz_Mz2_ml_ml[6];
    complex B0p_Mz_Mz2_mq_mq[6];    
    complex Bfp_Mz_Mz2_ml_ml[6];
    complex Bfp_Mz_Mz2_mq_mq[6];    
    
    /* One-loop functions in vertex corrections */    
    complex B0_Mw_Mz2_Mw_Mw;
    complex B0_Mw_Mz2_Mt_Mt;
    complex C0_Mz2_0_Mz_0; 
    complex C0_Mz2_0_Mw_0;
    complex C0_Mz2_Mt_Mw_Mt;
    complex C0_Mz2_Mw_0_Mw;    
    complex C0_Mz2_Mw_Mt_Mw;    
    
    /* One-loop functions in the WF renormalization of the W boson */
    double A0_Mw_Mw;
    double A0_Mw_Mz;
    double A0_Mw_mh;
    complex B0_Mw_Mw2_Mz_Mw;
    complex B0_Mw_Mw2_0_Mw;
    complex B0_Mw_Mw2_mh_Mw;
    complex B0p_Mw_Mw2_Mz_Mw;
    complex B0p_Mw_Mw2_0_Mw;
    complex B0p_Mw_Mw2_mh_Mw;
    complex Bf_Mw_Mw2_mlprime_ml[3];
    complex Bf_Mw_Mw2_mqprime_mq[3];              
    complex Bfp_Mw_Mw2_mlprime_ml[3];
    complex Bfp_Mw_Mw2_mqprime_mq[3];            
    complex B1p_Mw_Mw2_ml_mlprime[3];
    complex B1p_Mw_Mw2_mq_mqprime[3];
    complex B1p_Mw_Mw2_mlprime_ml[3];
    complex B1p_Mw_Mw2_mqprime_mq[3];
    complex C0_Mw2_0_Mz_0;
    complex C0_Mw2_Mw_0_Mz;
    
};

#endif	/* EWSMCOMMON_H */

