/* 
 * File:   OneLoopEW.h
 * Author: mishima
 */

#ifndef ONELOOPEW_H
#define	ONELOOPEW_H

#include <StandardModel.h>
#include "EWSMcommon.h"

using namespace gslpp;


class OneLoopEW {
    
public:

    /**
     * @brief OneLoopEW constructor
     * @param[in] EWSMC_i reference to an EWSMcommon object
     */
    OneLoopEW(const EWSMcommon& EWSMC_i);


    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief leptonic contribution to alpha
     * @return Delta alpha_{lept}^{alpha}
     */
    double DeltaAlpha_l() const;

    /**
     * @brief top-quark contribution to alpha
     * @return Delta alpha_{top}^{alpha}
     */
    double DeltaAlpha_t() const;

    /**
     * @brief leading contribution to Delta r
     * @return Delta rho^{alpha}
     */
    double DeltaRho() const;

    /**
     * @brief remainder contribution to Delta r
     * @return Delta r_{rem}^{alpha}
     */
    double DeltaR_rem() const;

    /**
     * @brief remainder contribution for rho_Z^f and kappa_Z^f
     * @return Delta rbar_{rem}^{alpha}
     */
    double DeltaRbar_rem() const;

    /**
     * @brief remainder contribution to rho_Z^f for a given u_f
     * @param[in] uf a combination of the unified form factors
     * @return delta rho_{rem}^{f, alpha}(u_f)
     */
    complex deltaRho_rem_tmp(const complex u_f) const;
    
    /**
     * @brief remainder contribution to rho_Z^l
     * @param[in] l name of a lepton 
     * @return delta rho_{rem}^{l, alpha}
     */
    complex deltaRho_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to rho_Z^q
     * @param[in] q name of a quark 
     * @return delta rho_{rem}^{q, alpha}
     */
    complex deltaRho_rem_q(const StandardModel::quark q) const;

    /**
     * @brief remainder contribution to kappa_Z^l for given delta_f and u_f
     * @param[in] deltaf a combination of the effective couplings 
     * @param[in] uf a combination of the unified form factors
     * @return delta kappa_{rem}^{f, alpha}(delta_f, u_f)
     */
    complex deltaKappa_rem_tmp(const double deltaf, const complex uf) const;    
    
    /**
     * @brief remainder contribution to kappa_Z^l
     * @param[in] l name of a lepton 
     * @return delta kappa_{rem}^{l, alpha}
     */
    complex deltaKappa_rem_l(const StandardModel::lepton l) const;

    /**
     * @brief remainder contribution to kappa_Z^q
     * @param[in] q name of a quark 
     * @return delta kappa_{rem}^{q, alpha}
     */
    complex deltaKappa_rem_q(const StandardModel::quark q) const;

    /**
     * @param[in] Qi the electric charge of f_i
     * @param[in] Qj the electric charge of f_j
     * @return O(alpha) contribution to the width of W -> f_i bar{f}_j for given Q_i and Q_j
     * @attention masses for virtual fermions are neglected. 
     */
    double rho_GammaW_tmp(const double Qi, const double Qj) const;    
    
    /**
     * @param[in] li name of lepton
     * @param[in] lj name of lepton
     * @return O(alpha) contribution to the width of W -> l_i bar{l}_j
     */
    double rho_GammaW_l(const StandardModel::lepton li, 
                        const StandardModel::lepton lj) const;

     /**
     * @param[in] qi name of quark
     * @param[in] qj name of quark
     * @return O(alpha) contribution to the width of W -> q_i bar{q}_j
     */
    double rho_GammaW_q(const StandardModel::quark qi, 
                        const StandardModel::quark qj) const;

    
    ////////////////////////////////////////////////////////////////////////    
       
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @return the bosonic contribution to the self-energy function of the W boson
     */
    complex SigmaWW_bos(const double mu, const double s) const;
 
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @return the fermionic contribution to the self-energy function of the W boson
     */
    complex SigmaWW_fer(const double mu, const double s) const;
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @return the bosonic contribution to the self-energy function of the Z boson
     */
    complex SigmaZZ_bos(const double mu, const double s) const;    
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @return the fermionic contribution to the self-energy function of the Z boson
     */
    complex SigmaZZ_fer(const double mu, const double s) const;
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @return the bosonic contribution to the self-energy function of the photon
     */
    complex PiGammaGamma_bos(const double mu, const double s) const;

     /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] l name of a lepton
     * @return the lepton "l" contribution to the self-energy function of the photon
     */
    complex PiGammaGamma_fer(const double mu, const double s, 
                             const StandardModel::lepton l) const;
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @param[in] q name of a quark 
     * @return the quark "q" contribution to the self-energy function of the photon
     */
    complex PiGammaGamma_fer(const double mu, const double s, 
                             const StandardModel::quark q) const;
    
    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @return the fermionic contribution to the self-energy function of the photon
     */
    complex PiGammaGamma_fer(const double mu, const double s) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @return the bosonic contribution to the self-energy function of the Z-gamma mixing
     */
    complex PiZgamma_bos(const double mu, const double s) const;

    /**
     * @param[in] mu renormalization scale
     * @param[in] s momentum-squared
     * @return the fermionic contribution to the self-energy function of the Z-gamma mixing
     */    
    complex PiZgamma_fer(const double mu, const double s) const;

    
    ////////////////////////////////////////////////////////////////////////   
    
    /**
     * @param[in] mu renormalization scale
     * @return the bosonic contribution to the wave-function renormalization of the W boson with s=Mw^2
     */   
    complex SigmaPrime_WW_bos_Mw2(const double mu) const;
    
     /**
     * @param[in] mu renormalization scale
     * @return the fermionic contribution to the wave-function renormalization of the W boson with s=Mw^2
     */   
    complex SigmaPrime_WW_fer_Mw2(const double mu) const;   
    
     /**
     * @param[in] mu renormalization scale
     * @return the bosonic contribution to the wave-function renormalization of the Z boson with s=Mz^2
     */   
    complex SigmaPrime_ZZ_bos_Mz2(const double mu) const;
    
    /**
     * @param[in] mu renormalization scale
     * @return the fermionic contribution to the wave-function renormalization of the Z boson with s=Mz^2
     */   
    complex SigmaPrime_ZZ_fer_Mz2(const double mu) const;     
    
    ////////////////////////////////////////////////////////////////////////       
    
    /**
     * @param[in] mu renormalization scale
     * @return Delta Rhobar^{F} for the renormalization scale mu
     */
    double DeltaRhobar(const double mu) const;
    
    /**
     * @param[in] mu renormalization scale
     * @return Delta Rhobar_W^{F} for the renormalization scale mu
     */
    double DeltaRhobarW(const double mu) const;
    
    
    ////////////////////////////////////////////////////////////////////////   
    
    /**
     * @brief TEST function
     * @return Delta Rhobar^{bos,F} for mu=Mw
     * @attention The renormalization scale is fixed to be mu=Mw.
     */
    double TEST_DeltaRhobar_bos() const;

    /**
     * @brief TEST function
     * @return Delta Rhobar_W^{bos,F}
     * @attention The renormalization scale is fixed to be mu=Mw.
     */
    double TEST_DeltaRhobarW_bos() const;    


    ////////////////////////////////////////////////////////////////////////    

    /**
     * @param[in] s momentum-squared
     * @return abelian-type vertex corrections with virtual Z in the chiral limit
     */
    complex FZa_0(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @return abelian-type vertex corrections with virtual W in the chiral limit
     */    
    complex FWa_0(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @return non-abelian-type vertex corrections with virtual W in the chiral limit
     */     
    complex FbarWa_0(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @return non-abelian-type vertex corrections with virtual W in the chiral limit
     */        
    complex FWn_0(const double s) const;
    
    /**
     * @param[in] s momentum-squared
     * @return abelian-type vertex corrections with virtual W and top-quark for Zbb
     */  
    complex FWa_t(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @return abelian-type vertex corrections with virtual W and top-quark for Zbb
     */      
    complex FbarWa_t(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @return non-abelian-type vertex corrections with virtual W and top-quark for Zbb
     */     
    complex FWn_t(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @return Unified form factor F_Z
     */      
    complex FZ(const double s) const;
    
    /**
     * @param[in] s momentum-squared
     * @param[in] l name of a lepton
     * @return Unified form factor F_W for the l-lbar channel
     */  
    complex FW(const double s, const StandardModel::lepton l) const;
    
    /**
     * @param[in] s momentum-squared
     * @param[in] q name of a quark, except for TOP
     * @return Unified form factor F_W for the q-qbar channel
     */  
    complex FW(const double s, const StandardModel::quark q) const;
    
    
    ////////////////////////////////////////////////////////////////////////        
    
    /**
     * @brief TEST function
     * @param[in] s momentum-squared
     * @param[in] mf the mass of the fermion in the loop
     * @return FWn_0 + FWn_t   
     */
    complex TEST_FWn(const double s, const double mf) const;
    
    
    ////////////////////////////////////////////////////////////////////////    

private:
    const EWSMcommon& EWSMC;    
    
    
};

#endif	/* ONELOOPEW_H */

