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

    /**
     * @brief OneLoopEW copy constructor
     * @param[in] orig reference to a OneLoopEW object
     */
    //OneLoopEW(const OneLoopEW& orig);

    /**
     * @brief OneLoopEW destructor
     */
    virtual ~OneLoopEW();


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
    
    
    ////////////////////////////////////////////////////////////////////////    
    
private:
    const EWSMcommon& EWSMC;

    
    ////////////////////////////////////////////////////////////////////////    
       
    /**
     * @param[in] s momentum-squared
     * @return the bosonic contribution to the self-energy function of the W boson for mu=Mz
     */
    complex SigmaWW_bos(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @return the bosonic contribution to the self-energy function of the Z boson for mu=Mz
     */
    complex SigmaZZ_bos(const double s) const;    

    /**
     * @param[in] s momentum-squared
     * @return the bosonic contribution to the self-energy function of the photon for mu=Mz
     */
    complex PiGammaGamma_bos(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @return the bosonic contribution to the self-energy function of the Z-gamma mixing for mu=Mz
     */
    complex PiZgamma_bos(const double s) const;
 
    /**
     * @param[in] s momentum-squared
     * @return the fermionic contribution to the self-energy function of the W boson for mu=Mz
     */
    complex SigmaWW_fer(const double s) const;
    
    /**
     * @param[in] s momentum-squared
     * @return the fermionic contribution to the self-energy function of the Z boson for mu=Mz
     */
    complex SigmaZZ_fer(const double s) const;
  
    /**
     * @param[in] s momentum-squared
     * @return the fermionic contribution to the self-energy function of the photon for mu=Mz
     */
    complex PiGammaGamma_fer(const double s) const;

    /**
     * @param[in] s momentum-squared
     * @return the fermionic contribution to the self-energy function of the Z-gamma mixing for mu=Mz
     */    
    complex PiZgamma_fer(const double s) const;

    
    ////////////////////////////////////////////////////////////////////////   
    
    /**
     * @brief TEST function
     * @return Delta Rhobar^{bos,F} for mu=Mw
     * @attention The renormalization scale is fixed to be mu=Mw.
     */
    double TEST_DeltaRhobar_bos() const;

    
    /**
     * @return Delta Rhobar_W^{bos,F}
     */
    double DeltaRhobarW_bos() const;    


};

#endif	/* ONELOOPEW_H */

