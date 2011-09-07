/* 
 * File:   EW.h
 * Author: mishima
 */

#ifndef EW_H
#define	EW_H

#include <ThObsType.h>
#include <StandardModel.h>
#include <EWSM.h>
#include <ZFitter.h>

using namespace gslpp;


class EW : public ThObsType {
public:

    /**
     * @brief EW constructor
     * @param[in] SM_i an object of StandardModel class
     */
    EW(const StandardModel& SM_i);

    /**
     * @brief EW copy constructor
     * @param[in] orig reference to an EW object
     */
    //EW(const EW& orig);
    
    /**
     * @brief EW destructor
     */
    virtual ~EW();

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @return a reference to the EWSM object in EW class
     */
    const EWSM& getEWSM() const { return this->myEWSM; }
    
    /**
     * @return a reference to the ZFitter object in EW class
     */
    const ZFitter& getZFitter() const { return this->myZFitter; }
    
    /**
     * @brief electromagnetic coupling alpha at Mz
     * @return alpha(Mz)
     */
    double getAlphaMZ() const {
        return alphaMZ;
    }

    /**
     * @brief the W bosn mass
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
     * @brief computes M_W and the effective weak couplings with EWSM class
     * @param[in] schemeMw resummation scheme for Mw
     * @param[in] schemeRhoZ resummation scheme for rho_Z^f
     * @param[in] schemeKappaZ resummation scheme for kappa_Z^f
     * @param[in] flag_order
     */
    void ComputeEWSM(const EWSM::schemes_EW schemeMw, 
                     const EWSM::schemes_EW schemeRhoZ,
                     const EWSM::schemes_EW schemeKappaZ,
                     const bool flag_order[EWSM::orders_EW_size]);
    
    /**
     * @brief computes M_W and the effective weak couplings with ZFITTER class
     * @param[in] schemeMw resummation scheme for Mw
     * @param[in] schemeRhoZ resummation scheme for rho_Z^f
     * @param[in] schemeKappaZ resummation scheme for kappa_Z^f
     * @param[in] flag_order
     */
    void ComputeZFitter(const EWSM::schemes_EW schemeMw, 
                        const EWSM::schemes_EW schemeRhoZ,
                        const EWSM::schemes_EW schemeKappaZ,
                        const bool flag_order[EWSM::orders_EW_size]);
    

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
     * @return the effective mixing angle for lepton "l"
     */
    double sin2thetaEff(const StandardModel::lepton l) const;
    
     /**
     * @param[in] q name of a quark
     * @return the effective mixing angle for quark "q"
     */
    double sin2thetaEff(const StandardModel::quark q) const;   
    
    /**
     * @param[in] l name of a lepton
     * @return the partial width of Z decay into an l\bar{l} pair 
     */
    double Gamma_l(const StandardModel::lepton l) const;
        
    /**
     * @param[in] q name of a quark
     * @return the partial width of Z decay into a q\bar{q} pair 
     */
    double Gamma_q(const StandardModel::quark q) const;
        
    /**
     * @return the partial width of Z decay into neutrinos
     */
    double Gamma_inv() const;

    /**
     * @return the hadronic width of the Z boson
     */
    double Gamma_had() const;

    /**
     * @return the total width of the Z boson
     */
    double Gamma_Z() const;
    
    /**
     * @param[in] l name of a lepton
     * @return the cross section for e^+e^- -> Z -> l\bar{l}
     */
    double sigma0_l(const StandardModel::lepton l) const;

    /**
     * @return the cross section e^+e^- -> Z -> hadrons
     */
    double sigma0_had() const; 
 
    /**
     * @param[in] l name of a lepton
     * @return asymmetry parameter for Z->l\bar{l}
     */
    double A_l(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of a quark
     * @return asymmetry parameter for Z->q\bar{q}
     */
    double A_q(const StandardModel::quark q) const;
    
    
    ////////////////////////////////////////////////////////////////////////     
    
private:
    EWSM myEWSM;
    ZFitter myZFitter;
    
    double alphaMZ;
    double Mw;
    complex rhoZ_l[6], rhoZ_q[6];
    complex kappaZ_l[6], kappaZ_q[6];
    
    double mcMz, mbMz;
    
};

#endif	/* EW_H */

