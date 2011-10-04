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
     * @brief schemes for the resummations in Mw, rho_Z^f and kappa_Z^f
     * 
     * APPROXIMATEFORMULA for the use of approximate formulae
     */
    enum schemes_EW {NORESUM=0, OMSI, INTERMEDIATE, OMSII, APPROXIMATEFORMULA}; 

    
    //////////////////////////////////////////////////////////////////////// 
    
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
    double getAlphaMz() const {
        return alphaMz;
    }

    /**
     * @return the total radiative corrections to alpha at Mz
     */
    double getDeltaAlpha() const {
        return DeltaAlpha;
    }

    /**
     * @return the leptonic and 5-flavour hadronic corrections to alpha at Mz
     */
    double getDeltaAlpha_l5q() const {
        return DeltaAlpha_l5q;
    }     
    
    /**
     * @brief the W bosn mass
     * @return M_W
     */
    double getMw() const {
        return Mw;
    }    
    
    /**
     * @return s_W^2 = 1 - Mw^2/Mz^2
     */    
    double getSW2() const {
        return sW2;
    }
    
    /**
     * @return c_W^2 = Mw^2/Mz^2 
     */
    double getCW2() const {
        return cW2;
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
    virtual void ComputeEWSM(const schemes_EW schemeMw, 
                             const schemes_EW schemeRhoZ,
                             const schemes_EW schemeKappaZ,
                             const bool flag_order[EWSM::orders_EW_size]);
    
    /**
     * @brief computes M_W and the effective weak couplings with ZFITTER class
     * @param[in] schemeMw resummation scheme for Mw
     * @param[in] schemeRhoZ resummation scheme for rho_Z^f
     * @param[in] schemeKappaZ resummation scheme for kappa_Z^f
     * @param[in] flag_order
     */
    virtual void ComputeZFitter(const schemes_EW schemeMw, 
                                const schemes_EW schemeRhoZ,
                                const schemes_EW schemeKappaZ,
                                const bool flag_order[EWSM::orders_EW_size]);
    

    //////////////////////////////////////////////////////////////////////// 
    
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
    ApproximateFormulae* myApproximateFormulae;
        
    double Mw_error;    
    
    double alphaMz, DeltaAlpha, DeltaAlpha_l5q;
    double Mw, sW2, cW2;
    complex rhoZ_l[6], rhoZ_q[6];
    complex kappaZ_l[6], kappaZ_q[6];
    
    double mcMz, mbMz;

    
    ////////////////////////////////////////////////////////////////////////     

    /**
     * @param[in] schemeMw
     * @return resummed Mw
     */
    double resumMw(const schemes_EW schemeMw);
    
    /**
     * @param[in] schemeRhoZ resummation scheme for rho_Z^f
     * @param[in] deltaRho_rem 
     * @return resummed Re[rho_Z^f]
     */
    double resumRhoZ(const schemes_EW schemeRhoZ, 
                     const double deltaRho_rem[EWSM::orders_EW_size]);
    
    /**
     * @param[in] schemeMw resummation scheme for Mw
     * @param[in] deltaKappa_rem 
     * @return resummed Re[kappa_Z^f]
     */
    double resumKappaZ(const schemes_EW schemeKappaZ,
                       const double deltaKappa_rem[EWSM::orders_EW_size]);
        
    /**
     * @brief sets flags for ZFITTER
     * @param[in] schemeMw resummation scheme for Mw
     * @param[in] schemeRhoZ resummation scheme for rho_Z^f
     * @param[in] schemeKappaZ resummation scheme for kappa_Z^f
     * @param[in] flag_order
     */
    void SetZFitterFlags(const schemes_EW schemeMw, 
                         const schemes_EW schemeRhoZ,
                         const schemes_EW schemeKappaZ,
                         const bool flag_order[EWSM::orders_EW_size]);
    
};

#endif	/* EW_H */

