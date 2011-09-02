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
     * @brief scheme for the resummations in Delta r, rho_Z^f and kappa_Z^f
     */
    enum schemes_EW {NORESUM=0, OMSI, INTERMEDIATE, OMSII, schemes_EW_size};    
    
    
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
     * @return a reference to an EWSM object
     */
    const EWSM& getEWSM() const { return this->myEWSM; }
    
    /**
     * @return a reference to a ZFitter object
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
     * @brief computes M_W and the effective weak couplings
     */
    void ComputeEWSM();
    
    /**
     * @brief computes M_W and the effective weak couplings with ZFITTER
     */
    void ComputeZFitter();
    
    
    ////////////////////////////////////////////////////////////////////////     
    
private:
    EWSM myEWSM;
    ZFitter myZFitter;
    
    double alphaMZ;
    double Mw;
    complex rhoZ_l[6], rhoZ_q[6];
    complex kappaZ_l[6], kappaZ_q[6];
    
};

#endif	/* EW_H */

