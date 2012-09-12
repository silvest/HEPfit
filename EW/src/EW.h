/* 
 * File:   EW.h
 * Author: mishima
 */

#ifndef EW_H
#define	EW_H

#include <ThObsType.h>
#include <StandardModel.h>

using namespace gslpp;


class EW : public ThObsType {
public:
    
    /**
     * @brief EW constructor
     * @param[in] SM_i an object of StandardModel class
     * @param[in] bDebug_i boolean value for debugging (true for debugging)
     */
    EW(const StandardModel& SM_i, const bool bDebug_i=false);

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @return boolean: true for the case where S, T and U are employed for the current model
     */
    bool checkModelForSTU() const;
    
    /**
     * @return a reference to the StandardModel object in the current class
     */
    const StandardModel& getSM() const {
        return SM;
    } 
        
    /**
     * @param[in] l lepton
     * @return electric charge of a lepton "l"
     */
    double Ql(const StandardModel::lepton l) const;    

    /**
     * @param[in] q quark
     * @return electric charge of a quark "q"
     */
    double Qq(const StandardModel::quark q) const;    

    /**
     * @return sin^2\theta_W without weak corrections, but with \alpha(Mz^2)
     */
    double s2() const;

    /**
     * @return cos^2\theta_W without weak corrections, but with \alpha(Mz^2)
     */
    double c2() const;

    /**
     * @return the oblique parameters S
     */
    double S() const {
        return ( SM.obliqueS() );
    }
    
    /**
     * @return the oblique parameters T
     */    
    double T() const {
        return ( SM.obliqueT() );
    }
    
    /**
     * @return the oblique parameters U
     */    
    double U() const {
        return ( SM.obliqueU() );
    }

    /**
     * @return the oblique parameters V
     */    
    double V() const {
        return ( SM.obliqueV() );
    }
    
    /**
     * @return the oblique parameters W
     */    
    double W() const {
        return ( SM.obliqueW() );
    }
    
    /**
     * @return the oblique parameters X
     */    
    double X() const {
        return ( SM.obliqueX() );
    }
    
    /**
     * @return the oblique parameters Y
     */    
    double Y() const {
        return ( SM.obliqueY() );
    }    

    
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

    
    double dsigma_lLEP2(const StandardModel::lepton l,const double s,const double W,
                               const double X,const double Y, const double cos_theta) const;
    
    double dsigma_qLEP2(const QCD::quark q,const double s,const double W,
                               const double X,const double Y, const double cos_theta) const;
   
    ////////////////////////////////////////////////////////////////////////     
    
private:
    bool bDebug; // true for debugging    
    
    const StandardModel& SM;
    
};

#endif	/* EW_H */

