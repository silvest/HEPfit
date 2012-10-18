/* 
 * File:   LEP2TwoFermions.h
 * Author: mishima
 */

#ifndef LEP2TWOFERMIONS_H
#define	LEP2TWOFERMIONS_H

#include <StandardModel.h>
#include "EW.h"
using namespace gslpp;


/**
 * @class LEP2TwoFermions
 * @brief Cross sections and forward-backward asymmetries for e^+e^- -> f fbar at LEP-II
 */
    class LEP2TwoFermions : public EW {
public:

    // Radiative Corrections for the LEP-II observables
    enum LEP2RCs {Weak=0, WeakBox, ISR, QEDFSR, QCDFSR, NUMofLEP2RCs};
    
    /**
     * @brief LEP2TwoFermions constructor
     * @param[in] SM_i an object of StandardModel class
     */    
    LEP2TwoFermions(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] flags to control radiative corrections
     * @return the total cross section for e^+ e^- -> l lbar in GeV^{-2}
     */
    double sigma_l(const StandardModel::lepton l, const double s, 
                   const double Mw, const double GammaZ, const bool bRCs[]) const;
    
    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] flags to control radiative corrections
     * @return the total cross section for e^+ e^- -> q qbar in GeV^{-2}
     */
    double sigma_q(const StandardModel::quark q, const double s, 
                   const double Mw, const double GammaZ, const bool bRCs[]) const;

    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] flags to control radiative corrections
     * @return the forward-backward asymmetry for e^+ e^- -> l lbar
     */
    double AFB_l(const StandardModel::lepton l, const double s, 
                 const double Mw, const double GammaZ, const bool bRCs[]) const;
    
    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] flags to control radiative corrections
     * @return the forward-backward asymmetry for e^+ e^- -> q qbar
     */
    double AFB_q(const StandardModel::quark q, const double s, 
                 const double Mw, const double GammaZ, const bool bRCs[]) const;

    /**
     * @param x s'=(1-x)s
     * @param s the invariant mass squared of the initial-state e^+ e^- pair
     * @return the additive radiator of initial-state radiations in cross sections
     */
    double H_ISR(const double x, const double s) const;

    /**
     * @param x s'=(1-x)s
     * @param s the invariant mass squared of the initial-state e^+ e^- pair
     * @return the additive radiator of initial-state radiations in forward-backward asysmmetries
     */
    double H_ISR_FB(const double x, const double s) const;

    /**
     * @param[in] l name of a lepton
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] flags to control radiative corrections
     * @return the form factor beta_f^2*G_3(s) for e^+ e^- -> l lbar
     */
    double G_3prime_l(const StandardModel::lepton l, const double s, 
                      const double Mw, const double GammaZ, const bool bRCs[]) const;
    
    /**
     * @param[in] q name of a quark
     * @param[in] s invariant mass squared of the initial-state e^+ e^- pair
     * @param[in] Mw the W-boson mass 
     * @param[in] GammaZ the Z-boson decay width
     * @param[in] flags to control radiative corrections
     * @return the form factor beta_f^2*G_3(s) for e^+ e^- -> q qbar
     */
    double G_3prime_q(const StandardModel::quark q, const double s, 
                      const double Mw, const double GammaZ, const bool bRCs[]) const;    

    
    ////////////////////////////////////////////////////////////////////////     
private:

    /**
     * @param s invariant mass squared
     * @return the electromagnetic coupling at s
     */
    double alpha_at_s(const double s) const;
    
    /**
     * @param[in] l name of lepton
     * @return mass of lepton
     */
    double ml(const StandardModel::lepton l) const;

    /**
     * @param[in] q name of quark
     * @param[in] mu renormalization scale
     * @param[in] order (=LO, NLO, NNLO, FULLNLO[defalut], FULLNNLO)
     * @return the MSbar mass of u, d, s, c, b or the pole mass of t
     */
    double mq(const StandardModel::quark q, const double mu, 
              const orders order=FULLNLO) const;
    
    double sigma(const double s, const double Mw, const double GammaZ, 
                 const double I3f, const double Qf, const double mfp,
                 const double mf, const double Ncf, const bool bRCs[]) const;
    
    double AFB(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mfp,
               const double mf, const bool bRCs[]) const;
    
};

#endif	/* LEP2TWOFERMIONS_H */

