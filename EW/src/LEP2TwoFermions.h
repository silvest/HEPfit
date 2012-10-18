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

    
    ////////////////////////////////////////////////////////////////////////     
private:

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
              const orders order=FULLNNLO) const;
    
    double sigma(const double s, const double Mw, const double GammaZ, 
                 const double I3f, const double Qf, const double mfp,
                 const double mf, const double Ncf, const bool bRCs[]) const;
    
    double AFB(const double s, const double Mw, const double GammaZ, 
               const double I3f, const double Qf, const double mfp,
               const double mf, const bool bRCs[]) const;
    
};

#endif	/* LEP2TWOFERMIONS_H */

