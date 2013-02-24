/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EW_CHMN_H
#define	EW_CHMN_H

#include <gslpp.h>
#include <StandardModel.h>
using namespace gslpp;

/**
 * @class EW_CHMN
 * @ingroup EW 
 * @brief A test class for the electroweak precision observables. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details A class for the electroweak precision observables with a simple 
 * parametrization in JHEP1111 (2011) 068 [arXiv:1104.1769] by G-C. Cho, 
 * K. Hagiwara, Y. Matsumoto and D. Nomura, which is designed to be valid 
 * in the region 100 GeV <= mHl <= 1 TeV. 
 * The NP contributions to the Fermi constant G_F, to the effective couplings 
 * g_L^f and g_R^f except for f=b, to R_Z, and to R_W are assumed to be 0. 
 * See also NPB574 (2000) 623 [hep-ph/9912260] by G-C. Cho and K. Hagiwara. 
 * This class is used for a test. 
 */
class EW_CHMN {
public:

    /**
     * @brief EW_CHMN constructor
     * @param[in] SM_i an object of StandardModel class
     */
    EW_CHMN(const StandardModel& SM_i);

    
    //////////////////////////////////////////////////////////////////////// 
    
    double Mw() const;
    double GammaW() const;
    
    double DeltaSz() const;
    double DeltaTz() const;
    
    // Effective couplings
    double gL_l(const StandardModel::lepton l) const;
    double gR_l(const StandardModel::lepton l) const;
    double gL_q(const StandardModel::quark q) const;
    double gR_q(const StandardModel::quark q) const;

    // Z-boson partial width
    double GammaZ_l(StandardModel::lepton l) const;
    double GammaZ_q(StandardModel::quark q) const;
    
    // Z-boson hadronic width
    double GammaZ_had() const;
    
    // Z-boson total width
    double GammaZ() const;    
    
    double R_l(const StandardModel::lepton l) const;
    double R_c() const;    
    double R_b() const;    

    // hadronic peak cross section
    double sigma0_had() const;
    
    // left-right asymmetry parameter
    double A_l(const StandardModel::lepton l) const;
    double A_q(const StandardModel::quark q) const;
    
    // forward-backward asymmetry
    double AFB_l(const StandardModel::lepton l) const;
    double AFB_q(const StandardModel::quark q) const;
    
    // effective weak mixing angle
    double sin2thetaEff() const;
    
    /**
     * @return the oblique parameters S
     */
    double S() const 
    {
        return ( SM.obliqueS() );
    }
    
    /**
     * @return the oblique parameters T
     */    
    double T() const 
    {
        return ( SM.obliqueT() );
    }
    
    /**
     * @return the oblique parameters U
     */    
    double U() const 
    {
        return ( SM.obliqueU() );
    }
    
    
    ////////////////////////////////////////////////////////////////////////     

private:
    const StandardModel& SM;
    
    double x_alpha() const;
    double x_t() const;
    double x_h() const;
    double x_s() const;
    
    double DeltaRz() const;
    double DeltaRw() const;    
    
    double DeltaMw() const;

    double Delta_gbarZ2() const;
    double Delta_sbar2() const;
    
    // SM contributions to S_Z, T_Z, M_W, R_Z and R_W
    double DeltaSz_SM() const;
    double DeltaTz_SM() const;
    double DeltaMw_SM() const;
    double DeltaRz_SM() const;
    double DeltaRw_SM() const;
    
    // the color factors, including mass and QCD corrections
    double CV_l(StandardModel::lepton l) const;
    double CV_q(StandardModel::quark q) const;
    double CA_l(StandardModel::lepton l) const;
    double CA_q(StandardModel::quark q) const;
    
    // corrections from the imaginary part of loop-induced mixing between the photon and the Z boson
    double deltaImKappa_l(StandardModel::lepton l) const;
    double deltaImKappa_q(StandardModel::quark q) const;
    
    // non-factorizable mixed EW/QCD corrections
    double DeltaEWQCD_l(StandardModel::lepton l) const;
    double DeltaEWQCD_q(StandardModel::quark q) const;
    
    
};

#endif	/* EW_CHMN_H */

