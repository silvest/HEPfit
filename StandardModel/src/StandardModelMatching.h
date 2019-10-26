/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef STANDARDMODELMATCHING_H
#define	STANDARDMODELMATCHING_H

#include "ModelMatching.h"
#include "gslpp.h"

#define LEPS 1.e-5 // tolerance in the limit of S(x,y) to S(x)

class StandardModel;

/**
 * @class StandardModelMatching
 * @ingroup StandardModel
 * @brief A class for the matching in the Standard Model. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class StandardModelMatching : public ModelMatching {
public:

    StandardModelMatching(const StandardModel & SM_i);
    
    virtual ~StandardModelMatching();
    
    /**
     *
     * @brief Updates to new Standard %Model parameter sets.
     */
    void updateSMParameters();
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{d} \f$ 
     * @return return the vector of SM Wilson coefficients
     */
    virtual  std::vector<WilsonCoefficient>& CMdbd2() ;
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{s} \f$ 
     * @return return the vector of SM Wilson coefficients
     */
    virtual   std::vector<WilsonCoefficient>& CMdbs2() ;

    /**
     * 
     * @brief \f$ \Delta C = 2 \f$,
     * @return return the vector of SM Wilson coefficients
     */
    virtual   std::vector<WilsonCoefficient>& CMdd2() ;
    
    /**
     * 
     * @brief \f$ \Delta S = 2 \f$ 
     * @return return the vector of SM Wilson coefficients
     */
    virtual   std::vector<WilsonCoefficient>& CMdk2() ;
        
    /**
     * 
     * @brief operator basis: 
     * @brief current current; qcd penguins; magnetic and chromomagnetic penguins; semileptonic  
     * @param[in] a a=0 deltaS=0 deltaC=0;  a=1 deltaS=1 deltaC=0;
     * @return Wilson coefficients Buras base for non-leptonic B decays 
     */
    virtual   std::vector<WilsonCoefficient>& CMbnlep( int a) ;
    
    /**
     * 
     * @brief operator basis: - current current opertors  
     * @param[in] a a=0 deltaS=0 deltaC=0;  a=1 1,0 ;  a=2 0,1 ; a=3 1,1
     * @return Wilson coefficients, Buras basis, for non-leptonic B decays 
     */
    virtual   std::vector<WilsonCoefficient>& CMbnlepCC( int a) ;
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbsg() ;
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    virtual   std::vector<WilsonCoefficient>& CMprimebsg() ;
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow K^* l^{+} l^{-} \f$
     */
    virtual   std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton) ;
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow K^* l^{+} l^{-} \f$
     */
    virtual   std::vector<WilsonCoefficient>& CMprimeBMll(QCD::lepton lepton) ;
    
    /** 
     * 
     * @brief operator basis: 
     * @brief current-current; qcd penguins;
     * @brief hep/ph 9512380  
     * @return Wilson coefficients Buras base for \f$ K \rightarrow \pi \pi \f$ decays
     */
    virtual  std::vector<WilsonCoefficient>& CMK() ;
    
    /** 
     * 
     * @brief operator basis: 
     * @brief current-current (open up - open charm) 
     * @brief hep/ph 9512380 
     * @return Wilson coefficients Buras base for \f$ B \rightarrow \pi \pi \f$ decays
     */
    virtual  std::vector<WilsonCoefficient>& CMKCC() ;
    
    /**
     * 
     * @brief current-current oerators, Misiak basis
     * @return Wilson coefficients for \f$ D^{0} \rightarrow \pi \pi , K K \f$
     */
    virtual   std::vector<WilsonCoefficient>& CMd1() ;
    
    /**
     * 
     * @brief current-current oerators, Buras basis
     * @return Wilson coefficients for \f$ D^{0} \rightarrow \pi \pi , K K \f$
     */
    virtual   std::vector<WilsonCoefficient>& CMd1Buras() ;
    
    /**
     * 
     * @return Wilson coefficients for \f$ K_{L} \rightarrow \pi \nu \nu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMkpnn() ;
    /**
     * 
     * @return Wilson coefficient for \f$ K^{+} \rightarrow \mu \mu \f$, short distance top contribution
     */
    virtual  std::vector<WilsonCoefficient>& CMkmm() ;
    
    /**
     * @brief Wilson coefficient for \f$ B_{s} \rightarrow \mu \mu \f$
     * @return Wilson coefficient for \f$ B_{s} \rightarrow \mu \mu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbsmm() ;
    
    /**
     * @brief Wilson coefficient for \f$ B_{d} \rightarrow \mu \mu \f$
     * @return Wilson coefficient for \f$ B_{d} \rightarrow \mu \mu \f$
     */
    virtual   std::vector<WilsonCoefficient>& CMbdmm() ;
    
    /**
     * 
     * @return Wilson coefficient for \f$ B \rightarrow \tau \nu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbtaunu(QCD::meson meson_i) ;
    
    /**
     * 
     * @return Wilson coefficients for \f$ B_{s} \rightarrow X_{s} \nu \nu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMBXsnn() ;
    
    /**
     * 
     * @return Wilson coefficients for \f$ B_{d} \rightarrow X_{d} \nu \nu \f$
     */
    virtual   std::vector<WilsonCoefficient>& CMBXdnn() ;
    
     /**
     * 
     * @return Wilson coefficients for \f$ \ell_i \to \ell_j\f$
     */
    virtual  std::vector<WilsonCoefficient>& CMDLij(int li_lj) ;

     /**
     * 
     * @return Wilson coefficients for \f$ \ell_i \to \ell_j\ell_j\ell_j\f$
     */
    virtual  std::vector<WilsonCoefficient>& CMDLi3j(int li_lj) ;

     /**
     * 
     * @return Wilson coefficients for \f$ \mu \to e\f$ conversion in nuclei
     */
    virtual  std::vector<WilsonCoefficient>& CMmueconv() ;

     /**
     * 
     * @return Wilson coefficients for \f$ (g-2)_\mu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMgminus2mu() ;

    double x_c(const double mu, const orders order = FULLNNLO) const;
    double x_t(const double mu, const orders order = FULLNNLO) const;
    double mt2omh2(const double mu, const orders order = FULLNNLO) const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double A0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the semileptonic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double B0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double C0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double D0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the chromomagnetic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double F0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the chromomagnetic operator
     * in the Misiak basis, Chetyrkin et al hep-ph/9612313 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double E0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the semileptonic operator
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/9901278v1
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double A1t(double x, double mu) const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the semileptonic operator
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/9901278v1
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double B1t(double x, double mu)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/9901278v1
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double C1t(double x, double mu)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/9901278v1
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double D1t(double x, double mu)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the semileptonic operator
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/9901278v1
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double F1t(double x, double mu) const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/9910220
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double E1t(double x, double mu) const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/9910220
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double G1t(double x, double mu) const;
    
     /**
     * 
     * @brief power expansion appearing in the Wilson coefficient 
     * C7 at NNLO, see Misiak and Steinhauser hep-ph/0401041
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double C7c_3L_at_mW(double x) const;
    
     /**
     * 
     * @brief power expansion appearing in the Wilson coefficient 
     * C7 at NNLO, see Misiak and Steinhauser hep-ph/0401041
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double C7t_3L_at_mt(double x) const;
    
     /**
     * 
     * @brief loop function which appear in the Wilson coefficient 
     * C7 at NNLO, see Misiak and Steinhauser hep-ph/0401041
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double C7t_3L_func(double x, double mu) const;
    
     /**
     * 
     * @brief power expansion appearing in the Wilson coefficient 
     * C8 at NNLO, see Misiak and Steinhauser hep-ph/0401041
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double C8c_3L_at_mW(double x) const;
    
     /**
     * 
     * @brief power expansion appearing in the Wilson coefficient 
     * C8 at NNLO, see Misiak and Steinhauser hep-ph/0401041
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double C8t_3L_at_mt(double x) const;
    
     /**
     * 
     * @brief loop function which appear in the Wilson coefficient 
     * C8 at NNLO, see Misiak and Steinhauser hep-ph/0401041
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double C8t_3L_func(double x, double mu) const;

    /**
     * 
     * @brief loop function which appear in the Wilson coefficient
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/9910220
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double Tt(double x) const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/0512066
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double Wt(double x) const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/0512066
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double Eet(double x) const;
    
    /**
     * 
     * @brief approximation of two-loops EW correction for Q_10 operator
     * in the non-effective Misiak basis, Misiak and Urban hep-ph/1311.1348
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double Rest(double x, double mu) const;

    /**
     * 
     * @brief approximation of two-loops EW correction for Q_10 operator from 1311.1348
     * @param[in] x the ratio \f$ \frac{m_t^2}{M_W^2} \f$
     * @param[in] mu the matching scale
     * @return return the fit to the 22 matching of C10 in the OS1 scheme
     */
    double C10_OS1(double x, double mu);
    
    /**
     *
     * @param[in] x the square of the ratio between top mass and W mass
     * @return 0th order loop function for the top contribution to K_L -> mu mu decays
     */
    double Y0(double x)const;
    
    /**
     *
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     * @return first order in QCD loop function for the top contribution to K_L -> mu mu decays
     */
    double Y1(double x, double mu)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the effective Misiak basis, LO term, Chetyrkin et al hep-ph/9612313 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double C7LOeff(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the effective Misiak basis, NLO term, Chetyrkin et al hep-ph/9612313 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double C7NLOeff(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the chromomagnetic operator
     * in the effective Misiak basis, LO term, Chetyrkin et al hep-ph/9612313 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double C8LOeff(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the chromomagnetic operator
     * in the effective Misiak basis, LO term, Chetyrkin et al hep-ph/9612313 
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double C8NLOeff(double x)const;

    /**
     * 
     * @brief loop functions for non-leptonic B decays, Buiras Basis
     * Buras et al, hep-ph/9512380v1
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double B0b(double x)const;
    
    /**
     * 
     * @brief loop functions for non-leptonic B decays, Buiras Basis
     * Buras et al, hep-ph/9512380v1
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double C0b(double x)const;
    
    /**
     * 
     * @brief loop functions for non-leptonic B decays, Buiras Basis
     * Buras et al, hep-ph/9512380v1
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double D0b(double x)const;
    
    /**
     * 
     * @brief loop functions for non-leptonic B decays, Buras Basis
     * Buras et al, hep-ph/9512380
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double D0b_tilde(double x)const;
    
    /**
     * 
     * @brief loop functions for non-leptonic B decays, Buiras Basis
     * Buras et al, hep-ph/9512380v1
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double E0b(double x)const;
    
    /**
     * 
     * @brief loop functions for QED corrections of Delta F = 1 decays
     * Buras, Gambino, Haisch, hep-ph/9911250
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double B1d(double x, double mu)const;
    
    /**
     * 
     * @brief loop functions for QED corrections of Delta F = 1 decays
     * Buras, Gambino, Haisch, hep-ph/9911250
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double B1d_tilde(double x, double mu)const;
    
    /**
     * 
     * @brief loop functions for QED corrections of Delta F = 1 decays
     * Buras, Gambino, Haisch, hep-ph/9911250
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double B1u(double x, double mu)const;
    
    /**
     * 
     * @brief loop functions for QED corrections of Delta F = 1 decays
     * Buras, Gambino, Haisch, hep-ph/9911250
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double B1u_tilde(double x, double mu)const;
    
    /**
     * 
     * @brief loop functions for QED corrections of Delta F = 1 decays
     * Buras, Gambino, Haisch, hep-ph/9911250
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double C1ew(double x)const;
    
    /**
     * 
     * @brief auxiliary loop function for QED corrections of Delta F = 1 decays
     * Buras, Gambino, Haisch, hep-ph/9911250
     * @param[in] xt the square of the ratio between top mass and W mass
     * @param[in] xz the square of the ratio between Z mass and W mass
     */
    double Zew(double xt, double xz)const;
        
    /**
     * 
     * @brief loop functions for QED corrections of Delta F = 1 decays
     * Buras, Gambino, Haisch, hep-ph/9911250
     * @param[in] xt the square of the ratio between top mass and W mass
     * @param[in] xz the square of the ratio between Z mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double Gew(double xt, double xz, double mu)const;
    
    /**
     * 
     * @brief loop functions for QED corrections of Delta F = 1 decays
     * Buras, Gambino, Haisch, hep-ph/9911250
     * @param[in] xt the square of the ratio between top mass and W mass
     * @param[in] xz the square of the ratio between Z mass and W mass
     * @param[in] mu the matching scale of the Wilson coefficients
     */
    double Hew(double xt, double xz, double mu)const;
    
    /**
     *  
     * @brief hep-ph/9512380
     * @return the loop function for the charm-charm contribution to the Delta S = 2 effective hamiltonian multiplied by the CKM element 
     */
    gslpp::complex S0c() const;
    
    /**
     *  
     * @brief hep-ph/9512380
     * @return the loop function for the charm-top contribution to the Delta S = 2 effective hamiltonian multiplied by the CKM element 
     */
    gslpp::complex S0ct() const;
    
    /**
     *  
     * @brief hep-ph/9512380v1
     * @return the loop function for the top-top contribution to the Delta S = 2 effective hamiltonian
     */
    gslpp::complex S0tt() const;
    
    /**
     *
     * @brief hep-ph/9512380v1
     * @param[in] x the square of the ratio between top mass and W mass
     * @return 0th order loop function for the top contribution to K -> pi nu nu decays
     */
    double X0t(double x)const;
    
    /**
     *
     * @brief hep-ph/1009.0947v2
     * @param[in] x the square of the ratio between top mass and W mass
     * @return first order in QCD loop function for the top contribution to K -> pi nu nu decays
     */
    double X1t(double x)const;
    
    /**
     *
     * @brief hep-ph/1009.0947v2
     * @param[in] x the square of the ratio between top mass and W mass
     * @param[in] a
     * @param[in] mu
     * @return two loop EW loop function for the top contribution to K -> pi nu nu decays
     * in the limit of small theta_W
     */
    double Xewt(double x, double a, double mu)const;

    /**
     *
     * @param[in] x the square of the ratio between top mass and W mass
     * @return
     */

    double S1(double x) const;

    virtual double S0(double, double) const;

    /*
     * Wilson coefficients Misiak basis
     * Operator block: L (2)
     * Normalization: 4 G_F / sqrt(2) x CKM
     */
    virtual std::vector<WilsonCoefficientNew>& CMDF1(std::string blocks, unsigned int nops);

    double getMt_mut() const {
        return Mt_mut;
    }
     
protected:
    std::vector<WilsonCoefficient> vmcdb, vmcds, vmcd2, vmck2, vmck, vmckcc;
    std::vector<WilsonCoefficient> vmcbsg, vmcprimebsg, vmcBMll, vmcprimeBMll, vmcbnlep, vmcbnlepCC, vmcd1, vmcd1Buras;
    std::vector<WilsonCoefficient> vmckpnn, vmckmm, vmcbsnn, vmcbdnn, vmcbsmm, vmcbdmm, vmcbtaunu;
    std::vector<WilsonCoefficient> vmcDLij, vmcDLi3j, vmcmueconv, vmcgminus2mu;
    std::vector<WilsonCoefficientNew> vmcDF1;

    
private:    
    const StandardModel & SM; ///< An object of the %StandardModel class.
    double S0(double) const;
    double S0p(double x) const;
    double S11(double x) const;
    double S18(double x) const;
    WilsonCoefficient mcdbd2, mcdbs2, mcdd2, mcdk2, mck, mckcc;
    WilsonCoefficient mcbsg, mcprimebsg, mcBMll, mcprimeBMll, mcbnlep, mcbnlepCC, mcd1, mcd1Buras;
    WilsonCoefficient mckpnn, mckmm, mcbsnn, mcbdnn, mcbsmm, mcbdmm, mcbtaunu;
    WilsonCoefficient mcDLij, mcDLi3j, mcmueconv, mcgminus2mu;
    WilsonCoefficientNew mcC, mcP, mcM, mcL, mcQ, mcB;
    
    double Mut, Muw, Ale, GF, Mw_tree, Nc, CF, Mt_muw, Mt_mut;
    double gamma0, J5, BtNDR, Mw, sW2, mu_b;
    double L, Lz;
    double alstilde, aletilde;
    
    gslpp::complex lam_t;
    gslpp::matrix<gslpp::complex> Vckm;

    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order
     * @return return the value of the wilson coefficients for \f$ B_{s} \rightarrow  l^{+} l^{-} \f$
     */
    double setWCBsmm(int i, double x, orders order);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order_qed
     * @return return the electroweak value of the wilson coefficients for \f$ B_{s} \rightarrow  l^{+} l^{-} \f$
     */
    double setWCBsmmEW(int i, double x, orders_qed order_qed);
    
     /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order
     * @return return the value of the wilson coefficients for \f$ B_{d} \rightarrow  l^{+} l^{-} \f$
     */
    double setWCBdmm(int i, double x, orders order);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order_qed
     * @return return the electroweak value of the wilson coefficients for \f$ B_{d} \rightarrow  l^{+} l^{-} \f$
     */
    double setWCBdmmEW(int i, double x, orders_qed order_qed);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order
     * @return return the value of the wilson coefficients for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    double setWCbsg (int i, double x, orders order);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order
     * @return return the value of the wilson coefficients for \f$ B \rightarrow k^* l^{+} l^{-} \f$
     */
    double setWCBMll (int i, double x, orders order);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order
     * @return return the value of the QCD contribution to the Wilson coefficients for non-leptonic B decays
     */
    double setWCbnlep (int i, double x, orders order);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @return return the value of the electroweak contribution to the Wilson coefficients for non-leptonic B decays
     */
    double setWCbnlepEW (int i, double x);

    /**
     * 
     * @param z
     * @return two loop EW loop functions for K-> P nu nu, hep-ph/1009.0947v2 
     */
    double phi1(double z) const;
    
    /**
     * 
     * @param x
     * @param y
     * @return two loop EW loop functions for K-> P nu nu, hep-ph/1009.0947v2 
     */
    double phi2 (double x, double y) const;
    
    double CWbsgArrayLO[8], CWbsgArrayNLO[8], CWbsgArrayNNLO[8];
    double CWprimebsgArrayLO[8], CWprimebsgArrayNLO[8];
    double CWBMllArrayLO[19], CWBMllArrayNLO[19];
    double CWD1ArrayLO[10], CWD1ArrayNLO[10];
    double CWbnlepArrayLOqcd[10], CWbnlepArrayNLOqcd[10];
    double CWbnlepArrayLOew[10], CWbnlepArrayNLOew[10];
    
    double CWBsmmArrayNNLOqcd[8], CWBsmmArrayNLOqcd[8], CWBsmmArrayLOqcd[8];
    double CWBsmmArrayNLOewt4[8], CWBsmmArrayNLOewt2[8], CWBsmmArrayNLOew[8];
    
    double CWBdmmArrayNNLOqcd[8], CWBdmmArrayNLOqcd[8], CWBdmmArrayLOqcd[8];
    double CWBdmmArrayNLOewt4[8], CWBdmmArrayNLOewt2[8], CWBdmmArrayNLOew[8];
    
    double sw, swa, swb, swc, swd, swe, swf; //sen(theta_W) tree level
    double xcachea, xcacheb, xcachec, xcached, xcachee, xcachef; // caching

    
    /**
     * 
     * @param z
     * @return two loop EW loop functions, hep-ph/9911250 
     */
    double phi_z(double z) const;
    
    /**
     * 
     * @param x
     * @param y
     * @return two loop EW loop functions, hep-ph/9911250 
     */
    double phi_xy (double x, double y) const;
    
    /**
     * 
     * @param x \f$ \frac{m_t^2}{M_W^2} \f$
     * @return function appearing in the matching condition for \f$ C_3 \f$ 
     */
    double C3funNNLO(double x);
    
    /**
     * 
     * @param x \f$ \frac{m_t^2}{M_W^2} \f$
     * @param ord NLO or NNLO
     * @return function appearing in the matching condition for \f$ C_4 \f$ 
     */
    double C4fun(double x, orders ord);

    /**
     * 
     * @param x \f$ \frac{m_t^2}{M_W^2} \f$
     * @return function appearing in the matching condition for \f$ C_5 \f$ 
     */
    double C5funNNLO(double x);

    /**
     * 
     * @param x \f$ \frac{m_t^2}{M_W^2} \f$
     * @return function appearing in the matching condition for \f$ C_6 \f$ 
     */
    double C6funNNLO(double x);
    
    /**
     * 
     * @param x \f$ \frac{m_t^2}{M_W^2} \f$
     * @return function appearing in the matching condition for \f$ C_7 \f$ 
     */
    double C7funLO(double x);
    
    /**
     * 
     * @param x \f$ \frac{m_t^2}{M_W^2} \f$
     * @return function appearing in the matching condition for \f$ C_8 \f$ 
     */
    double C8funLO(double x);

    unsigned int setCMDF1(WilsonCoefficientNew& CMDF1, WilsonCoefficientNew& DF1block, unsigned int tot, schemes scheme, qcd_orders order_qcd, qed_orders order_qed);

    /**
     * 
     * @brief auxiliary function for mc_L() 
     * @param[in] x \f$ \frac{M_H^2}{m_t^2} \f$
     * @return return an approximation of the function \f$ f(x) \f$ from hep-ph/9707243 in the range x=(0,1)
     */    
    double fbb(double x);
    
    /**
     * 
     * @brief auxiliary function for mc_L() 
     * @param[in] x \f$ \frac{M_H^2}{m_t^2} \f$
     * @return return the function  \f$ g(x) \f$ from hep-ph/9707243 
     */    
    double gbb(double x);

    /**
     * 
     * @brief auxiliary function for mc_L() 
     * @param[in] x \f$ \frac{M_H^2}{m_t^2} \f$
     * @return return the function  \f$ \tau_b^{(2)}(x) \f$ from hep-ph/9707243 
     */    
    double taub2(double x);

    /**
     * 
     * @brief auxiliary function for mc_L()
     * @param[in] mu matching scale 
     * @param[in] x \f$ \frac{M_H^2}{m_t^2} \f$
     * @return return the function  \f$ \Delta_t(\mu, x) \f$ from hep-ph/9707243 
     */    
    double Delta_t(double mu, double x);

    /*
     * Wilson coefficients Misiak basis
     * Operator block: C (2)
     * Normalization: 4 G_F / sqrt(2) x CKM
     */
    WilsonCoefficientNew& mc_C();

    /*
     * Wilson coefficients Misiak basis
     * Operator block: P (4)
     * Normalization: 4 G_F / sqrt(2) x CKM
     */
    WilsonCoefficientNew& mc_P();

    /*
     * Wilson coefficients Misiak basis
     * Operator block: M (2)
     * Normalization: 4 G_F / sqrt(2) x CKM
     * QED only available at NLO and in approximate formulas
     * QED ref.: Gambino, Haisch, JHEP 0110, 020, hep-ph/0109058
     */
    WilsonCoefficientNew& mc_M();

    /*
     * Wilson coefficients Misiak basis
     * Operator block: L (2)
     * Normalization: 4 G_F / sqrt(2) x CKM
     */
    WilsonCoefficientNew& mc_L();

    /*
     * Wilson coefficients Misiak basis
     * Operator block: Q (4)
     * Normalization: 4 G_F / sqrt(2) x CKM
     * QED_NLO ref.: Gambino, Haisch, JHEP 0110, 020, hep-ph/0109058 - COULD BE CHANGED TO X,Y,W
     */
    WilsonCoefficientNew& mc_Q();

    /*
     * Wilson coefficients Misiak basis
     * Operator block: B (1)
     * Normalization: 4 G_F / sqrt(2) x CKM
     */
    WilsonCoefficientNew& mc_B();

    friend double gslpp_special_functions::dilog(double x);
    friend double gslpp_special_functions::clausen(double x);
    friend double gslpp_special_functions::zeta(int i);
    
};

#endif	/* STANDARDMODELMATCHING_H */
