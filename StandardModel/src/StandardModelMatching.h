/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef STANDARDMODELMATCHING_H
#define	STANDARDMODELMATCHING_H

#include "ModelMatching.h"
#include <gslpp.h>

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
    
    /**
     *
     * @brief Updates to new Standard Model parameter sets.
     */
    
    void updateSMParameters();
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{d} \f$ 
     * @return return the vector of SM Wilson coefficients
     */
    virtual  std::vector<WilsonCoefficient>& CMdbd2();
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{s} \f$ 
     * @return return the vector of SM Wilson coefficients
     */
    virtual  std::vector<WilsonCoefficient>& CMdbs2();
    
    /**
     * 
     * @brief \f$ \Delta C = 2 \f$,
     * @return return the vector of SM Wilson coefficients
     */
    virtual  std::vector<WilsonCoefficient>& CMdd2();
    
    /**
     * 
     * @brief \f$ \Delta S = 2 \f$ 
     * @return return the vector of SM Wilson coefficients
     */
    virtual  std::vector<WilsonCoefficient>& CMdk2();
        
    /** 
     * 
     * @brief operator basis: 
     * @brief current current; qcd penguins; magnetic and chromomagnetic penguins; semileptonic  
     * @param[in] a a=0 deltaS=0 deltaC=0;  a=1 deltaS=1 deltaC=0;
     * @return Wilson coefficients Buras base for non-leptonic B decays 
     */
    virtual  std::vector<WilsonCoefficient>& CMbnlep( int a);
    
    /**
     * 
     * @brief operator basis: - current current opertors  
     * @param[in] a a=0 deltaS=0 deltaC=0;  a=1 1,0 ;  a=2 0,1 ; a=3 1,1
     * @return Wilson coefficients, Buras basis, for non-leptonic B decays 
     */
    virtual  std::vector<WilsonCoefficient>& CMbnlepCC( int a);
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l{-} \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbsg();
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l{-} \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMprimebsg();
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow K^* l^{+} l{-} \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMBMll();
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow K^* l^{+} l{-} \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMprimeBMll();
    
    /** 
     * 
     * @brief operator basis: 
     * @brief current-current; qcd penguins;
     * @brief hep/ph 9512380  
     * @return Wilson coefficients Buras base for \f$ K \rightarrow \pi \pi \f$ decays
     */
    virtual  std::vector<WilsonCoefficient>& CMK();
    
    /** 
     * 
     * @brief operator basis: 
     * @brief current-current (open up - open charm) 
     * @brief hep/ph 9512380 
     * @return Wilson coefficients Buras base for \f$ B \rightarrow \pi \pi \f$ decays
     */
    virtual  std::vector<WilsonCoefficient>& CMKCC();
    
    /**
     * 
     * @brief current-current oerators, Misiak basis
     * @return Wilson coefficients for \f$ D^{0} \rightarrow \pi \pi , K K \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMd1();
    
    /**
     * 
     * @brief current-current oerators, Buras basis
     * @return Wilson coefficients for \f$ D^{0} \rightarrow \pi \pi , K K \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMd1Buras();
    
    /**
     * 
     * @return Wilson coefficients for \f$ K_{L} \rightarrow \pi \nu \nu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMkpnn();
    /**
     * 
     * @return Wilson coefficient for \f$ K^{+} \rightarrow \mu \mu \f$, short distance top contribution
     */
    virtual  std::vector<WilsonCoefficient>& CMkmm();
    
    /**
     * 
     * @return Wilson coefficient for \f$ B_{s} \rightarrow \mu \mu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbsmm();
    
    /**
     * 
     * @return Wilson coefficient for \f$ B_{d} \rightarrow \mu \mu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbdmm();
    
    /**
     * 
     * @return Wilson coefficient for \f$ B \rightarrow \tau \nu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbtaunu();
    
    /**
     * 
     * @return Wilson coefficients for \f$ B_{s} \rightarrow X_{s} \nu \nu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMBXsnn();
    
    /**
     * 
     * @return Wilson coefficients for \f$ B_{d} \rightarrow X_{d} \nu \nu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMBXdnn();
    
     /**
     * 
     * @return Wilson coefficients for \f$ \ell_i \to \ell_j\f$
     */
    virtual  std::vector<WilsonCoefficient>& CMDLij(int li_lj);

     /**
     * 
     * @return Wilson coefficients for \f$ \ell_i \to \ell_j\ell_j\ell_j\f$
     */
    virtual  std::vector<WilsonCoefficient>& CMDLi3j(int li_lj);

     /**
     * 
     * @return Wilson coefficients for \f$ \mu \to e\f$ conversion in nuclei
     */
    virtual  std::vector<WilsonCoefficient>& CMmueconv();

     /**
     * 
     * @return Wilson coefficients for \f$ (g-2)_\mu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMgminus2mu();

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
     * @brief loop functions for non-leptonic B decays, Buiras Basis
     * Buras et al, hep-ph/9512380v1
     * @param[in] x the square of the ratio between top mass and W mass
     */
    double E0b(double x)const;
    
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

    double S0(double, double) const;

    /**
     * double penguin contribution to Kaon mixing - double top contribution
     * @return 
     */
    gslpp::complex ZDPtt() const;
    
    /**
     * double penguin contribution to Kaon mixing - charm top contribution
     * @return 
     */
    gslpp::complex ZDPct() const;
    
    
    
protected:
    std::vector<WilsonCoefficient> vmcdb, vmcds, vmcd2, vmck2, vmck, vmckcc;
    std::vector<WilsonCoefficient> vmcbsg, vmcprimebsg, vmcBMll, vmcprimeBMll, vmcbnlep, vmcbnlepCC, vmcd1, vmcd1Buras;
    std::vector<WilsonCoefficient> vmckpnn, vmckmm, vmcbsnn, vmcbdnn, vmcbsmm, vmcbdmm, vmcbtaunu;
    std::vector<WilsonCoefficient> vmcDLij, vmcDLi3j, vmcmueconv, vmcgminus2mu;
    
private:
    
    const StandardModel & SM;
    double S0(double) const;
    double S0p(double x) const;
    double S11(double x) const;
    double S18(double x) const;
    double ZDP(const double x, const double y) const;
    WilsonCoefficient mcdbd2, mcdbs2, mcdd2, mcdk2, mck, mckcc;
    WilsonCoefficient mcbsg, mcprimebsg, mcBMll, mcprimeBMll, mcbnlep, mcbnlepCC, mcd1, mcd1Buras;
    WilsonCoefficient mckpnn, mckmm, mcbsnn, mcbdnn, mcbsmm, mcbdmm, mcbtaunu;
    WilsonCoefficient mcDLij, mcDLi3j, mcmueconv, mcgminus2mu;
    
    double Mut;
    double Muw;
    double Ale;
    double GF;
    double Mw_tree;
    double Nc;
    double CF;
    double gamma0;
    double J5;
    double BtNDR;
    double Mw;
    double sW2;
    double mu_b;
    //double MM;
    gslpp::matrix<gslpp::complex> Vckm;
    gslpp::complex lam_t;
    double L;
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order
     * @return return the value of the wilson coefficients for \f$ B_{s} \rightarrow  l^{+} l{-} \f$
     */
    double setWCBsmm(int i, double x, orders order);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order_ew
     * @return return the electroweak value of the wilson coefficients for \f$ B_{s} \rightarrow  l^{+} l{-} \f$
     */
    double setWCBsmmEW(int i, double x, orders_ew order_ew);
    
     /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order
     * @return return the value of the wilson coefficients for \f$ B_{d} \rightarrow  l^{+} l{-} \f$
     */
    double setWCBdmm(int i, double x, orders order);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order_ew
     * @return return the electroweak value of the wilson coefficients for \f$ B_{d} \rightarrow  l^{+} l{-} \f$
     */
    double setWCBdmmEW(int i, double x, orders_ew order_ew);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order
     * @return return the value of the wilson coefficients for \f$ B \rightarrow X_{s} \gamma, l^{+} l{-} \f$
     */
    double setWCbsg (int i, double x, orders order);
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and W mass
     * @param order
     * @return return the value of the wilson coefficients for \f$ B \rightarrow k^* l^{+} l{-} \f$
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
    
    double CWbsgArrayLO[8], CWbsgArrayNLO[8];
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
    
    
};

#endif	/* STANDARDMODELMATCHING_H */

