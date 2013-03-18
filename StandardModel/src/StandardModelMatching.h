/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef STANDARDMODELMATCHING_H
#define	STANDARDMODELMATCHING_H

//#include "StandardModel.h"
#include "ModelMatching.h"

#define LEPS 1.e-5 // tolerance in the limit of S(x,y) to S(x)

class StandardModel;

/**
 * @class StandardModelMatching
 * @ingroup StandardModel
 * @brief A class for the matching in the Standard Model. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class StandardModelMatching : public ModelMatching {
public:
    StandardModelMatching(const StandardModel & SM_i);
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{d} \f$ 
     * @return return the vector of SM Wilson coefficients
     */
    virtual const std::vector<WilsonCoefficient>& CMdbd2();
    
    /**
     * 
     * @brief \f$ \Delta B = 2 \f$, \f$ B_{s} \f$ 
     * @return return the vector of SM Wilson coefficients
     */
    virtual const std::vector<WilsonCoefficient>& CMdbs2();
    
    /**
     * 
     * @brief \f$ \Delta C = 2 \f$,
     * @return return the vector of SM Wilson coefficients
     */
    virtual const std::vector<WilsonCoefficient>& CMdd2();
    
    /**
     * 
     * @brief \f$ \Delta S = 2 \f$ 
     * @return return the vector of SM Wilson coefficients
     */
    virtual const std::vector<WilsonCoefficient>& CMdk2();
    
    /** 
     * 
     * @brief operator basis: 
     * @brief current current; qcd penguins; magnetic and chromomagnetic penguins; semileptonic  
     * @param a, a=0 deltaS=0 deltaC=0;  a=1 deltaS=1 deltaC=0;
     * @return Wilson coefficients Buras base for non-leptonic B decays 
     */
    virtual const std::vector<WilsonCoefficient>& CMbnlep(const int& a);
    
    /**
     * 
     * @brief operator basis: - current current opertors  
     * @param a, a=0 deltaS=0 deltaC=0;  a=1 1,0 ;  a=2 0,1 ; a=3 1,1 
     * @return Wilson coefficients, Buras basis, for non-leptonic B decays 
     */
    virtual const std::vector<WilsonCoefficient>& CMbnlepCC(const int& a);
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l{-} \f$
     */
    virtual const std::vector<WilsonCoefficient>& CMbsg();
    
    /** 
     * 
     * @brief operator basis: 
     * @brief current-current; qcd penguins;
     * @brief hep/ph 9512380  
     * @return Wilson coefficients Buras base for \f$ K \rightarrow \pi \pi \f decays
     */
    virtual const std::vector<WilsonCoefficient>& CMK();
    
    /** 
     * 
     * @brief operator basis: 
     * @brief current-current (open up - open charm) 
     * @brief hep/ph 9512380 
     * @return Wilson coefficients Buras base for \f$ B \rightarrow \pi \pi \f decays
     */
    virtual const std::vector<WilsonCoefficient>& CMKCC();
    
    /**
     * 
     * @brief current-current oerators, Misiak basis
     * @return Wilson coefficients for \f$ D^{0} \rightarrow \pi \pi , K K \f$
     */
    virtual const std::vector<WilsonCoefficient>& CMd1();
    
    /**
     * 
     * @brief current-current oerators, Buras basis
     * @return Wilson coefficients for \f$ D^{0} \rightarrow \pi \pi , K K \f$
     */
    virtual const std::vector<WilsonCoefficient>& CMd1Buras();
    
    /**
     * 
     * @return Wilson coefficients for \f$ K_{L} \rightarrow \pi \nu \nu \f$
     */
    virtual const std::vector<WilsonCoefficient>& CMkp0nn();
    
    /**
     * 
     * @return Wilson coefficients for \f$ K^{+} \rightarrow \pi^{+} \nu \nu \f$
     */
    virtual const std::vector<WilsonCoefficient>& CMkppnn();
    
    /**
     * 
     * @return Wilson coefficient for \f$ K^{+} \rightarrow \mu \mu \f$, short distance top contribution
     */
    virtual const std::vector<WilsonCoefficient>& CMkmm();
    
    /**
     * 
     * @return Wilson coefficient for \f$ B_{s} \rightarrow \mu \mu \f$
     */
    virtual const std::vector<WilsonCoefficient>& CMbsmm();
    
    /**
     * 
     * @return Wilson coefficient for \f$ B_{d} \rightarrow \mu \mu \f$
     */
    virtual const std::vector<WilsonCoefficient>& CMbdmm();
    
    /**
     * 
     * @return Wilson coefficients for \f$ B_{s} \rightarrow X_{s} \nu \nu \f$
     */
    virtual const std::vector<WilsonCoefficient>& CMBXsnn();
    
    /**
     * 
     * @return Wilson coefficients for \f$ B_{d} \rightarrow X_{s} \nu \nu \f$
     */
    virtual const std::vector<WilsonCoefficient>& CMBXdnn();

    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param x the square ratio between top mass and W mass
     */
    double A0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the semileptonic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param x the square ratio between top mass and W mass
     */
    double B0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param x the square ratio between top mass and W mass
     */
    double C0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param x the square ratio between top mass and W mass
     */
    double D0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the chromomagnetic operator
     * in the non-effective Misiak basis, Bobeth et al hep-ph/9910220 
     * @param x the square ratio between top mass and W mass
     */
    double F0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the chromomagnetic operator
     * in the Misiak basis, Chetyrkin et al hep-ph/9612313 
     * @param x the square ratio between top mass and W mass
     */
    double E0t(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the effective Misiak basis, LO term, Chetyrkin et al hep-ph/9612313 
     * @param x the square ratio between top mass and W mass
     */
    double C7LOeff(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the magnetic operator
     * in the effective Misiak basis, NLO term, Chetyrkin et al hep-ph/9612313 
     * @param x the square ratio between top mass and W mass
     */
    double C7NLOeff(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the chromomagnetic operator
     * in the effective Misiak basis, LO term, Chetyrkin et al hep-ph/9612313 
     * @param x the square ratio between top mass and W mass
     */
    double C8LOeff(double x)const;
    
    /**
     * 
     * @brief loop function which appear in the Wilson coefficient for the chromomagnetic operator
     * in the effective Misiak basis, LO term, Chetyrkin et al hep-ph/9612313 
     * @param x the square ratio between top mass and W mass
     */
    double C8NLOeff(double x)const;

    /**
     * 
     * @brief loop functions for non-leptonic B decays, Buiras Basis
     * Buras et al, hep-ph/9512380v1
     * @param x the square ratio between top mass and W mass
     */
    double B0b(double x)const;
    
    /**
     * 
     * @brief loop functions for non-leptonic B decays, Buiras Basis
     * Buras et al, hep-ph/9512380v1
     * @param x the square ratio between top mass and W mass
     */
    double C0b(double x)const;
    
    /**
     * 
     * @brief loop functions for non-leptonic B decays, Buiras Basis
     * Buras et al, hep-ph/9512380v1
     * @param x the square ratio between top mass and W mass
     */
    double D0b(double x)const;
    
    /**
     * 
     * @brief loop functions for non-leptonic B decays, Buiras Basis
     * Buras et al, hep-ph/9512380v1
     * @param x the square ratio between top mass and W mass
     */
    double E0b(double x)const;
    
    /**
     *  
     * @brief hep-ph/9512380
     * @return the loop function for the charm-charm contribution to the Delta S = 2 effective hamiltonian multiplied by the CKM element 
     */
    complex S0c() const;
    
    /**
     *  
     * @brief hep-ph/9512380
     * @return the loop function for the charm-top contribution to the Delta S = 2 effective hamiltonian multiplied by the CKM element 
     */
    complex S0ct() const;
    
    /**
     *  
     * @brief hep-ph/9512380
     * @return the loop function for the top-top contribution to the Delta S = 2 effective hamiltonian
     */
    complex S0tt() const;
    
    /**
     *
     * @brief hep-ph/1009.0947v2
     * @return 0th order loop function for the top contribution to K -> pi nu nu decays
     */
    double X0t(double x)const;
    
    /**
     *
     * @brief hep-ph/1009.0947v2
     * @return first order in QCD loop function for the top contribution to K -> pi nu nu decays
     */
    double X1t(double x)const;
    
    /**
     *
     * @brief hep-ph/1009.0947v2
     * @return two loop EW loop function for the top contribution to K -> pi nu nu decays
     * @return in the limit of small theta_W
     */
    double Xewt(double x, double a, double mu)const;
    
    /**
     * 
     * @param x
     * @return 0th order loop function for the top contribution to K_L -> mu mu decays
     */
    double Y0t(double x)const;
    
    /**
     * 
     * @param x
     * @return first order in QCD loop function for the top contribution to K_L -> mu mu decays
     */
    double Y1t(double x)const;
    
protected:
    std::vector<WilsonCoefficient> vmcdb, vmcds, vmcd2, vmck2, vmck, vmckcc;
    std::vector<WilsonCoefficient> vmcbsg, vmcbnlep, vmcbnlepCC, vmcd1, vmcd1Buras;
    std::vector<WilsonCoefficient> vmckpnn, vmckppnn, vmckmm, vmcbsnn, vmcbdnn, 
                                   vmcbsmm, vmcbdmm;
    
    
private:
    const StandardModel & SM;
    double S0(double, double) const;
    double S0(double) const;
    double S0p(double x) const;
    double S11(double x) const;
    double S18(double x) const;
    double S1(double x) const;
    WilsonCoefficient mcdbd2, mcdbs2, mcdd2, mcdk2, mck, mckcc;  
    WilsonCoefficient mcbsg, mcbnlep, mcbnlepCC, mcd1, mcd1Buras;
    WilsonCoefficient mckpnn, mckppnn, mckmm, mcbsnn, mcbdnn, mcbsmm, mcbdmm;
    
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
     * @param z
     * @return two loop EW loop functions for K-> P nu nu, hep-ph/1009.0947v2 
     */
    double phi2 (double x, double y) const;
    
    double CWbsgArrayLO[10], CWbsgArrayNLO[10];
    double CWD1ArrayLO[10], CWD1ArrayNLO[10];
    double CWbnlepArrayLOqcd[10], CWbnlepArrayNLOqcd[10];
    double CWbnlepArrayLOew[10], CWbnlepArrayNLOew[10];
    
    double sw, swa, swb, swc; //sen(theta_W) tree level
    double xcachea, xcacheb, xcachec; // caching
};

#endif	/* STANDARDMODELMATCHING_H */


