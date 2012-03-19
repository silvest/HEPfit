/* 
 * File:   StandardModelMatching.h
 * Author: silvest
 *
 * Created on June 9, 2011, 2:16 PM
 */

#ifndef STANDARDMODELMATCHING_H
#define	STANDARDMODELMATCHING_H

#include "StandardModel.h"
#include "ModelMatching.h"

#define LEPS 1.e-5 // tolerance in the limit of S(x,y) to S(x)

class StandardModelMatching : public ModelMatching{
public:
    StandardModelMatching(const StandardModel& SM_i);
    
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
    
protected:
    std::vector<WilsonCoefficient> vmcdb, vmcds, vmcd2, vmck2;
    std::vector<WilsonCoefficient> vmcbsg, vmcbnlep, vmcbnlepCC, vmcd1, vmcd1Buras;
    
    
private:
    const StandardModel& SM;
    double S0(double, double) const;
    double S0(double) const;
    double S0p(double x) const;
    double S11(double x) const;
    double S18(double x) const;
    double S1(double x) const;
    WilsonCoefficient mcdbd2, mcdbs2, mcdd2, mcdk2;  
    WilsonCoefficient mcbsg, mcbnlep, mcbnlepCC, mcd1, mcd1Buras;
    
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
    
    double CWbsgArrayLO[10], CWbsgArrayNLO[10];
    double CWD1ArrayLO[10], CWD1ArrayNLO[10];
    double CWbnlepArrayLOqcd[10], CWbnlepArrayNLOqcd[10];
    double CWbnlepArrayLOew[10], CWbnlepArrayNLOew[10];
    
    double sw, swa, swb, swc; //sen(theta_W) tree level
    double xcachea, xcacheb, xcachec; // caching
};

#endif	/* STANDARDMODELMATCHING_H */


