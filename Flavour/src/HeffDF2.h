/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDF2_H
#define	HEFFDF2_H

#include <StandardModel.h>
#include <WilsonCoefficient.h>
#include "EvolDF2.h"

/**
 * @class HeffDF2
 * @ingroup Flavour
 * @brief A class for the @f$\Delta F = 2@f$ effective Hamiltonian.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class is used to compute the contributions to the
 * @f$\Delta F = 2@f$ by taking the values of the Wilson coefficients
 * from the model matching classes and runs them down to the relevant scale
 * for the observables.
 */
class HeffDF2 {
public:
    /**
     * @brief constructor
     * @param SM 
     * @param SM_Matching 
     */
    HeffDF2(const StandardModel& SM);
    
    /**
     * @brief destructor
     */
    virtual ~HeffDF2();
    
    /**
     * @brief change scheme for a Wilson Coefficient
     * @param schout is the renormalization scheme in output
     * @param c_in is the Wilson Coefficient to be converted to scheme schout
     * @param order
     */ 
    void ChangeScheme(schemes schout, WilsonCoefficient& c_in, orders order);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme (Default: NDR)
     * @return the effective hamiltonian at the scale mu for B_d oscillations
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBd(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme (Default: NDR)
     * @return the effective hamiltonian at the scale mu for B_s oscillations
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBs(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme (Default: NDR)
     * @return the effective hamiltonian at the scale mu for D oscillations
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffdd(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme (Default: NDR)
     * @return the effective hamiltonian at the scale mu for K oscillations
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffK(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme (Default: NDR)
     * @brief for Delta M_K the SM contribution is set to zero
     * @return the effective hamiltonian at the scale mu for Delta M_K
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffmK(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param order
     * @param nf is the number of active flavours
     * @return the anomalous dimension for DeltaF=2 processes
     */
    gslpp::matrix<double> AnomalousDimension(orders order, unsigned int nf = 0) const;

    WilsonCoefficient getCoeffBd() const {
        return coeffbd;
    }

    WilsonCoefficient getCoeffBs() const {
        return coeffbs;
    }
    
    WilsonCoefficient getCoeffDD() const {
        return coeffDd;
    }
    
    WilsonCoefficient getCoeffK() const {
        return coeffk;
    }
    
     WilsonCoefficient getCoeffmK() const {
        return coeffmk;
    }
    
    EvolDF2 getUDF2() const {
        return evolDF2;
    }


private:
    gslpp::complex S0tt(double mu) const;
    const StandardModel& model;
    gslpp::matrix<double> drNDRLRI;
    WilsonCoefficient coeffbd;
    WilsonCoefficient coeffbs;
    WilsonCoefficient coeffDd;
    WilsonCoefficient coeffk;
    WilsonCoefficient coeffmk;
    
    EvolDF2 evolDF2;
};

#endif	/* HEFFDF2_H */

