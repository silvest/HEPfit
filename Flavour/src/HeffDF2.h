/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDF2_H
#define	HEFFDF2_H

#include <StandardModel.h>
#include <SUSYMassInsertion.h>
#include <SUSYMassInsertionMatching.h>
#include <WilsonCoefficient.h>
#include "EvolDF2.h"

using namespace gslpp;

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
     * @brief change scheme from NDR to LRI
     * @param schout is the renormalization scheme in output
     * @param schin is the renormalization scheme in output
     * @param order
     * @param j, j=0 B_s, j=1 B_d, j=2 K, j=3 D.
     */ 
    void ChangeScheme(schemes schout, schemes schin, orders order, int j);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B_d oscillations
     */
    vector<complex>** ComputeCoeffBd(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B_s oscillations
     */
    vector<complex>** ComputeCoeffBs(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for D oscillations
     */
    vector<complex>** ComputeCoeffdd(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for K oscillations
     */
    vector<complex>** ComputeCoeffK(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @brief for Delta M_K the SM contribution is set to zero
     * @return the effective hamiltonian at the scale mu for Delta M_K
     */
    vector<complex>** ComputeCoeffmK(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param order
     * @param nf is the number of active flavours
     * @return the anomalous dimension for DeltaF=2 processes
     */
    matrix<double> AnomalousDimension(orders order, unsigned int nf = 0) const;

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
    complex S0tt(double mu) const;
    const StandardModel& model;
    matrix<double> drNDRLRI;
    WilsonCoefficient coeffbd, coeffbs, coeffDd, coeffk, coeffmk;
    
    EvolDF2 evolDF2;
};

#endif	/* HEFFDF2_H */

