/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDB1_H
#define	HEFFDB1_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>
#include "EvolDF1nlep.h"
#include "EvolDF1bsg.h"
#include <sstream>

using namespace gslpp;

class HeffDB1 {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDB1(const StandardModel & SM);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDB1();
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the reonrmalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 0 | \f$, \f$ |\Delta S = 0 | \f$
     */
    vector<complex>** ComputeCoeffBnlep00(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 1 | \f$, \f$ |\Delta S = 0 | \f$
     */
    vector<complex>** ComputeCoeffBnlep10(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 0 | \f$, \f$ |\Delta S = 1 | \f$
     */
    vector<complex>** ComputeCoeffBnlep01(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 1 | \f$, \f$ |\Delta S = 1 | \f$
     */
    vector<complex>** ComputeCoeffBnlep11(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B_{s} \rightarrow \mu \bar{\mu} \f$
     */
    vector<complex>** ComputeCoeffsmumu();
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B_{d} \rightarrow \mu \bar{\mu} \f$
     */
    vector<complex>** ComputeCoeffdmumu();
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B_{s} \rightarrow \nu \bar{\nu} \f$
     */
    vector<complex>** ComputeCoeffsnunu();
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B_{d} \rightarrow \nu \bar{\nu} \f$
     */
    vector<complex>** ComputeCoeffdnunu();
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> K^*ll decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    vector<complex>** ComputeCoeffBKstarll(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> K^*ll decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    vector<complex>** ComputeCoeffprimeBKstarll(double mu, schemes scheme = NDR);
    
    WilsonCoefficient getCoeffnlep00() const {
        return coeffnlep00;
    }
    
    WilsonCoefficient getCoeffnlep10() const {
        return coeffnlep01;
    }
    
    WilsonCoefficient getCoeffnlep01() const {
        return coeffnlep10;
    }
    
    WilsonCoefficient getCoeffnlep11() const {
        return coeffnlep11;
    }
    
    WilsonCoefficient getCoeffsmumu() const {
        return coeffsmumu;
    }
    
    WilsonCoefficient getCoeffdmumu() const {
        return coeffdmumu;
    }
    
    WilsonCoefficient getCoeffsnunu() const {
        return coeffsnunu;
    }
    
    WilsonCoefficient getCoeffdnunu() const {
        return coeffdnunu;
    }
    
    EvolDF1nlep getUDF1() const {
        return u;
    }
    
    EvolDF1bsg getUDF1BKstarll() const {
        return evolDF1BKstarll;
    }

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;
    
    WilsonCoefficient coeffnlep00qcd, coeffnlep00;
    WilsonCoefficient coeffnlep10qcd, coeffnlep10;
    WilsonCoefficient coeffnlep01, coeffnlep01A, coeffnlep01B, coeffnlep00CC;
    WilsonCoefficient coeffnlep11, coeffnlep11A, coeffnlep11B, coeffnlep10CC;
    WilsonCoefficient coeffsmumu, coeffdmumu;
    WilsonCoefficient coeffsnunu, coeffdnunu;
    WilsonCoefficient coeffBKstarll, coeffprimeBKstarll;
    EvolDF1bsg evolDF1BKstarll;
    EvolDF1nlep u;
    
    //StandardModelMatching& standardmodelmatching;

    
    gslpp::vector<complex> nlep, nlep2, nlepCC;
};

#endif	/* HEFFDB1_H */
