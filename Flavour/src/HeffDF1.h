/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDF1_H
#define HEFFDF1_H

#include "StandardModel.h"
#include "StandardModelMatching.h"
#include "WilsonCoefficientNew.h"
#include "EvolDF1.h"
#include <map>

#define N_OPS 15 // number of operators in the basis

class HeffDF1 {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDF1(std::string blocks, const StandardModel & SM, qcd_orders order_qcd = QCD1, qed_orders order_qed = QED0);

    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDF1() {
    };

    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> K^*ll decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    virtual Expanded<gslpp::vector<gslpp::complex> > ComputeCoeff(double mu, schemes scheme = NDR);

    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> K^*ll decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    //gslpp::vector<gslpp::complex>** ComputeCoeffprime(double mu, schemes scheme = NDR);

    /**
     * 
     * @param Coeff vector of Wilson coefficient
     * @param nm order of the expansion
     * @return the coefficient of the expansion in low-energy coupling constants as defined in eq. (68) of Huber et al., hep-ph/0512066
     */
    gslpp::vector<gslpp::complex> LowScaleCoeff(qcd_orders order_qcd, qed_orders order_qed);

    EvolDF1 getEvol() const {
        return evolDF1;
    }

    const StandardModel& GetModel() const {
        return model;
    }
    
protected:
    WilsonCoefficientNew coeff;

private:
    const StandardModel& model;

    EvolDF1 evolDF1;

    std::string blocks;
    unsigned int nops;
    double mu_cache;
    schemes scheme_cache;
    std::vector<double> Vmu_cache;
    std::vector<WilsonCoefficientNew> WC_cache;
};


class HeffDF1_NP : public HeffDF1 {
public:
    /**
     * @brief constructor
     * @param SM
     * @param modelmatching
     */
    HeffDF1_NP(std::string blocks, const StandardModel & SM, qcd_orders order_qcd, qed_orders order_qed);

    /**
     * 
     * @brief destructor
     */
    virtual ~HeffDF1_NP() {
    };
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> K^*ll decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    Expanded<gslpp::vector<gslpp::complex> > ComputeCoeff(double mu, schemes scheme = NDR);
    
private:
    gslpp::vector<gslpp::complex> coeffNP;
};

#endif /* HEFFDF1_H */
