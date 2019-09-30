/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef HEFFDB1_H
#define	HEFFDB1_H

class StandardModel;
class StandardModelMatching;
class EvolDB1Mll;
class EvolDB1bsg;
class EvolBsmm;
#include "QCD.h"
#include "WilsonCoefficient.h"
#include <memory>

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
    gslpp::vector<gslpp::complex>** ComputeCoeffBnlep00(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 1 | \f$, \f$ |\Delta S = 0 | \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBnlep10(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 0 | \f$, \f$ |\Delta S = 1 | \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBnlep01(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu for B decays, \f$ |\Delta C = 1 | \f$, \f$ |\Delta S = 1 | \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBnlep11(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B_{s} \rightarrow \mu \bar{\mu} \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffsmumu(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B_{s} \rightarrow \mu \bar{\mu} \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffsmumuStandardNorm(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B_{d} \rightarrow \mu \bar{\mu} \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffdmumu(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B \rightarrow \tau \nu \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffbtaunu(QCD::meson meson_i);
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B_{s} \rightarrow \nu \bar{\nu} \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffsnunu();
    
    /**
     * 
     * @param scheme
     * @return short distance contribution to the rare decay \f$ B_{d} \rightarrow \nu \bar{\nu} \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffdnunu();
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return short distance contribution to the rare decay \f$ b \rightarrow s \gamma \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffsgamma(double mu, bool noSM = false, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return short distance contribution to the rare decay \f$ b \rightarrow s \gamma \f$
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffprimesgamma(double mu, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> K^*ll decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffBMll(double mu, QCD::lepton lepton, bool noSM = false, schemes scheme = NDR);
    
    /**
     * 
     * @param mu is the low energy scale
     * @param scheme indicates the renormalization scheme
     * @return the effective hamiltonian at the scale mu B -> K^*ll decay, Misiak basis, Chetyrkin et al hep-ph/9612313
     */
    gslpp::vector<gslpp::complex>** ComputeCoeffprimeBMll(double mu, QCD::lepton lepton, schemes scheme = NDR);
    
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
    
    WilsonCoefficient getCoeffbtaunu() const {
        return coeffbtaunu;
    }
    
    WilsonCoefficient getCoeffsnunu() const {
        return coeffsnunu;
    }
    
    WilsonCoefficient getCoeffdnunu() const {
        return coeffdnunu;
    }
    
    WilsonCoefficient getCoeffsgamma() const {
        return coeffsgamma;
    }
    
    WilsonCoefficient getCoeffprimesgamma() const {
        return coeffprimesgamma;
    }
    
    EvolBsmm& getUBsmm() const {
        return *evolbs;
    }
    
    EvolBsmm& getUBdmm() const {
        return *evolbd;
    }
    
    EvolDB1bsg& getUDB1bsg() const {
        return *evolDB1bsg;
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
    WilsonCoefficient coeffbtaunu;
    WilsonCoefficient coeffsnunu, coeffdnunu;
    WilsonCoefficient coeffsgamma, coeffprimesgamma;
    WilsonCoefficient coeffBMll, coeffprimeBMll;
    std::unique_ptr<EvolDB1Mll> evolDF1BMll;
    std::unique_ptr<EvolDB1bsg> evolDB1bsg;
    std::unique_ptr<EvolBsmm> evolbs;
    std::unique_ptr<EvolBsmm> evolbd;
    
    
    double Bsgamma_mu_cache;
    unsigned int Bsgamma_order_ini_cache;
    schemes Bsgamma_scheme_cache;
    std::vector<double> Bsgamma_Mu_cache;
    std::vector<double> Bpsgamma_Mu_cache;
    std::vector<WilsonCoefficient> Bsgamma_WC_cache;
    std::vector<WilsonCoefficient> Bpsgamma_WC_cache;
    
    double BMll_mu_cache;
    unsigned int BMll_order_ini_cache;
    schemes BMll_scheme_cache;
    std::vector<double> BMll_Mu_cache;
    std::vector<WilsonCoefficient> BMll_WC_cache;
    
    double BMllprime_mu_cache;
    schemes BMllprime_scheme_cache;
    std::vector<double> BMllprime_Mu_cache;
    std::vector<WilsonCoefficient> BMllprime_WC_cache;
    
    double Bsmumu_mu_cache;
    schemes Bsmumu_scheme_cache;
    std::vector<double> Bsmumu_Mu_cache;
    std::vector<WilsonCoefficient> Bsmumu_WC_cache;
    
    double Bdmumu_mu_cache;
    schemes Bdmumu_scheme_cache;
    std::vector<double> Bdmumu_Mu_cache;
    std::vector<WilsonCoefficient> Bdmumu_WC_cache;
    
    gslpp::vector<gslpp::complex> nlep, nlep2, nlepCC;
};

#endif	/* HEFFDB1_H */
