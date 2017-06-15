/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EVOLDF1_H
#define EVOLDF1_H

#include <string>
#include "RGEvolutor.h"
#include "StandardModel.h"
#include "gslpp_special_functions.h"
#include <map>
#include "boost/multi_array.hpp"


class EvolDF1 : public RGEvolutor {
    /**
     * @class EvolDC1Buras
     * @brief \f$ |\Delta F = 1 | \f$ Evolutor Class
     * @details This evolutor is properly written to study \f$ |\Delta F = 1 | \f$ 
     * processes such as radiative and semileptonic weak decays of the B meson; 
     * it is implemented for the evolution of the 10 Wilson coefficients (related to the 
     * 2 current x current + 4 QCD penguins + em & chromo penguins + 2 semilptonic ones)
     * in the Chetyrkin, Misiak and Munz basis at the NLO in \f$ \alpha_{strong} \f$;
     * principal reference: hep-ph/0612329 
     */
public:    
    /**
     * @brief EvolDF1bsg constructor
     * @param dim an unsigned integer  for the dimension of the evolutor 
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @param order an enum "orders" for the order of perturbation theory of the evolutor
     * @param model an object of StandardModel class
     */
    EvolDF1(unsigned int nops, std::string reqblocks, schemes scheme, orders order, const StandardModel& model);
    /**
     * @brief EvolDF1bsg destructor
     */
    virtual ~EvolDF1();
    /**
     * @brief a method returning the anomalous dimension matrix given in the Misiak basis - qcd corrections
     * @param order an enum "orders" for the order of perturbation theory of the ADM
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM at the order LO/NLO in the Misiak basis - qcd corrections
     */
    gslpp::matrix<double> AnomalousDimension_s(orders order, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief a method returning the anomalous dimension matrix given in the Misiak basis - qed corrections
     * @param order an enum "orders" for the order of perturbation theory of the ADM
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM at the order LO/NLO in the Misiak basis - qed corrections
     */
    gslpp::matrix<double> AnomalousDimension_e(orders order, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief a method returning the evolutor related to the high scale \f$ M \f$ and the low scale \f$ \mu \f$
     * @param mu a double for the low scale of the evolution
     * @param M a double for the high scale of the evolution
     * @param order an enum "orders" for the order of perturbation theory of the evolutor
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @return the evolutor \f$ U (\mu , M) \f$
     */
    gslpp::matrix<double>& DF1Evol(double mu, double M, orders order, schemes scheme = NDR);
    /**
     * @brief a method returning the anomalous dimension in the Chetyrkin, Misiak and Munz operator basis 
     * @param order an enum "orders" for the order of perturbation theory of the evolutor
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM at the order LO/NLO in the Chetyrkin, Misiak and Munz basis
     */
    //gslpp::matrix<double> ToRescaleBasis(orders order, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief a method returning the anomalous dimension for the evolution of the effective Wilson coefficients
     * @param mat a temporary variable of gslpp::matrix type
     * @return the ADM at the order LO/NLO for the effective Wilson coefficients
     */
    //gslpp::matrix<double> ToEffectiveBasis(gslpp::matrix<double> mat)const;

    std::map<std::string,unsigned int> blocks_nops;

//    {{"C",2},{"CP",6},{"CPM",8},{"L",2},{"CPML",10},{"CPQB",11},{"CPMQB",13},{"CPMLQB",15}};
    std::map<std::string,orders> blocks_ord;// = {{"C",NNLO},{"CP",NNLO},{"CPM",NNLO},{"L",NNLO},{"CPML",NNLO},{"CPQB",NLO},{"CPMQB",NLO},{"CPMLQB",NLO}};

private:
   
    /**
     * @param nu an unsigned integer for the up-type number of d.o.f.
     * @param nu an unsigned integer for the down-type number of d.o.f.
     */
    int nu, nd;
    /**
     * @param a array of double for the magic numbers of the evolutor ( LO evolution )
     * @param b array of double for the magic numbers of the evolutor ( LO evolution ) 
     * @param c array of double for the magic numbers of the evolutor ( NLO evolution, associated to \f$ \alpha_{strong}(\mu) \f$ )
     * @param d array of double for the magic numbers of the evolutor ( NLO evolution, associated to \f$ \alpha_{strong}(M) \f$ )
     */
    //double a[4][13],
    //double b[4][13][13][13], c[4][13][13][13], d[4][13][13][13];

    typedef boost::multi_array<double, 4> array_type4;
    typedef boost::multi_array<double, 2> array_type2;

    array_type2 a;
    array_type4 b,c,d;

    const StandardModel& model;
    /**
     * @brief a void type method storing properly the magic numbers for the implementation of the evolutor
     * @param mu a double for the low scale of the evolution
     * @param M a double for the high scale of the evolution
     * @param nf a double for the active number of flavors
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     */
    void DF1Evol(double mu, double M, double nf, schemes scheme);
    
    friend double gslpp_special_functions::zeta(int i);

    gslpp::matrix<double> GammaCC_s(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaCP_s(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaCM_s(orders order, unsigned int n_u, unsigned int n_d) const;

    gslpp::matrix<double> GammaPP_s(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaPM_s(orders order, unsigned int n_u, unsigned int n_d) const;

    gslpp::matrix<double> GammaMM_s(orders order, unsigned int n_u, unsigned int n_d) const;

    gslpp::matrix<double> GammaQP_s(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaQM_s(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaQQ_s(orders order, unsigned int n_u, unsigned int n_d) const;

    gslpp::matrix<double> GammaBP_s(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaBB_s(orders order, unsigned int n_u, unsigned int n_d) const;

    gslpp::matrix<double> GammaCC_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaCP_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaCM_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaCL_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaCQ_e(orders order, unsigned int n_u, unsigned int n_d) const;

    gslpp::matrix<double> GammaPP_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaPM_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaPL_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaPQ_e(orders order, unsigned int n_u, unsigned int n_d) const;

    gslpp::matrix<double> GammaMM_e(orders order, unsigned int n_u, unsigned int n_d) const;
    
//    gslpp::matrix<double> GammaLP_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaLL_e(orders order, unsigned int n_u, unsigned int n_d) const;
//    gslpp::matrix<double> GammaLQ_e(orders order, unsigned int n_u, unsigned int n_d) const;

    gslpp::matrix<double> GammaQP_e(orders order, unsigned int n_u, unsigned int n_d) const;
//    gslpp::matrix<double> GammaQM_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaQL_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaQQ_e(orders order, unsigned int n_u, unsigned int n_d) const;

    gslpp::matrix<double> GammaBP_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaBL_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaBQ_e(orders order, unsigned int n_u, unsigned int n_d) const;
    gslpp::matrix<double> GammaBB_e(orders order, unsigned int n_u, unsigned int n_d) const;

    /**
     * @brief QCD beta function coefficients - eq. (36) hep-ph/0512066
     * @param powers of alpha_s and alpha_e as an integer
     * @param number of active flavor
     * @return coefficient of the beta function
     */
    double Beta_s(int i, double nf);

    /**
     * @brief QED beta function coefficients - eq. (36) hep-ph/0512066
     * @param powers of alpha_s and alpha_e as an integer
     * @param number of active flavor
     * @return coefficient of the beta function
     */
    double Beta_e(int i, double nf);
  
    // operators to include         {C, P, M, L, Q, b}            Huber et al., hep-ph/0512066
    unsigned int nops;
    std::string blocks;

    gslpp::matrix<gslpp::complex> v, vi, js, h, gg, s_s, jssv, jss, jv, vij;
    gslpp::vector<gslpp::complex> e;
    double alsMZ_cache;
    double Mz_cache;
};

#endif /* EVOLDF1_H */

