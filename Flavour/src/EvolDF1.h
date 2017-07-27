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

typedef unsigned int indices;

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
     * @brief EvolDF1 constructor
     * @param dim an unsigned integer  for the dimension of the evolutor 
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @param order an enum "orders" for the order \f$ \alpha_s\f$ in the evolutor
     * @param order_qed an enum "orders_qed" for the order \f$ \alpha_e\f$ in the evolutor
     * @param model an object of StandardModel class
     */
    EvolDF1(unsigned int nops, std::string reqblocks, schemes scheme, const StandardModel& model, orders order, orders_qed order_qed);
    /**
     * @brief EvolDF1 destructor
     */
    virtual ~EvolDF1();
    /**
     * @brief a method returning the anomalous dimension matrix given in the Misiak basis
     * @param nm indices nm corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM in the Misiak basis
     */
    gslpp::matrix<double> AnomalousDimension(indices nm, unsigned int n_u, unsigned int n_d) const;

    /**
     * @brief a method returning the evolutor related to the high scale \f$ M \f$ and the low scale \f$ \mu \f$
     * @param mu a double for the low scale of the evolution
     * @param M a double for the high scale of the evolution
     * @param order an enum "orders" for the order of perturbation theory of the evolutor
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @return the evolutor \f$ U (\mu , M) \f$
     */
    gslpp::matrix<double>& DF1Evol(double mu, double M, orders order, orders_qed order_qed, schemes scheme = NDR);
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

    double alphatilde_s(double mu) { return 0.; } // to be fixed
        
private:
   
    /**
     * @param a array of double for the magic numbers of the evolutor ( LO evolution )
     * @param b array of double for the magic numbers of the evolutor ( LO evolution ) 
     * @param c array of double for the magic numbers of the evolutor ( NLO evolution, associated to \f$ \alpha_{strong}(\mu) \f$ )
     * @param d array of double for the magic numbers of the evolutor ( NLO evolution, associated to \f$ \alpha_{strong}(M) \f$ )
     */
    //double a[4][13],
    //double b[4][13][13][13], c[4][13][13][13], d[4][13][13][13];

//    typedef boost::multi_array<double, 4> array_type4;
//    typedef boost::multi_array<double, 2> array_type2;
//    array_type2 mn_a;
//    array_type4 mn_b,mn_c,mn_d;
    std::map< std::vector<int>, double > ai;
    std::map< std::vector<int>, gslpp::complex > vM0vi, vM1vi, vM2vi, vM11vi, vM3vi, vM4vi, vM5vi, vM6vi, vM33vi, vM31vi, vM13vi,
    vM34vi, vM43vi, vM23vi, vM32vi, vM14vi, vM41vi, vM113vi, vM131vi, vM311vi, vM133vi, vM313vi, vM331vi;
    std::vector<int> index;

    const StandardModel& model;

    /**
     * @brief map indexes to a vector<int>
     * @param nf number of active flavours
     * @param i,j,... indices in eqs. (49)-(53) in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @return a reference to the vector index
     */    
    std::vector<int>& idx(int nf, int a, int b=0, int i=0, int j=0, int k=0, int l=0, int p=0, int m=0, int q=0);

    /**
     * @brief Check if anomalous dimension indices and Nf match
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param nf number of active flavours
     */    
    void CheckNf(indices nm, unsigned int nf) const;
    
    /**
     * @brief a void type method storing properly the magic numbers for the implementation of the evolutor
     * @param mu a double for the low scale of the evolution
     * @param M a double for the high scale of the evolution
     * @param nf a double for the active number of flavors
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     */
//    void DF1Evol(double mu, double M, double nf, schemes scheme);
    void DF1Evol(double mu, double M, int nf, schemes scheme);

    friend double gslpp_special_functions::zeta(int i);

    /**
     * @brief CC block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM CC block in the Misiak basis
     */    
    gslpp::matrix<double> GammaCC(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief CP block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM CP block in the Misiak basis
     */    
    gslpp::matrix<double> GammaCP(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief CM block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM CM block in the Misiak basis
     */    
    gslpp::matrix<double> GammaCM(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief CL block of the QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM CL block in the Misiak basis
     */    
    gslpp::matrix<double> GammaCL(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief CQ block of the QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM CQ block in the Misiak basis
     */    
    gslpp::matrix<double> GammaCQ(indices nm, unsigned int n_u, unsigned int n_d) const;

    /**
     * @brief PP block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM PP block in the Misiak basis
     */    
    gslpp::matrix<double> GammaPP(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief PM block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM PM block in the Misiak basis
     */    
    gslpp::matrix<double> GammaPM(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief PL block of the QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM PL block in the Misiak basis
     */    
    gslpp::matrix<double> GammaPL(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief PQ block of the QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM PQ block in the Misiak basis
     */    
    gslpp::matrix<double> GammaPQ(indices nm, unsigned int n_u, unsigned int n_d) const;

    /**
     * @brief MM block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM MM block in the Misiak basis
     */    
    gslpp::matrix<double> GammaMM(indices nm, unsigned int n_u, unsigned int n_d) const;

    /**
     * @brief LL block of the QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM LL block in the Misiak basis
     */    
    gslpp::matrix<double> GammaLL(indices nm, unsigned int n_u, unsigned int n_d) const;

    /**
     * @brief QP block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM QP block in the Misiak basis
     */    
    gslpp::matrix<double> GammaQP(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief QM block of the QCD anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM QM block in the Misiak basis
     */    
    gslpp::matrix<double> GammaQM(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief QL block of the QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM QL block in the Misiak basis
     */    
    gslpp::matrix<double> GammaQL(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief QQ block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM QQ block in the Misiak basis
     */    
    gslpp::matrix<double> GammaQQ(indices nm, unsigned int n_u, unsigned int n_d) const;

    /**
     * @brief BP block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM BP block in the Misiak basis
     */    
    gslpp::matrix<double> GammaBP(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief BL block of the QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM BL block in the Misiak basis
     */    
    gslpp::matrix<double> GammaBL(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief BQ block of the QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM BQ block in the Misiak basis
     */    
    gslpp::matrix<double> GammaBQ(indices nm, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief BB block of the QCD+QED anomalous dimension
     * @param nm indices corresponding to powers of alpha_s and alpha_em as in Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM BB block in the Misiak basis
     */    
    gslpp::matrix<double> GammaBB(indices nm, unsigned int n_u, unsigned int n_d) const;

    /**
     * @brief QCD beta function coefficients including QED corrections - eq. (36) hep-ph/0512066
     * @param powers of alpha_s and alpha_e as an integer
     * @param number of active flavor
     * @return coefficient of the beta function
     */
    double Beta_s(int i, unsigned int nf);
    
    /**
     * @brief QED beta function coefficients - eq. (36) hep-ph/0512066
     * @param powers of alpha_s and alpha_e as an integer
     * @param number of active flavor
     * @return coefficient of the beta function
     */
    double Beta_e(int i, unsigned int nf);


    /**
     * @brief auxiliary function f - eq. (50) of Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param i matrix index
     * @param j matrix index
     * @param k order index
     * @param eta als(M)/als(mu)
     * @return function value
     */
    double f_f(unsigned int nf, unsigned int i, unsigned int j, int k, double eta);

    /**
     * @brief auxiliary function r - eq. (51) of Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param i matrix index
     * @param j matrix index
     * @param k order index
     * @param eta als(M)/als(mu)
     * @return function value
     */
    double f_r(unsigned int nf, unsigned int i, unsigned int j, int k, double eta);
    
    /**
     * @brief auxiliary function g - eq. (52) of Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param i matrix index
     * @param p matrix index
     * @param j matrix index
     * @param k order index
     * @param l order index
     * @param eta als(M)/als(mu)
     * @return function value
     */
    double f_g(unsigned int nf, unsigned int i, unsigned int p, unsigned int j, int k, int l, double eta);
  
    /**
     * @brief auxiliary function h - eq. (53) of Huber, Lunghi, Misiak, Wyler, hep-ph/0512066
     * @param i matrix index
     * @param p matrix index
     * @param q matrix index
     * @param j matrix index
     * @param k order index
     * @param l order index
     * @param m order index
     * @param eta als(M)/als(mu)
     * @return function value
     */
    double f_h(unsigned int nf, unsigned int i, unsigned int p, unsigned int q, unsigned int j, int k, int l, int m, double eta);
    
    // operators to include         {C, P, M, L, Q, b}            Huber et al., hep-ph/0512066
    unsigned int nops, nfmin, nfmax;
    std::string blocks;

    gslpp::matrix<gslpp::complex> evec, evec_i, js, h, gg, s_s, jssv, jss, jv, vij;
    gslpp::vector<gslpp::complex> eval;
    double alsMZ_cache, Mz_cache;
};

#endif /* EVOLDF1_H */
