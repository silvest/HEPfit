/* 
 * File:   EvolDF1nlep.h
 * Author: Mauro_87
 *
 * Created on 11 ottobre 2011, 15.48
 */

#ifndef EVOLDF1NLEP_H
#define	EVOLDF1NLEP_H


#include <RGEvolutor.h>
#include <StandardModel.h>
#include <sstream>

using namespace gslpp;

class EvolDF1nlep : public RGEvolutor {
/**
 * @class EvolDF1nlep
 * @brief \f$ |\Delta F = 1 | \f$ Evolutor Class
 * @details This evolutor is properly written to study \f$ |\Delta F = 1 | \f$ 
 * processes such as no leptonic weak decays of the B meson; 
 * it is implemented for the evolution of the 10 Wilson coefficients (related to the 
 * 2 current x current + 4 QCD penguins + 4 em penguins)in the standard basis 
 * at the NLO in \f$ \alpha_{strong} \f$ for QCD corrections and at the NLO 
 * in \f$ \alpha_{em} \f$ for QED ones; 
 * principal reference: hep-ph/9512380v1
 */
    public:
    /**
     * @brief EvolDF1nlep constructor
     * @param dim an unsigned integer  for the dimension of the evolutor
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @param order an enum "orders" for the order of QCD perturbation theory of the evolutor
     * @param order_ew an enum "orders_ew" for the order of QED perturbation theory of the evolutor 
     * @param model an object of StandardModel class 
     */
    EvolDF1nlep(unsigned int dim,  schemes scheme,orders order, orders_ew order_ew, const StandardModel& model);
    /**
     * @brief EvolDF1nlep destructor 
     */
    virtual ~EvolDF1nlep();
    /**
     * @brief a method returning the anomalous dimension matrix given in the standard basis
     * @param order an enum "orders" for the order of QCD perturbation theory of the ADM
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM related to QCD corrections at the order LO/NLO in the standard basis 
     */
    matrix<double> AnomalousDimension_nlep_S(orders order, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief a method returning the anomalous dimension matrix given in the standard basis
     * @param order an enum "orders" for the order of QED perturbation theory of the ADM
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM related to QED corrections at the order LO/NLO in the standard basis 
     */
    matrix<double> AnomalousDimension_nlep_EM(orders order, unsigned int n_u, unsigned int n_d) const;
    /**
     * @brief a method returning the evolutor related to the high scale \f$ M \f$ and the low scale \f$ \mu \f$
     * @param mu a double for the low scale of the evolution
     * @param M a double for the high scale of the evolution
     * @param order an enum "orders" for the order of QCD perturbation theory of the evolutor 
     * @param order_ew an enum "orders_ew" for the order of QED perturbation theory of the evolutor
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @return the evolutor \f$ U (\mu , M) \f$
     */
    matrix<double>& Df1Evolnlep(double mu, double M, orders order, orders_ew order_ew, schemes scheme = NDR);
    /**
     * @brief a method returning the matrix threshold for the QCD penguins at the NLO
     * @param nf a double for the active number of flavors
     * @return matrix threshold for QCD penguin operators
     */
    matrix<double> Df1threshold_deltarsT(double nf) const;
    /**
     * @brief a method returning the matrix threshold for the QED penguins at the NLO
     * @param nf a double for the active number of flavors
     * @return QED matrix threshold for QED penguin operators
     */
    matrix<double> Df1threshold_deltareT(double nf) const;    
    
    private:
    /**
     * @param nu an unsigned integer for the up-type number of d.o.f.
     * @param nu an unsigned integer for the down-type number of d.o.f.
     */
    int nu,nd;
    /**
     * @param a array of double for the magic numbers of the evolutor ( LO evolution )
     * @param b array of double for the magic numbers of the evolutor ( LO evolution ) 
     * @param c array of double for the magic numbers of the evolutor ( QED corrections proportional to \f$ \alpha_{em} / \alpha_{strong}(\mu) \f$ )
     * @param d array of double for the magic numbers of the evolutor ( QED corrections proportional to \f$ \alpha_{em} / \alpha_{strong}(M) \f$ )
     * @param m array of double for the magic numbers of the evolutor ( NLO evolution, associated to \f$ \alpha_{strong}(\mu) \f$ )
     * @param n array of double for the magic numbers of the evolutor ( NLO evolution, associated to \f$ \alpha_{strong}(M) \f$ )
     * @param o array of double for the magic numbers of the evolutor ( QED corrections proportional to \f$ \alpha_{em} \f$ )
     * @param p array of double for the magic numbers of the evolutor ( QED corrections proportional to \f$ \alpha_{em} \f$ )
     * @param q array of double for the magic numbers of the evolutor ( QED corrections proportional to \f$ \alpha_{strong}(M) / \alpha_{strong}(\mu) * \alpha_{em} \f$ )
     * @param r array of double for the magic numbers of the evolutor ( QED corrections proportional to \f$ \alpha_{em} \f$ )
     * @param s array of double for the magic numbers of the evolutor ( QED corrections proportional to \f$ \alpha_{em} \f$ )
     * @param t array of double for the magic numbers of the evolutor ( QED corrections proportional to \f$ \alpha_{strong}(\mu) / \alpha_{strong}(M) * \alpha_{em} \f$ )
     */
    double a[3][10], b[3][10][10][10], c[3][10][10][10], d[3][10][10][10],
           m[3][10][10][10], n[3][10][10][10], o[3][10][10][10], 
           p[3][10][10][10], q[3][10][10][10], r[3][10][10][10], 
           s[3][10][10][10], t[3][10][10][10];
    const StandardModel& model;
    /**
     * @brief a void type method storing properly the magic numbers for the implementation of the evolutor
     * @param mu a double for the low scale of the evolution
     * @param M a double for the high scale of the evolution
     * @param nf a double for the active number of flavors
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     */
    void Df1Evolnlep(double mu, double M, double nf, schemes scheme);
    /**
     * @brief a void type method for the implementation of the NLO threshold effects in the evolutor
     * @param M a double for the high scale of the evolution
     * @param nf a double for the active number of flavors
     */
    void Df1threshold_nlep(double M, double nf);
    gslpp::matrix <gslpp::complex> v, vi, js, h, gg, jv, vij, g_0, k_0, vk_0, k_0vi, g_1,
                       k_11, vk_11, k_11vi, k_12, vh, k_12vi, k12s, s_svi, k_12s, k_12svi,
                       k_13s, vk_13s, vs_s, Gamma_T, Gamma_ew, k_13, vk_13, vg_1, vg_0h, 
                       hvi, vhg_0, vg_0, Gamma_1, s_s, jss, jssv; 
    gslpp::vector<gslpp::complex> e;
    

 };


#endif	/* EVOLDF1NLEP_H */

