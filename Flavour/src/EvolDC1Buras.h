/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EVOLDC1BURAS_H
#define	EVOLDC1BURAS_H


#include <RGEvolutor.h>
#include <StandardModel.h>
#include <sstream>

class EvolDC1Buras : public RGEvolutor {
/**
 * @class EvolDC1Buras
 * @brief \f$ |\Delta C = 1 | \f$ Evolutor Class
 * @details This evolutor is properly written for Charm Physics studies;
 * it is implemented for the evolution of the 6 Wilson coefficients
 * (related to the 2 current x current + 4 QCD penguin operators)
 * in the standard (Buras) basis at the NLO in \f$ \alpha_{strong} \f$;
 * principal reference: hep-ph/9512380v1
 */
    public:
    /**
     * @brief EvolDC1Buras constructor
     * @param dim an unsigned integer  for the dimension of the evolutor
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @param order an enum "orders" for the order of perturbation theory of the evolutor
     * @param model an object of StandardModel class 
     */    
    EvolDC1Buras(unsigned int dim_i,  schemes scheme, orders order, const StandardModel& model);
    /**
     * @brief EvolDC1Buras destructor
     */
    virtual ~EvolDC1Buras();
    /**
     * @brief a method returning the evolutor related to the high scale \f$ M \f$ and the low scale \f$ \mu \f$
     * @param mu a double for the low scale of the evolution
     * @param M a double for the high scale of the evolution
     * @param order order an enum "orders" for the order of perturbation theory of the evolutor
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     * @return the evolutor \f$ U (\mu , M) \f$
     */
    gslpp::matrix<double>& DC1EvolBuras(double mu, double M, orders order,  schemes scheme = NDR);
    /**
     * @brief a method returning the anomalous dimension matrix given in the standard basis
     * @param order an enum "orders" for the order of perturbation theory of the ADM
     * @param n_u an unsigned integer for the up-type number of d.o.f.
     * @param n_d an unsigned integer for the down-type number of d.o.f.
     * @return the ADM at the order LO/NLO in the standard basis
     */
    gslpp::matrix<double> AnomalousDimension_DC1_Buras(orders order, unsigned int n_u,
        unsigned int n_d) const;
    /**
     * @brief a method returning the matrix threshold for the QCD penguins at the NLO
     * @return matrix threshold for QCD penguin operators
     */    
    gslpp::matrix<double> StrongThresholds() const;
    
    private:
    /**
     * @param nu an unsigned integer for the up-type number of d.o.f.
     * @param nu an unsigned integer for the down-type number of d.o.f.
     */
    int nu,nd;    
    /**
     * @param a array of double for the magic numbers of the evolutor ( LO evolution )
     * @param b array of double for the magic numbers of the evolutor ( LO evolution )
     * @param c array of double for the magic numbers of the evolutor ( NLO evolution, associated to \f$ \alpha_{strong}(\mu) \f$ )
     * @param d array of double for the magic numbers of the evolutor ( NLO evolution, associated to \f$ \alpha_{strong}(M) \f$ )
     */
    double a[3][10], b[3][10][10][10], c[3][10][10][10], d[3][10][10][10];
    const StandardModel& model;
    /**
     * @brief a void type method storing properly the magic numbers for the implementation of the evolutor  
     * @param mu a double for the low scale of the evolution
     * @param M a double for the high scale of the evolution
     * @param nf a double for the active number of flavors
     * @param scheme an enum "schemes" for the regularization scheme of the evolutor
     */
    void DC1EvolBuras(double mu, double M, double nf, schemes scheme);
    /**
     * @brief a void type method for the implementation of the NLO threshold effects in the evolutor 
     * @param M a double for the intermidiate scale of the threshold
     * @param order an enum "orders" for the order of perturbation theory of the threshold (LO is trivial)
     */
    void DC1PenguinThresholds(double M, orders order);
    gslpp::matrix<gslpp::complex> v, vi, js, h, gg, s_s, jssv, jss, jv, vij;
    gslpp::vector<gslpp::complex> e;
    unsigned int dim;
    double alsMZ_cache;
    double Mz_cache;
    
 };




#endif	/* EVOLDC1BURAS_H */

