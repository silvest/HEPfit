/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EVOLDF2_H
#define	EVOLDF2_H

#include <RGEvolutor.h>
#include <StandardModel.h>
#include <gsl/gsl_sf_dilog.h>

class EvolDF2 : public RGEvolutor {
public:
    /**
     * 
     * @brief constructor
     * @param dim 
     * @param scheme
     * @param order
     * @param model 
     */
    EvolDF2(unsigned int dim_i, schemes scheme, orders order, const StandardModel& model);
    
    /**
     * 
     * @brief destructor
     */
    virtual ~EvolDF2();
    
    /**
     * 
     * @brief ADM in the basis used in Ciuchini et.al. hep-ph/9711402 (basis = 0, default) or in the basis (QVLL, QLR, QSLL) used in Buras et.al. hep-ph/0005183 (basis = 1)
     * @param order
     * @param nf number of active flavours
     * @param basis basis identifier (0: ciuchini, 1: buras)
     * @return Anomalous dimension for DeltaF=2 processes
     */
    gslpp::matrix<double> AnomalousDimension(orders order, unsigned int nf, int basis = 0) const;

    /**
     * 
     * @param mu low energy scale
     * @param M matching scale
     * @param order
     * @param scheme
     * @return the Wilson coefficients evolved from the scale M to the scale mu
     */
    gslpp::matrix<double>& Df2Evol(double mu, double M, orders order, 
            schemes scheme = NDR);
    
    /**
     * 
     * @brief Buras et al, hep-ph/9512380
     * @param mu
     * @return the NLO corrective factor for the charm-charm contribution to the kaon oscillations
     */
    double etacc(double mu) const;
    
    /**
     * 
     * @brief Buras et al, hep-ph/9512380
     * @param mu
     * @return the NLO corrective factor for the charm-top contribution to the kaon oscillations
     */
    double etact(double mu) const;
    
    /**
     * 
     * @brief Buras et al, hep-ph/9512380
     * @param mu
     * @return the NLO corrective factor for the top-top contribution to the kaon oscillations
     */
    double etatt(double mu) const;    
    
private:
    // c and d are the coefficient of als(mu) e als(M)
    // first index number of flavours
    // double b[5][5][5], c[3][5][5][5], d[3][5][5][5];
    
    //double S1tt() const;
    void Df2Evol(double mu, double M, double nf, schemes scheme);
    double a[5];
    double b[5][5][5];
    double c[3][5][5][5];
    double d[3][5][5][5];
    const StandardModel& model;
    unsigned int dim;
    double alsMZ_cache;
    double Mz_cache;
};

#endif	/* EVOLDF2_H */

