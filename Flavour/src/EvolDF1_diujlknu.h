/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EVOLDF1_DIUJLKNU_H
#define EVOLDF1_DIUJLKNU_H

class StandardModel;

#include "RGEvolutor.h"

class EvolDF1_diujlknu : public RGEvolutor {
public:
    /**
     * 
     * @brief constructor
     * @param dim 
     * @param scheme
     * @param order
     * @param model 
     */
    EvolDF1_diujlknu(unsigned int dim_i, schemes scheme, orders order, const StandardModel& model);
    virtual ~EvolDF1_diujlknu();
    
    /**
     * 
     * @brief ADM in the JMS basis, in the order CnueduVLLkkij, CnueduVLRkkij, CnueduSRRkkij, CnueduSRLkkij, CnueduTRRkkij, from 1706.00410
     * @param order
     * @param nf number of active flavours
     * @return Anomalous dimension for \f$ \bar{d}_i u_j \bar{\nu} \ell_k \f$ operators in the JMS basis ordered as CnueduVLLkkij, CnueduVLRkkij, CnueduSRRkkij, CnueduSRLkkij, CnueduTRRkkij
     */
    gslpp::matrix<double> AnomalousDimension(orders order, unsigned int nf) const;

    /**
     * 
     * @brief QED ADM in the JMS basis, in the order CnueduVLLkkij, CnueduVLRkkij, CnueduSRRkkij, CnueduSRLkkij, CnueduTRRkkij, from 1706.00410
     * @param order
     * @param nf number of active flavours
     * @return QED Anomalous dimension for \f$ \bar{d}_i u_j \bar{\nu} \ell_k \f$ operators in the JMS basis ordered as CnueduVLLkkij, CnueduVLRkkij, CnueduSRRkkij, CnueduSRLkkij, CnueduTRRkkij
     */
    gslpp::matrix<double> AnomalousDimension(orders_qed order_qed, unsigned int nf) const;

    /**
     * 
     * @param mu low energy scale
     * @param M matching scale
     * @param order
     * @param scheme
     * @return the Wilson coefficients evolved from the scale M to the scale mu
     */
    gslpp::matrix<double>& Df1diujlknuEvol(double mu, double M, orders order, 
            schemes scheme = NDR);
    

    
private:
    void Df1diujlknuEvol(double mu, double M, double nf, schemes scheme);
    double a[5];
    double b[5][5][5];
    double c[3][5][5][5];
    double d[3][5][5][5];
    const StandardModel& model;
    unsigned int dim;
    double alsMZ_cache;
    double Mz_cache;

};

#endif /* EVOLDF1_DIUJLKNU_H */

