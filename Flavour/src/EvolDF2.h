/* 
 * File:   EvolDF2.h
 * Author: marco
 *
 * Created on May 20, 2011, 3:55 PM
 */

#ifndef EVOLDF2_H
#define	EVOLDF2_H

#include <RGEvolutor.h>

class EvolDF2 : public RGEvolutor {
public:
    EvolDF2(unsigned int dim, schemes scheme, orders order,
            const QCD& model);
    virtual ~EvolDF2();
    matrix<double> AnomalousDimension(orders order, unsigned int nf) const;
    matrix<double>& EvolDF2::Df2Evol(double mu, double M, orders order, 
            schemes scheme = NDR);
private:
    double b[5][5][5], c[3][5][5][5], d[3][5][5][5];
    const StandardModel& model;
};

#endif	/* EVOLDF2_H */

