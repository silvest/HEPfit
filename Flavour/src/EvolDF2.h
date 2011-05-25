/* 
 * File:   EvolDF2.h
 * Author: marco
 *
 * Created on May 20, 2011, 3:55 PM
 */

#ifndef EVOLDF2_H
#define	EVOLDF2_H

#include <RGEvolutor.h>
#include <StandardModel.h>
#include <sstream>

class EvolDF2 : public RGEvolutor {
public:
    EvolDF2(unsigned int dim, schemes scheme, orders order,
            const StandardModel& model);
    virtual ~EvolDF2();
    matrix<double> AnomalousDimension(orders order, unsigned int nf) const;
    matrix<double>& Df2Evol(double mu, double M, orders order, 
            schemes scheme = NDR);
private:
        matrix<double> Df2Evol(double mu, double M, double nf, orders order, 
        schemes scheme);
        double a[5], b[5][5][5], c[3][5][5][5], d[3][5][5][5];
    const StandardModel& model;
};

#endif	/* EVOLDF2_H */

