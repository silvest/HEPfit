/* 
 * File:   HeffDF2.h
 * Author: silvest
 *
 * Created on April 28, 2011, 4:34 PM
 */

#ifndef HEFFDF2_H
#define	HEFFDF2_H

#include <StandardModel.h>
#include <WilsonCoefficient.h>
#include <EvolDF2.h>

using namespace gslpp;

class HeffDF2 {
public:
    HeffDF2(const StandardModel & SM);
    virtual ~HeffDF2();
//    vector<complex>& getCoeff(double mu, schemes scheme = NDR, orders order = NLO); 
    void ChangeScheme(schemes schout, schemes schin, orders order);
    vector<complex>** Coeff(double mu, schemes scheme);
    matrix<double> AnomalousDimension(orders order, unsigned int nf = 0) const;

private:
    // Magic Number
    // c and d are the coefficient of als(mu) e als(M)
    // first index number of flavours
    double b[5][5][5], c[3][5][5][5], d[3][5][5][5];
    matrix<double> drNDRLRI;
    const StandardModel& model;
    std::vector<WilsonCoefficient> mcDF2;
    WilsonCoefficient coeffDF2;
    EvolDF2 uDF2;
};

#endif	/* HEFFDF2_H */
