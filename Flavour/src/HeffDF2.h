/* 
 * File:   HeffDF2.h
 * Author: silvest
 *
 * Created on April 28, 2011, 4:34 PM
 */

#ifndef HEFFDF2_H
#define	HEFFDF2_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>
#include <EvolDF2.h>

using namespace gslpp;

class HeffDF2 {
public:
    HeffDF2(const StandardModel & SM, StandardModelMatching & SM_Matching);
    virtual ~HeffDF2();
//    vector<complex>& getCoeff(double mu, schemes scheme = NDR, orders order = NLO); 

    void ChangeScheme(schemes schout, schemes schin, orders order);
    
    vector<complex>** ComputeCoeff(double mu, schemes scheme = NDR);
    
    matrix<double> AnomalousDimension(orders order, unsigned int nf = 0) const;

    WilsonCoefficient getCoeff() const {
        return coeff;
    }

    EvolDF2 getUDF2() const {
        return u;
    }

    const StandardModel& GetModel() const {
        return model;
    }

private:
    // Magic Number
    // c and d are the coefficient of als(mu) e als(M)
    // first index number of flavours
    double b[5][5][5], c[3][5][5][5], d[3][5][5][5];
    matrix<double> drNDRLRI;
    const StandardModel& model;
    StandardModelMatching& modelmatching;
    WilsonCoefficient coeff;
    EvolDF2 u;
};

#endif	/* HEFFDF2_H */
