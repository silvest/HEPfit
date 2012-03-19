/* 
 * File:   HeffDF1bnlep.h
 * Author: stefano
 *
 * Created on 11 ottobre 2011, 15.19
 */

#ifndef HEFFDF1BNLEP_H
#define	HEFFDF1BNLEP_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>
#include "EvolDF1nlep.h"
#include <sstream>
//#include <QCD.h>

using namespace gslpp;

class HeffDF1bnlep {
public:
    HeffDF1bnlep(const StandardModel & SM, StandardModelMatching& modelmatching);
    virtual ~HeffDF1bnlep();
    

    gslpp::vector<complex>** ComputeCoeffBnlep00(double mu, schemes scheme = NDR);
    gslpp::vector<complex>** ComputeCoeffBnlep10(double mu, schemes scheme = NDR);
    gslpp::vector<complex>** ComputeCoeffBnlep01(double mu, schemes scheme = NDR);
    gslpp::vector<complex>** ComputeCoeffBnlep11(double mu, schemes scheme = NDR);
    
    WilsonCoefficient getCoeffbnlep00() const {
        return coeffbnlep00;
    }
    
    WilsonCoefficient getCoeffbnlep10() const {
        return coeffbnlep01;
    }
    
    WilsonCoefficient getCoeffbnlep01() const {
        return coeffbnlep10;
    }
    
    WilsonCoefficient getCoeffbnlep11() const {
        return coeffbnlep11;
    }
    
    EvolDF1nlep getUDF1() const {
        return u;
    }

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;
    ModelMatching& modelmatching;
    
    WilsonCoefficient coeffbnlep00qcd, coeffbnlep00;
    WilsonCoefficient coeffbnlep10qcd, coeffbnlep10;
    WilsonCoefficient coeffbnlep01, coeffbnlep01A, coeffbnlep01B, coeffbnlep00CC;
    WilsonCoefficient coeffbnlep11, coeffbnlep11A, coeffbnlep11B, coeffbnlep10CC;
    EvolDF1nlep u;
    
    gslpp::vector<complex> bnlep, bnlep2, bnlepCC;
};

#endif	/* HEFFDF1BNLEP_H */


