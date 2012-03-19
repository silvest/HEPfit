/* 
 * File:   HeffDF1bsg.h
 * Author: stefano
 *
 * Created on 15 settembre 2011, 15.36
 */

#ifndef HEFFDF1_H
#define	HEFFDF1_H

#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <WilsonCoefficient.h>
#include "EvolDF1bsg.h"
#include <sstream>

using namespace gslpp;

class HeffDF1bsg {
public:
    HeffDF1bsg(const StandardModel & SM, StandardModelMatching & SM_Matching);
    virtual ~HeffDF1bsg();
    
    gslpp::vector<complex>** ComputeCoeffBsg(double mu, schemes scheme = NDR);
    
    WilsonCoefficient getCoeffbsg() const {
        return coeffbsg;
    }
    
    EvolDF1bsg getUDF1() const {
        return u;
    }

    const StandardModel& GetModel() const {
        return model;
    }
    
private :
    const StandardModel& model;
    ModelMatching& modelmatching;
    WilsonCoefficient coeffbsg;
    EvolDF1bsg u;
    
};

#endif	/* HEFFDF1_H */

