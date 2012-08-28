/* 
 * File:   Flavour.h
 * Author: silvest
 *
 * Created on March 29, 2011, 12:49 PM
 */

#ifndef FLAVOUR_H
#define	FLAVOUR_H

#include <ThObsType.h>
#include <StandardModel.h>
#include "HeffDF2.h"

using namespace gslpp;

class Flavour : public ThObsType {
public:

    Flavour(const StandardModel& SM_i) : ThObsType(SM_i), 
            HDF2(SM_i) {
        
        if(!SM_i.IsMatchingInitialized())
            throw "Matching not initialized ";
    };

    const HeffDF2& getHDF2() const {
        return HDF2;
    }
    
    vector<complex>** ComputeCoeffBd(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffBd(mu, scheme);
    }

    vector<complex>** ComputeCoeffBs(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffBs(mu, scheme);
    }

    vector<complex>** ComputeCoeffdd(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffdd(mu, scheme);
    }
    
    vector<complex>** ComputeCoeffK(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffK(mu, scheme);
    }
    
private:
    HeffDF2 HDF2;
};

#endif	/* FLAVOUR_H */
