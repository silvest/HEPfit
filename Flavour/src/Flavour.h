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
#include <HeffDF2.h>

using namespace gslpp;

class Flavour : public ThObsType {
public:

    Flavour(const StandardModel& SM_i, StandardModelMatching& SM_Matching_i) : ThObsType(SM_i), 
            HDF2(SM_i, SM_Matching_i) {
    };

    const HeffDF2& getHDF2() const {
        return HDF2;
    }
    
    vector<complex>** ComputeCoeff(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeff(mu, scheme);
    }

private:
    HeffDF2 HDF2;
};

#endif	/* FLAVOUR_H */
