/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef FLAVOUR_H
#define	FLAVOUR_H

#include <StandardModel.h>
#include "HeffDF2.h"
#include "HeffDS1.h"
#include "HeffDB1.h"

using namespace gslpp;

class Flavour {
public:

    Flavour(const StandardModel& SM_i) : HDF2(SM_i), HDS1(SM_i), HDB1(SM_i) {
    };

    const HeffDF2& getHDF2() const {
        return HDF2;
    }
    
    const HeffDS1& getHDS1() const {
        return HDS1;
    }
    
    const HeffDB1& getHDB1() const {
        return HDB1;
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
    
    vector<complex>** ComputeCoeffmK(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffmK(mu, scheme);
    }
    
    vector<complex>** ComputeCoeffDS1PP(double mu, schemes scheme = NDR) {
        return HDS1.ComputeCoeffDS1PP(mu, scheme);
    }
    
    vector<complex>** ComputeCoeffDS1pnunu() {
        return HDS1.ComputeCoeffDS1pnunu();
    }
    
    vector<complex>** ComputeCoeffDS1mumu() {
        return HDS1.ComputeCoeffDS1mumu();
    }
    
    vector<complex>** ComputeCoeffsmumu() {
        return HDB1.ComputeCoeffsmumu();
    }
    
    
    vector<complex>** ComputeCoeffdmumu() {
        return HDB1.ComputeCoeffdmumu();
    }
    
    vector<complex>** ComputeCoeffsnunu() {
        return HDB1.ComputeCoeffdmumu();
    }
    
    vector<complex>** ComputeCoeffdnunu() {
        return HDB1.ComputeCoeffdmumu();
    }
    
    vector<complex>** ComputeCoeffBKstarll(double mu, schemes scheme = NDR) {
        return HDB1.ComputeCoeffBKstarll(mu, scheme);
    }
    
    vector<complex>** ComputeCoeffprimeBKstarll(double mu, schemes scheme = NDR) {
        return HDB1.ComputeCoeffprimeBKstarll(mu, scheme);
    }
    
private:
    
    HeffDF2 HDF2;
    HeffDS1 HDS1;
    HeffDB1 HDB1;
};

#endif	/* FLAVOUR_H */