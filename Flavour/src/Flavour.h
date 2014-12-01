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
#include "MVll.h"
#include "MPll.h"

using namespace gslpp;

class Flavour {
public:

    Flavour(const StandardModel& SM_i) : 
        HDF2(SM_i), 
        HDB1(SM_i), 
        HDS1(SM_i)
    {
        myMVll_BdKstarmu = new MVll(SM_i, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
        myMVll_BdKstarel = new MVll(SM_i, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
        myMVll_Bsphimu = new MVll(SM_i, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
        myMVll_Bsphiel = new MVll(SM_i, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
        myMPll_BdKmu = new MPll(SM_i, StandardModel::B_D, StandardModel::K_0, StandardModel::MU);
        myMPll_BdKel = new MPll(SM_i, StandardModel::B_D, StandardModel::K_0, StandardModel::ELECTRON);
        updated = true;
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
    
    MVll* getMVll(StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) {
        if (meson_i == StandardModel::B_D && vector_i == StandardModel::K_star && lep_i == StandardModel::MU) return myMVll_BdKstarmu;
        if (meson_i == StandardModel::B_D && vector_i == StandardModel::K_star && lep_i == StandardModel::ELECTRON) return myMVll_BdKstarel;
        if (meson_i == StandardModel::B_S && vector_i == StandardModel::PHI && lep_i == StandardModel::MU) return myMVll_BdKstarmu;
        if (meson_i == StandardModel::B_S && vector_i == StandardModel::PHI && lep_i == StandardModel::ELECTRON) return myMVll_BdKstarel;
        else throw std::runtime_error("Flavour: Decay channel not implemented.");
    }
    
    MPll* getMPll(StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i) {
        if (meson_i == StandardModel::B_D && pseudoscalar_i == StandardModel::K_0 && lep_i == StandardModel::MU) return myMPll_BdKmu;
        if (meson_i == StandardModel::B_D && pseudoscalar_i == StandardModel::K_0 && lep_i == StandardModel::ELECTRON) return myMPll_BdKel;
        else throw std::runtime_error("Flavour: Decay channel not implemented.");
    }
    
    void setUpdateFlag(bool updated_i){
        updated = updated_i;
    }
    
    bool getUpdateFlag(){
        return updated;
    }
    
private:
    
    HeffDF2 HDF2;
    HeffDB1 HDB1;
    HeffDS1 HDS1;
    MVll* myMVll_BdKstarmu;
    MVll* myMVll_BdKstarel;
    MVll* myMVll_Bsphimu;
    MVll* myMVll_Bsphiel;
    MPll* myMPll_BdKmu;
    MPll* myMPll_BdKel;
    bool updated;
};

#endif	/* FLAVOUR_H */