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
#include <boost/tuple/tuple.hpp>

class Flavour {
public:

    Flavour(const StandardModel& SM_i) : 
        HDF2(SM_i), 
        HDB1(SM_i), 
        HDS1(SM_i)
    {
        myMVll_BdKstarmu = new MVll(SM_i, StandardModel::B_D, StandardModel::K_star, StandardModel::MU);
        myMVll_BdKstarel = new MVll(SM_i, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON);
        myMVll_Bsphimu = new MVll(SM_i, StandardModel::B_S, StandardModel::PHI, StandardModel::MU);
        myMVll_Bsphiel = new MVll(SM_i, StandardModel::B_S, StandardModel::PHI, StandardModel::ELECTRON);
        myMPll_BpKmu = new MPll(SM_i, StandardModel::B_P, StandardModel::K_P, StandardModel::MU);
        myMPll_BpKel = new MPll(SM_i, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON);
        
        update_BdKstarmu = true;
        update_BdKstarel = true;
        update_Bsphimu = true;
        update_Bsphiel = true;
        update_BpKmu = true;
        update_BpKel = true;
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
    
    gslpp::vector<gslpp::complex>** ComputeCoeffBd(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffBd(mu, scheme);
    }

    gslpp::vector<gslpp::complex>** ComputeCoeffBs(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffBs(mu, scheme);
    }

    gslpp::vector<gslpp::complex>** ComputeCoeffdd(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffdd(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffK(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffK(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffmK(double mu, schemes scheme = NDR) {
        return HDF2.ComputeCoeffmK(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffDS1PP(double mu, schemes scheme = NDR) {
        return HDS1.ComputeCoeffDS1PP(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffDS1pnunu() {
        return HDS1.ComputeCoeffDS1pnunu();
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffDS1mumu() {
        return HDS1.ComputeCoeffDS1mumu();
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffsmumu(double mu, schemes scheme = NDR) {
        return HDB1.ComputeCoeffsmumu(mu, scheme);
    }
    
    
    gslpp::vector<gslpp::complex>** ComputeCoeffdmumu(double mu, schemes scheme = NDR) {
        return HDB1.ComputeCoeffdmumu(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffbtaunu() {
        return HDB1.ComputeCoeffbtaunu();
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffsnunu() {
        return HDB1.ComputeCoeffsnunu();
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffdnunu() {
        return HDB1.ComputeCoeffdnunu();
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffsgamma(double mu, schemes scheme = NDR) {
        return HDB1.ComputeCoeffsgamma(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffprimesgamma(double mu, schemes scheme = NDR) {
        return HDB1.ComputeCoeffprimesgamma(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffBMll(double mu, schemes scheme = NDR) {
        return HDB1.ComputeCoeffBMll(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** ComputeCoeffprimeBMll(double mu, schemes scheme = NDR) {
        return HDB1.ComputeCoeffprimeBMll(mu, scheme);
    }
    
    MVll* getMVll(StandardModel::meson meson_i, StandardModel::meson vector_i, StandardModel::lepton lep_i) {
        if (meson_i == StandardModel::B_D && vector_i == StandardModel::K_star && lep_i == StandardModel::MU) return myMVll_BdKstarmu;
        if (meson_i == StandardModel::B_D && vector_i == StandardModel::K_star && lep_i == StandardModel::ELECTRON) return myMVll_BdKstarel;
        if (meson_i == StandardModel::B_S && vector_i == StandardModel::PHI && lep_i == StandardModel::MU) return myMVll_Bsphimu;
        if (meson_i == StandardModel::B_S && vector_i == StandardModel::PHI && lep_i == StandardModel::ELECTRON) return myMVll_Bsphiel;
        else throw std::runtime_error("Flavour: Decay channel not implemented.");
    }
    
    MPll* getMPll(StandardModel::meson meson_i, StandardModel::meson pseudoscalar_i, StandardModel::lepton lep_i) {
        if (meson_i == StandardModel::B_P && pseudoscalar_i == StandardModel::K_P && lep_i == StandardModel::MU) return myMPll_BpKmu;
        if (meson_i == StandardModel::B_P && pseudoscalar_i == StandardModel::K_P && lep_i == StandardModel::ELECTRON) return myMPll_BpKel;
        else throw std::runtime_error("Flavour: Decay channel not implemented.");
    }
    
    void setUpdateFlag(StandardModel::meson meson_i, StandardModel::meson meson_j, StandardModel::lepton lep_i, bool updated_i){
        if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::MU) {update_BdKstarmu = updated_i; return;}
        if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::ELECTRON) {update_BdKstarel = updated_i; return;}
        if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::MU) {update_Bsphimu = updated_i; return;}
        if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::ELECTRON) {update_Bsphiel = updated_i; return;}
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::MU) {update_BpKmu = updated_i; return;}
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::ELECTRON) {update_BpKel = updated_i; return;}
        else throw std::runtime_error("Flavour: Wrong update flag requested.");
    }
    
    bool getUpdateFlag(StandardModel::meson meson_i, StandardModel::meson meson_j, StandardModel::lepton lep_i){
        if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::MU) return update_BdKstarmu;
        if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::ELECTRON) return update_BdKstarel;
        if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::MU) return update_Bsphimu;
        if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::ELECTRON) return update_Bsphiel;
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::MU) return update_BpKmu;
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::ELECTRON) return update_BpKel;
        else throw std::runtime_error("Flavour: Wrong update flags requested.");
    }
    
    void setSMupdated(){
        update_BdKstarmu = true;
        update_BdKstarel = true;
        update_Bsphimu = true;
        update_Bsphiel = true;
        update_BpKmu = true;
        update_BpKel = true;
    }
    
private:
    
    HeffDF2 HDF2;
    HeffDB1 HDB1;
    HeffDS1 HDS1;
    MVll* myMVll_BdKstarmu;
    MVll* myMVll_BdKstarel;
    MVll* myMVll_Bsphimu;
    MVll* myMVll_Bsphiel;
    MPll* myMPll_BpKmu;
    MPll* myMPll_BpKel;
    bool update_BdKstarmu;
    bool update_BdKstarel;
    bool update_Bsphimu;
    bool update_Bsphiel;
    bool update_BpKmu;
    bool update_BpKel;
};

#endif	/* FLAVOUR_H */