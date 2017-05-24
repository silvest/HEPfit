/* 
 * Copyright (C) 2012 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModel.h"
#include "Flavour.h"

Flavour::Flavour(const StandardModel& SM_i) : HDF2(SM_i), HDB1(SM_i), HDS1(SM_i),
            MVll_BdKstarmu(SM_i, StandardModel::B_D, StandardModel::K_star, StandardModel::MU),
            MVll_BdKstarel(SM_i, StandardModel::B_D, StandardModel::K_star, StandardModel::ELECTRON),
            MVll_BpKstarmu(SM_i, StandardModel::B_P, StandardModel::K_star_P, StandardModel::MU),
            MVll_BpKstarel(SM_i, StandardModel::B_P, StandardModel::K_star_P, StandardModel::ELECTRON),
            MVll_Bsphimu(SM_i, StandardModel::B_S, StandardModel::PHI, StandardModel::MU),
            MVll_Bsphiel(SM_i, StandardModel::B_S, StandardModel::PHI, StandardModel::ELECTRON),
            MPll_BpKmu(SM_i, StandardModel::B_P, StandardModel::K_P, StandardModel::MU),
            MPll_BpKel(SM_i, StandardModel::B_P, StandardModel::K_P, StandardModel::ELECTRON)
    {
        update_BdKstarmu = true;
        update_BdKstarel = true;
        update_BpKstarmu = true;
        update_BpKstarel = true;
        update_Bsphimu = true;
        update_Bsphiel = true;
        update_BpKmu = true;
        update_BpKel = true;
    };
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffBd(double mu, schemes scheme) const {
        return HDF2.ComputeCoeffBd(mu, scheme);
    }

    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffBs(double mu, schemes scheme) const {
        return HDF2.ComputeCoeffBs(mu, scheme);
    }

    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffdd(double mu, schemes scheme) const {
        return HDF2.ComputeCoeffdd(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffK(double mu, schemes scheme) const {
        return HDF2.ComputeCoeffK(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffmK(double mu, schemes scheme) const {
        return HDF2.ComputeCoeffmK(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffDS1PP(double mu, schemes scheme) const {
        return HDS1.ComputeCoeffDS1PP(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffDS1pnunu() const {
        return HDS1.ComputeCoeffDS1pnunu();
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffDS1mumu() const {
        return HDS1.ComputeCoeffDS1mumu();
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffsmumu(double mu, schemes scheme) const {
        return HDB1.ComputeCoeffsmumu(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffdmumu(double mu, schemes scheme) const {
        return HDB1.ComputeCoeffdmumu(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffbtaunu() const {
        return HDB1.ComputeCoeffbtaunu();
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffsnunu() const {
        return HDB1.ComputeCoeffsnunu();
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffdnunu() const {
        return HDB1.ComputeCoeffdnunu();
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffsgamma(double mu, schemes scheme) const  {
        return HDB1.ComputeCoeffsgamma(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffprimesgamma(double mu, schemes scheme) const {
        return HDB1.ComputeCoeffprimesgamma(mu, scheme);
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffBMll(double mu, QCD::lepton lepton, schemes scheme) const {
        return HDB1.ComputeCoeffBMll(mu, lepton, scheme);
    }
    
    gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffprimeBMll(double mu, QCD::lepton lepton, schemes scheme) const {
        return HDB1.ComputeCoeffprimeBMll(mu, lepton, scheme);
    }
    
    MVll& Flavour::getMVll(unsigned int meson_i, unsigned int vector_i, unsigned int lep_i) const {
        if (meson_i == StandardModel::B_D && vector_i == StandardModel::K_star && lep_i == StandardModel::MU) return MVll_BdKstarmu;
        if (meson_i == StandardModel::B_D && vector_i == StandardModel::K_star && lep_i == StandardModel::ELECTRON) return MVll_BdKstarel;
        if (meson_i == StandardModel::B_P && vector_i == StandardModel::K_star_P && lep_i == StandardModel::MU) return MVll_BpKstarmu;
        if (meson_i == StandardModel::B_P && vector_i == StandardModel::K_star_P && lep_i == StandardModel::ELECTRON) return MVll_BpKstarel;
        if (meson_i == StandardModel::B_S && vector_i == StandardModel::PHI && lep_i == StandardModel::MU) return MVll_Bsphimu;
        if (meson_i == StandardModel::B_S && vector_i == StandardModel::PHI && lep_i == StandardModel::ELECTRON) return MVll_Bsphiel;
        else throw std::runtime_error("Flavour: Decay channel not implemented.");
    }
    
    MPll& Flavour::getMPll(unsigned int meson_i, unsigned int pseudoscalar_i, unsigned int lep_i) const {
        if (meson_i == StandardModel::B_P && pseudoscalar_i == StandardModel::K_P && lep_i == StandardModel::MU) return MPll_BpKmu;
        if (meson_i == StandardModel::B_P && pseudoscalar_i == StandardModel::K_P && lep_i == StandardModel::ELECTRON) return MPll_BpKel;
        else throw std::runtime_error("Flavour: Decay channel not implemented.");
    }
    
    void Flavour::setUpdateFlag(unsigned int meson_i, unsigned int meson_j, unsigned int lep_i, bool updated_i) const {
        if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::MU) {update_BdKstarmu = updated_i; return;}
        if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::ELECTRON) {update_BdKstarel = updated_i; return;}
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::MU) {update_BpKstarmu = updated_i; return;}
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::ELECTRON) {update_BpKstarel = updated_i; return;}
        if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::MU) {update_Bsphimu = updated_i; return;}
        if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::ELECTRON) {update_Bsphiel = updated_i; return;}
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::MU) {update_BpKmu = updated_i; return;}
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::ELECTRON) {update_BpKel = updated_i; return;}
        else throw std::runtime_error("Flavour: Wrong update flag requested.");
    }
    
    bool Flavour::getUpdateFlag(unsigned int meson_i, unsigned int meson_j, unsigned int lep_i) const {
        if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::MU) return update_BdKstarmu;
        if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::ELECTRON) return update_BdKstarel;
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::MU) return update_BpKstarmu;
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::ELECTRON) return update_BpKstarel;
        if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::MU) return update_Bsphimu;
        if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::ELECTRON) return update_Bsphiel;
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::MU) return update_BpKmu;
        if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::ELECTRON) return update_BpKel;
        else throw std::runtime_error("Flavour: Wrong update flags requested.");
    }
    
    void Flavour::setSMupdated() const {
        update_BdKstarmu = true;
        update_BdKstarel = true;
        update_BpKstarmu = true;
        update_BpKstarel = true;
        update_Bsphimu = true;
        update_Bsphiel = true;
        update_BpKmu = true;
        update_BpKel = true;
    }
    