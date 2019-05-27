/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <functional>
#include "StandardModel.h"
#include "Flavour.h"
#include "MVll.h"
#include "MPll.h"
#include "HeffDF2.h"
#include "HeffDS1.h"
#include "HeffDB1.h"
#include "MVgamma.h"
#include "MVlnu.h"
#include "MPlnu.h"

Flavour::Flavour(const StandardModel& SM_i)
: mySM(SM_i)
{
    dispersion = false;
    CLNflag = false;
    btocNPpmflag = false;
};

HeffDF2& Flavour::getHDF2() const
{
    return *getPtr<HeffDF2>(HDF2);
}

HeffDS1& Flavour::getHDS1() const
{
    return *getPtr<HeffDS1>(HDS1);
}

HeffDB1& Flavour::getHDB1() const
{
    return *getPtr<HeffDB1>(HDB1);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffBd(double mu, schemes scheme) const
{
    return getPtr<HeffDF2>(HDF2)->ComputeCoeffBd(mu, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffBs(double mu, schemes scheme, bool SM) const
{
    return getPtr<HeffDF2>(HDF2)->ComputeCoeffBs(mu, scheme, SM);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffdd(double mu, schemes scheme) const
{
    return getPtr<HeffDF2>(HDF2)->ComputeCoeffdd(mu, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffK(double mu, schemes scheme) const
{
    return getPtr<HeffDF2>(HDF2)->ComputeCoeffK(mu, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffmK(double mu, schemes scheme) const
{
    return getPtr<HeffDF2>(HDF2)->ComputeCoeffmK(mu, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffDS1PP(double mu, schemes scheme) const
{
    return getPtr<HeffDS1>(HDS1)->ComputeCoeffDS1PP(mu, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffDS1pnunu() const
{
    return getPtr<HeffDS1>(HDS1)->ComputeCoeffDS1pnunu();
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffDS1mumu() const
{
    return getPtr<HeffDS1>(HDS1)->ComputeCoeffDS1mumu();
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffsmumu(double mu, schemes scheme) const
{
    return getPtr<HeffDB1>(HDB1)->ComputeCoeffsmumu(mu, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffdmumu(double mu, schemes scheme) const
{
    return getPtr<HeffDB1>(HDB1)->ComputeCoeffdmumu(mu, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffbtaunu(QCD::meson meson_i) const
{
    return getPtr<HeffDB1>(HDB1)->ComputeCoeffbtaunu(meson_i);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffsnunu() const
{
    return getPtr<HeffDB1>(HDB1)->ComputeCoeffsnunu();
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffdnunu() const
{
    return getPtr<HeffDB1>(HDB1)->ComputeCoeffdnunu();
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffsgamma(double mu, bool noSM, schemes scheme) const
{
    return getPtr<HeffDB1>(HDB1)->ComputeCoeffsgamma(mu, noSM, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffprimesgamma(double mu, schemes scheme) const
{
    return getPtr<HeffDB1>(HDB1)->ComputeCoeffprimesgamma(mu, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffBMll(double mu, QCD::lepton lepton, bool noSM, schemes scheme) const
{
    return getPtr<HeffDB1>(HDB1)->ComputeCoeffBMll(mu, lepton, noSM, scheme);
}

gslpp::vector<gslpp::complex>** Flavour::ComputeCoeffprimeBMll(double mu, QCD::lepton lepton, schemes scheme) const
{
    return getPtr<HeffDB1>(HDB1)->ComputeCoeffprimeBMll(mu, lepton, scheme);
}

MVll& Flavour::getMVll(QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) const
{
//    std::reference_wrapper<std::shared_ptr<MVll> > x(MVll_BdKstarmu);
//    if (meson_i == StandardModel::B_D && vector_i == StandardModel::K_star && lep_i == StandardModel::MU) x = std::ref(MVll_BdKstarmu);
//    else if (meson_i == StandardModel::B_D && vector_i == StandardModel::K_star && lep_i == StandardModel::ELECTRON) x = std::ref(MVll_BdKstarel);
//    else if (meson_i == StandardModel::B_P && vector_i == StandardModel::K_star_P && lep_i == StandardModel::MU) x = std::ref(MVll_BpKstarmu);
//    else if (meson_i == StandardModel::B_P && vector_i == StandardModel::K_star_P && lep_i == StandardModel::ELECTRON) x = std::ref(MVll_BpKstarel);
//    else if (meson_i == StandardModel::B_S && vector_i == StandardModel::PHI && lep_i == StandardModel::MU) x = std::ref(MVll_Bsphimu);
//    else if (meson_i == StandardModel::B_S && vector_i == StandardModel::PHI && lep_i == StandardModel::ELECTRON) x = std::ref(MVll_Bsphiel);
//    else throw std::runtime_error("Flavour: Decay channel not implemented.");
//    return *getPtr<MVll>(x.get(), meson_i, vector_i, lep_i);
//    std::vector<int> key({meson_i, vector_i, lep_i});
//    if(MVllMap.find(key)==MVllMap.end()) {
//        MVllMap.insert(std::make_pair(key,std::shared_ptr<MVll>(new MVll(mySM,meson_i, vector_i, lep_i))));
//    }
//    return *MVllMap.at(key);
    return getM<MVll>(MVllMap, meson_i, vector_i, lep_i);
}

MVgamma& Flavour::getMVgamma(QCD::meson meson_i, QCD::meson vector_i) const
{
//    std::reference_wrapper<std::shared_ptr<MVgamma> > x(MVgamma_BdKstgamma);
//    if (meson_i == StandardModel::B_D && vector_i == StandardModel::K_star) x = std::ref(MVgamma_BdKstgamma);
//    else if (meson_i == StandardModel::B_P && vector_i == StandardModel::K_star_P) x = std::ref(MVgamma_BpKstgamma);
//    else if (meson_i == StandardModel::B_S && vector_i == StandardModel::PHI) x = std::ref(MVgamma_Bsphigamma);
//    else throw std::runtime_error("Flavour: Decay channel not implemented.");
//    return *getPtr<MVgamma>(x.get(), meson_i, vector_i);
//    std::vector<int> key({meson_i, vector_i});
//    if(MVgammaMap.find(key)==MVgammaMap.end()) {
//        MVgammaMap.insert(std::make_pair(key,std::shared_ptr<MVgamma>(new MVgamma(mySM,meson_i, vector_i))));
//    }
//    return *MVgammaMap.at(key);
    return getM<MVgamma>(MVgammaMap, meson_i, vector_i);
}

MVlnu& Flavour::getMVlnu(QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) const
{
//    std::reference_wrapper<std::shared_ptr<MVlnu> > x(MVlnu_BdbarDstartaunu);
//    if (meson_i == StandardModel::B_D && vector_i == StandardModel::D_star_P && lep_i == StandardModel::TAU) x = std::ref(MVlnu_BdbarDstartaunu);
//    else if (meson_i == StandardModel::B_D && vector_i == StandardModel::D_star_P && lep_i == StandardModel::MU) x = std::ref(MVlnu_BdbarDstarmunu);
//    else if (meson_i == StandardModel::B_D && vector_i == StandardModel::D_star_P && lep_i == StandardModel::ELECTRON) x = std::ref(MVlnu_BdbarDstarelnu);
//    else throw std::runtime_error("Flavour: Decay channel not implemented.");
//    return *getPtr<MVlnu>(x.get(), meson_i, vector_i, lep_i);
//    std::vector<int> key({meson_i, vector_i, lep_i});
//    if(MVlnuMap.find(key)==MVlnuMap.end()) {
//        MVlnuMap.insert(std::make_pair(key,std::shared_ptr<MVlnu>(new MVlnu(mySM,meson_i, vector_i, lep_i))));
//    }
//    return *MVlnuMap.at(key);
    return getM<MVlnu>(MVlnuMap, meson_i, vector_i, lep_i);

}

MPlnu& Flavour::getMPlnu(QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) const
{
//    std::reference_wrapper<std::shared_ptr<MPlnu> > x(MPlnu_BdbarDtaunu);
//    if (meson_i == StandardModel::B_D && vector_i == StandardModel::D_P && lep_i == StandardModel::TAU) x = std::ref(MPlnu_BdbarDtaunu);
//    else if (meson_i == StandardModel::B_D && vector_i == StandardModel::D_P && lep_i == StandardModel::MU) x = std::ref(MPlnu_BdbarDmunu);
//    else if (meson_i == StandardModel::B_D && vector_i == StandardModel::D_P && lep_i == StandardModel::ELECTRON) x = std::ref(MPlnu_BdbarDelnu);
//    else throw std::runtime_error("Flavour: Decay channel not implemented.");
//    return *getPtr<MPlnu>(x.get(), meson_i, vector_i, lep_i);
//    std::vector<int> key({meson_i, vector_i, lep_i});
//    if(MPlnuMap.find(key)==MPlnuMap.end()) {
//        MPlnuMap.insert(std::make_pair(key,std::shared_ptr<MPlnu>(new MPlnu(mySM,meson_i, vector_i, lep_i))));
//    }
//    return *MPlnuMap.at(key);
    return getM<MPlnu>(MPlnuMap, meson_i, vector_i, lep_i);
}

MPll& Flavour::getMPll(QCD::meson meson_i, QCD::meson vector_i, QCD::lepton lep_i) const
{
//    std::reference_wrapper<std::shared_ptr<MPll> > x(MPll_BpKmu);
//    if (meson_i == QCD::B_P && pseudoscalar_i == QCD::K_P && lep_i == QCD::MU) x = std::ref(MPll_BpKmu);
//    else if (meson_i == QCD::B_P && pseudoscalar_i == QCD::K_P && lep_i == QCD::ELECTRON) x = std::ref(MPll_BpKel);
//    else if (meson_i == QCD::B_D && pseudoscalar_i == QCD::K_0 && lep_i == QCD::MU) x = std::ref(MPll_B0Kmu);
//    else if (meson_i == QCD::B_D && pseudoscalar_i == QCD::K_0 && lep_i == QCD::ELECTRON) x = std::ref(MPll_B0Kel);
//    else throw std::runtime_error("Flavour: Decay channel not implemented.");
//    return *getPtr<MPll>(x.get(), meson_i, pseudoscalar_i, lep_i);
//    std::vector<int> key({meson_i, vector_i, lep_i});
//    if(MPllMap.find(key)==MPllMap.end()) {
//        MPllMap.insert(std::make_pair(key,std::shared_ptr<MPll>(new MPll(mySM,meson_i, vector_i, lep_i))));
//    }
//    return *MPllMap.at(key);
    return getM<MPll>(MPllMap, meson_i, vector_i, lep_i);
}

template <typename T, typename... Args>
T& Flavour::getM(std::map<std::vector<int>,std::shared_ptr<T> >& map, Args ... args) const
{
    std::vector<int> key({args...});
    if(map.find(key)==map.end()) {
        map.insert(std::make_pair(key,std::shared_ptr<T>(new T(mySM,args...))));
    }
    return *map.at(key);
}

void Flavour::setUpdateFlag(QCD::meson meson_i, QCD::meson meson_j, QCD::lepton lep_i, bool updated_i) const
{
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::MU) update_BdKstarmu = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::ELECTRON) update_BdKstarel = updated_i;
//    else if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::MU) update_BpKstarmu = updated_i;
//    else if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::ELECTRON) update_BpKstarel = updated_i;
//    else if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::MU) update_Bsphimu = updated_i;
//    else if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::ELECTRON) update_Bsphiel = updated_i;
//    else if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::MU) update_BpKmu = updated_i;
//    else if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::ELECTRON) update_BpKel = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_0 && lep_i == StandardModel::MU) update_B0Kmu = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_0 && lep_i == StandardModel::ELECTRON) update_B0Kel = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::NOLEPTON) update_BdKstgamma = updated_i;
//    else if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::NOLEPTON) update_BpKstgamma = updated_i;
//    else if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::NOLEPTON) update_Bsphigamma = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_star_P && lep_i == StandardModel::TAU) update_BdDstartaunu = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_star_P && lep_i == StandardModel::MU) update_BdDstarmunu = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_star_P && lep_i == StandardModel::ELECTRON) update_BdDstarelnu = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_P && lep_i == StandardModel::TAU) update_BdDtaunu = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_P && lep_i == StandardModel::MU) update_BdDmunu = updated_i;
//    else if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_P && lep_i == StandardModel::ELECTRON) update_BdDelnu = updated_i;
//    else throw std::runtime_error("Flavour: Wrong update flag requested.");
    flagUpdateMap[std::vector<int>({meson_i, meson_j, lep_i})] = updated_i;
}

bool Flavour::getUpdateFlag(QCD::meson meson_i, QCD::meson meson_j, QCD::lepton lep_i) const
{
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::MU) return update_BdKstarmu;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::ELECTRON) return update_BdKstarel;
//    if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::MU) return update_BpKstarmu;
//    if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::ELECTRON) return update_BpKstarel;
//    if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::MU) return update_Bsphimu;
//    if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::ELECTRON) return update_Bsphiel;
//    if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::MU) return update_BpKmu;
//    if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_P && lep_i == StandardModel::ELECTRON) return update_BpKel;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_0 && lep_i == StandardModel::MU) return update_B0Kmu;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_0 && lep_i == StandardModel::ELECTRON) return update_B0Kel;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::K_star && lep_i == StandardModel::NOLEPTON) return update_BdKstgamma;
//    if (meson_i == StandardModel::B_P && meson_j == StandardModel::K_star_P && lep_i == StandardModel::NOLEPTON) return update_BpKstgamma;
//    if (meson_i == StandardModel::B_S && meson_j == StandardModel::PHI && lep_i == StandardModel::NOLEPTON) return update_Bsphigamma;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_star_P && lep_i == StandardModel::TAU) return update_BdDstartaunu;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_star_P && lep_i == StandardModel::MU) return update_BdDstarmunu;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_star_P && lep_i == StandardModel::ELECTRON) return update_BdDstarelnu;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_P && lep_i == StandardModel::TAU) return update_BdDtaunu;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_P && lep_i == StandardModel::MU) return update_BdDmunu;
//    if (meson_i == StandardModel::B_D && meson_j == StandardModel::D_P && lep_i == StandardModel::ELECTRON) return update_BdDelnu;
//    else throw std::runtime_error("Flavour: Wrong update flags requested.");
//
    std::vector<int> key({meson_i, meson_j, lep_i});
    if(flagUpdateMap.find(key)==flagUpdateMap.end()) {
        // Default value
        flagUpdateMap[key]=true;
    }
    return flagUpdateMap.at(key);
}

template<typename T, typename... Args> std::shared_ptr<T>& Flavour::getPtr(std::shared_ptr<T>& x, Args ... args) const
{
    if (x.get() == nullptr)
        x.reset(new T(mySM, args...));
    return x;
}

void Flavour::setSMupdated() const
{
    for(auto &x : flagUpdateMap) x.second = true;
}
