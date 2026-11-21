// Link with HEPfit and run from the repository root.
#include "AmpDB2.h"
#include "Bdmumu.h"
#include "EvolBsmm.h"
#include "Bsmumu.h"
#include "FlavourWilsonCoefficient_DF2.h"
#include "Mll.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
using Parameters = std::map<std::string, double>;

// The production DF2 matching supplies five independent complex coefficients per channel.
class Model : public FlavourWilsonCoefficient_DF2 {
public:
    double phiD = 0., phiS = 0.;
    const double getPhiBd() const override { return phiD; }
    const double getPhiBs() const override { return phiS; }
};

Parameters inputs(bool withBags)
{
    Parameters p;
    const std::string dir = "Configurations/Temporary/bsll_2026/SMEFT/";
    for (auto file : {"StandardModel.conf", "Flavour.conf"}) {
        std::ifstream config(dir + file);
        if (!config) throw std::runtime_error("Missing " + dir + file);
        std::string line, kind, name;
        double value;
        while (std::getline(config, line)) {
            std::istringstream fields(line);
            if (fields >> kind >> name >> value && kind == "ModelParameter") p[name] = value;
        }
    }
    for (auto name : {"delSin2th_q", "delSin2th_b", "delsigma0H", "delR0l", "delR0c",
                      "delGammaWlv", "delGammaWqq"}) p.emplace(name, 0.);
    for (const auto& name : Model::FlavourWilsonCoefficient_DF2vars)
        p[name] = name.find("WCscale") == 0 ? p.at("mub") : 0.;
    p["BBsscale"] = p["BBdscale"] = p.at("mub");
    if (!withBags) {
        for (auto it = p.begin(); it != p.end();) {
            if (it->first.find("BB") == 0) it = p.erase(it);
            else ++it;
        }
    }
    return p;
}

void equal(double result, double expected, const char* quantity)
{
    if (!std::isfinite(result) || !std::isfinite(expected) ||
        std::abs(result - expected) > 2.e-10 * std::max(1., std::abs(expected))) {
        std::ostringstream error;
        error << quantity << " = " << result << ", expected " << expected;
        throw std::runtime_error(error.str());
    }
}

template<class Dilepton> void check(Dilepton& obs, double phase)
{
    equal(obs.computeAmumu(FULLNLO), std::cos(phase), "A_DeltaGamma");
    equal(obs.computeSmumu(FULLNLO), -std::sin(phase), "S_ll");
}

void run(bool withBags)
{
    Model model;
    if (!model.InitializeModel()) throw std::runtime_error("InitializeModel failed");
    model.setFlag("FlagCsi", false);
    Mll bs(model, 1, QCD::B_S, QCD::MU), bd(model, 1, QCD::B_D, QCD::MU);
    Mll bsBar(model, 2, QCD::B_S, QCD::MU), bdBar(model, 2, QCD::B_D, QCD::MU);
    Mll bsA(model, 3, QCD::B_S, QCD::MU), bdS(model, 4, QCD::B_D, QCD::MU);
    Bsmumu legacyS(model, 1);
    Bdmumu legacyD(model, 1);
    if (!model.Init(inputs(withBags))) throw std::runtime_error("Invalid test inputs");
    if (model.hasBBParameters() != withBags) throw std::runtime_error("Unexpected bag registration");
    const double promptS = bs.computeThValue(), promptD = bd.computeThValue();
    if (!(promptD > 0. && promptS > 0.)) throw std::runtime_error("Invalid prompt rates");

    auto verify = [&](double phaseD, double phaseS) {
        check(bs, phaseS); check(legacyS, phaseS);
        check(bd, phaseD); check(legacyD, phaseD);
        equal(bsA.computeThValue(), std::cos(phaseS), "factory A");
        equal(bdS.computeThValue(), -std::sin(phaseD), "factory S");
        equal(bs.computeThValue() / promptS, 1., "Bs prompt rate");
        equal(bd.computeThValue() / promptD, 1., "Bd prompt rate");
        equal(bdBar.computeThValue() / promptD, 1., "Bd integrated rate");
        const double y = model.getMesons(QCD::B_S).getDgamma_gamma() / 2.;
        equal(bsBar.computeThValue() / promptS, (1. + y * std::cos(phaseS)) / (1. - y*y),
              "Bs integrated rate");
    };
    verify(0., 0.);
    model.phiD = .12; model.phiS = .19;
    verify(-2.*model.phiD, 2.*model.phiS);

    // Rotate each full Q1 amplitude by a known angle at the low scale.
    Parameters wc;
    for (bool strange : {false, true}) {
        const double scale = model.getMub();
        auto c = strange ? model.getFlavour().ComputeCoeffBs(scale, NDR, true)
                         : model.getFlavour().ComputeCoeffBd(scale, NDR, true);
        const double angle = strange ? -.57 : .83;
        const gslpp::complex shift = (gslpp::complex(1., angle, true) - 1.)
                                     * ((*(c[LO]))(0) + (*(c[NLO]))(0));
        const std::string channel = strange ? "bs" : "bd";
        wc["reC1_" + channel] = shift.real();
        wc["imC1_" + channel] = shift.imag();
    }
    if (!model.Update(wc)) throw std::runtime_error("Q1 update failed");
    verify(.83 - 2.*model.phiD, -.57 + 2.*model.phiS);
    model.phiD = model.phiS = 0.;
    verify(.83, -.57);

    if (!model.Update({{"imC4_bs", 2.e-10}, {"reC2_bd", 1.e-11}}))
        throw std::runtime_error("Nonstandard operator update failed");
    if (withBags) {
        auto phase = [&](int channel, gslpp::complex lambda) {
            return (model.getFlavour().getDB2(channel).getM21(FULLNLO) / (lambda*lambda)).arg();
        };
        const double ps = phase(1, model.getCKM().computelamt_s());
        const double pd = phase(0, model.getCKM().computelamt_d());
        if (std::abs(ps + .57) < .1 || std::abs(pd - .83) < .1)
            throw std::runtime_error("Test operators did not change the mixing phases");
        verify(pd - 2.*model.phiD, ps + 2.*model.phiS);
    } else {
        if (!model.Update({{"DGs_Gs", 0.}})) throw std::runtime_error("Width update failed");
        equal(bsBar.computeThValue() / promptS, 1., "Bs integrated at zero width difference");
        equal(bs.computeThValue() / promptS, 1., "Bs prompt without mixing inputs");
        equal(bdBar.computeThValue() / promptD, 1., "Bd integrated without mixing inputs");
        try {
            bsA.computeThValue();
        } catch (const std::runtime_error& error) {
            if (std::string(error.what()).find("bag parameters") != std::string::npos) return;
            throw;
        }
        throw std::runtime_error("Missing matrix elements for Q4 were silently accepted");
    }
}
}

int main()
{
    run(false);
    run(true);
    std::cout << "Bd and Bs dilepton mixing tests passed.\n";
}
