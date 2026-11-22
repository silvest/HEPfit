#include "NPSMEFTd6General.h"
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
std::map<std::string, double> inputs(bool mwInput)
{
    std::map<std::string, double> p = {
        {"lambda", 0.2251}, {"A", 0.826}, {"rhob", 0.160}, {"etab", 0.350},
        {"Mz", 91.1875}, {"GF", 1.1663787e-5}, {"ale", 7.2973525698e-3},
        {"AlsMz", 0.1185}, {"mHl", 125.09}, {"mtop", 163.0}, {"mbottom", 4.18},
        {"mcharm", 1.3}, {"mstrange", 0.0938}, {"mup", 0.0023}, {"mdown", 0.0048},
        {"melectron", 5.109989e-4}, {"mmu", 0.10565837}, {"mtau", 1.77682},
        {"muc", 1.3}, {"mub", 4.8}, {"mut", 164.1}, {"muw", 80.385}, {"Lambda_NP", 1000.}
    };
    for (const auto& name : StandardModel::SMvars) p.emplace(name, 0.);
    for (const auto& name : NPSMEFTd6General::NPSMEFTd6GeneralVars) p.emplace(name, 0.);
    if (mwInput) {
        p.erase("dAle5Mz");
        p["Mw_inp"] = 80.379;
    } else {
        p["dAle5Mz"] = 0.0275;
    }
    return p;
}

// Vary v^2 and solve the tree-level input relations before evaluating the full loop amplitude.
double referenceWidth(const NPSMEFTd6General& sm, bool mwInput, double deltaGF)
{
    const double v2Ratio = 1. + deltaGF;
    const double c2 = std::pow(sm.getMw_tree() / sm.getMz(), 2);
    const double shiftedC2 = (1. + std::sqrt(1. - 4. * c2 * (1. - c2) * v2Ratio)) / 2.;
    const double mw = mwInput ? sm.getMw_tree() : sm.getMz() * std::sqrt(shiftedC2);
    const double alphaRatio = mwInput ? 1. / v2Ratio : 1.;
    auto tau = [&](double mass) { return std::pow(2. * mass / sm.getMHl(), 2); };
    gslpp::complex amplitude = sm.AH_W(tau(mw));
    for (auto q : {QCD::TOP, QCD::BOTTOM, QCD::CHARM, QCD::STRANGE}) {
        const Particle& quark = sm.getQuarks(q);
        const double mass = q == QCD::TOP ? sm.getMtpole() : quark.getMass();
        amplitude += 3. * std::pow(quark.getCharge(), 2) * sm.AH_f(tau(mass));
    }
    for (auto l : {QCD::TAU, QCD::MU})
        amplitude += sm.AH_f(tau(sm.getLeptons(l).getMass()));
    return alphaRatio * alphaRatio / v2Ratio * amplitude.abs2();
}

void checkScheme(bool mwInput)
{
    std::unique_ptr<NPSMEFTd6General> model(new NPSMEFTd6General());
    auto p = inputs(mwInput);
    if (!model->setFlag("MWinput", mwInput) || !model->setFlag("RGEci", false)
            || !model->InitializeModel() || !model->Init(p))
        throw std::runtime_error("Failed to initialize diphoton test");

    for (double mh : {125.09, 145., 180.}) {
        if (!model->Update({{"mHl", mh}, {"CHl3_11r_LNP", 0.}}))
            throw std::runtime_error("Failed to set the SM reference point");
        // Use the same positive and negative perturbations for the reference and the public observable.
        const double step = 1.e-5;
        const double reference = (referenceWidth(*model, mwInput, step) - referenceWidth(*model, mwInput, -step))
                / (2. * step * referenceWidth(*model, mwInput, 0.));
        const double wc = step * std::pow(p.at("Lambda_NP") / model->v(), 2);
        auto width = [&](double coefficient) {
            if (!model->Update({{"CHl3_11r_LNP", coefficient}}))
                throw std::runtime_error("Failed to update the electron-current coefficient");
            return model->GammaHgagaRatio();
        };
        const double plus = width(wc);
        const double minus = width(-wc);
        const double slope = (plus - minus) / (2. * step);
        if (!std::isfinite(slope) || !std::isfinite(reference) || std::abs(slope - reference) > 1.e-7) {
            std::cerr << "MWinput=" << mwInput << ", mH=" << mh
                      << ": expected " << reference << ", got " << slope << std::endl;
            throw std::runtime_error("Diphoton width slope mismatch");
        }
    }
}
}

int main()
{
    checkScheme(true);
    checkScheme(false);
}
