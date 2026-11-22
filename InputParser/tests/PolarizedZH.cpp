#include "NPSMEFTd6.h"
#include "NPSMEFTd6General.h"
#include "ThObservable.h"
#include "ThObsFactory.h"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

namespace {
const double lambdaNP = 1000.; // GeV

template <std::size_t N>
void checkModel(NPbase& model, const std::string (&names)[N],
        const std::string& electronCurrent, const std::string& fieldStrength)
{
    std::map<std::string, double> parameters = {
        {"lambda", 0.2251}, {"A", 0.826}, {"rhob", 0.160}, {"etab", 0.350},
        {"Mz", 91.1875}, {"GF", 1.1663787e-5}, {"ale", 7.2973525698e-3}, {"dAle5Mz", 0.0275},
        {"AlsMz", 0.1185}, {"mHl", 125.09}, {"mtop", 163.0}, {"mbottom", 4.18}, {"mcharm", 1.3},
        {"mstrange", 0.0938}, {"mup", 0.0023}, {"mdown", 0.0048}, {"melectron", 5.109989e-4},
        {"mmu", 0.10565837}, {"mtau", 1.77682}, {"muc", 1.3}, {"mub", 4.8}, {"mut", 164.1}, {"muw", 80.385}
    };
    for (const auto& name : StandardModel::SMvars) parameters.emplace(name, 0.);
    for (const auto& name : names) parameters.emplace(name, 0.);
    parameters["Lambda_NP"] = lambdaNP;
    parameters[electronCurrent] = 1.;
    if (!model.InitializeModel() || !model.Init(parameters))
        throw std::runtime_error("Failed to initialize " + model.getModelName());

    ThObsFactory factory;
    auto check = [&](const std::string& name, double expected) {
        std::unique_ptr<ThObservable> observable(factory.CreateThMethod(name, model));
        const double value = observable->computeThValue();
        if (!std::isfinite(value) || std::abs(value - expected) > 1.e-8 * std::abs(expected)) {
            std::cerr << model.getModelName() << " " << name << ": expected " << expected
                      << ", got " << value << std::endl;
            throw std::runtime_error("Polarized ZH coefficient mismatch");
        }
    };

    // Isolate the right-handed electron current for a, then Q_HW for b.
    for (const std::string prefix : {"aPsk", "bPsk"}) {
        if (prefix == "bPsk" && !model.Update({{electronCurrent, 0.}, {fieldStrength, 1.}}))
            throw std::runtime_error("Failed to update " + model.getModelName());
        const double sw2 = model.getTrueSM().s02();
        const double v2OverLambda2 = std::pow(model.v() / lambdaNP, 2);
        for (int energyGeV : {250, 350, 500}) {
            const double r = std::pow(model.getMz() / energyGeV, 2);
            for (double polElectron : {-0.8, 0.8}) {
                const double polPositron = polElectron < 0. ? 0.3 : -0.3;
                const double wL = (1. - polElectron) * (1. + polPositron) / 4.;
                const double wR = (1. + polElectron) * (1. - polPositron) / 4.;
                const std::string suffix = polElectron < 0. ? "_m80_p30" : "_p80_m30";
                const double expected = prefix == "aPsk"
                        ? -wR * v2OverLambda2 / (2. * sw2 * r)
                        : 2. * (1. - sw2) * v2OverLambda2
                                * (wL * (1. + sw2 * (1. - r) / (0.5 - sw2)) + wR * r);
                check(prefix + std::to_string(energyGeV) + suffix, expected);
            }
        }
    }
}
}

int main()
{
    std::unique_ptr<NPSMEFTd6> smeft(new NPSMEFTd6());
    smeft->setModelName("NPSMEFTd6");
    checkModel(*smeft, NPSMEFTd6::NPSMEFTd6Vars, "CHe_11", "CHW");
    std::unique_ptr<NPSMEFTd6General> general(new NPSMEFTd6General());
    general->setFlag("RGEci", false);
    checkModel(*general, NPSMEFTd6General::NPSMEFTd6GeneralVars, "CHe_11r_LNP", "CHW_LNP");
}
