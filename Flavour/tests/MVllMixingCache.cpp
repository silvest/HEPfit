// Link this test with MVll.cpp compiled with -Dgsl_integration_cquad=countMVllCquad.
// Run from the repository root. The wrapper counts only the outer bin integrals.
#include "StandardModel.h"
#include "MVll.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
const std::array<std::pair<double, double>, 2> bins = {{{2.13, 2.73}, {4.31, 4.71}}};
unsigned integrations = 0;
using Parameters = std::map<std::string, double>;

Parameters inputParameters()
{
    Parameters p;
    const std::string dir = "Configurations/Temporary/bsll_2024/config_LHCb_BKstmumu/";
    for (auto name : {"StandardModel.conf", "Flavour.conf", "FormFactors.conf"}) {
        std::ifstream config(dir + name);
        if (!config) throw std::runtime_error("Missing test input " + dir + name);
        std::string line, kind, key;
        double value;
        while (std::getline(config, line)) {
            std::istringstream fields(line);
            if (fields >> kind >> key >> value && kind == "ModelParameter") p[key] = value;
        }
    }
    for (auto key : {"delSin2th_q", "delSin2th_b", "delsigma0H", "delR0l", "delR0c",
                     "delGammaWlv", "delGammaWqq"}) p.emplace(key, 0.);
    return p;
}

enum class Change { cold, unchanged, mass, width, widthSign, phiFormFactor };
bool mixesInDelta(int i) { return i == 6 || i == 7 || i == 8 || i == 10 || i == 11; }

struct Point {
    StandardModel model;
    MVll& channel(bool phi)
    {
        return model.getFlavour().getMVll(phi ? QCD::B_S : QCD::B_D,
                                         phi ? QCD::PHI : QCD::K_star, QCD::MU);
    }
    explicit Point(Parameters p)
    {
        model.setModelName("StandardModel");
        if (!model.InitializeModel()) throw std::runtime_error("Model initialization failed");
        for (bool phi : {false, true}) {
            const auto names = channel(phi).initializeMVllParameters();
            model.addParameters(names);
            for (const auto& key : names) p.emplace(key, 0.);
        }
        if (!model.Init(p)) throw std::runtime_error("Invalid test parameters");
    }
};

void check(Point& warm, Point& fresh, Change change)
{
    for (const auto& bin : bins) for (bool phi : {false, true}) {
        for (bool delta : {false, true}) for (int i = 0; i < 12; ++i) {
            if (delta && (i == 4 || i == 5 || i == 9)) continue;
            auto integral = [&](Point& point) {
                auto& decay = point.channel(phi);
                return delta ? decay.integrateDelta(i, bin.first, bin.second)
                             : decay.integrateSigma(i, bin.first, bin.second);
            };
            const bool mixingTerm = delta ? mixesInDelta(i) : !mixesInDelta(i);
            const bool mixingChanged = delta ? change == Change::mass || change == Change::width
                                             : change == Change::width || change == Change::widthSign;
            const bool mustIntegrate = change == Change::cold ||
                (phi && (change == Change::phiFormFactor || (mixingTerm && mixingChanged)));
            const unsigned before = integrations;
            const double actual = integral(warm);
            const unsigned calls = integrations - before;
            const double expected = integral(fresh);
            if (calls != unsigned(mustIntegrate) || !std::isfinite(actual) || !std::isfinite(expected) ||
                std::abs(actual - expected) > 1.e-30 + 2.e-10 * std::abs(expected)) {
                std::ostringstream message;
                message << (phi ? "phi " : "K* ") << (delta ? "Delta " : "Sigma ") << i
                        << " at " << bin.first << ", step " << int(change) << ": " << calls
                        << " integrations, expected " << mustIntegrate
                        << ", values " << actual << " and " << expected;
                throw std::runtime_error(message.str());
            }
        }
    }
}
}

extern "C" int countMVllCquad(const gsl_function* f, double a, double b, double absolute,
                              double relative, gsl_integration_cquad_workspace* workspace,
                              double* result, double* error, size_t* evaluations)
{
    if (std::find(bins.begin(), bins.end(), std::make_pair(a,b)) != bins.end()) ++integrations;
    return gsl_integration_cquad(f, a, b, absolute, relative, workspace, result, error, evaluations);
}

int main()
{
    auto p = inputParameters();
    p["xs_phi"] = 25.;
    p["DGs_Gs"] = .12;
    p["imh_m_1"] = .00002;
    Point warm(p);
    auto step = [&](Change change, const Parameters& update) {
        for (const auto& v : update) p[v.first] = v.second;
        if (change != Change::cold && !warm.model.Update(update)) throw std::runtime_error("Update failed");
        Point fresh(p);
        check(warm, fresh, change);
    };
    step(Change::cold, {});
    step(Change::unchanged, {{"xs_phi", 25.}});
    step(Change::unchanged, {{"delGammaWqq", .01}});
    step(Change::mass, {{"xs_phi", 17.}});
    step(Change::width, {{"DGs_Gs", .20}});
    step(Change::widthSign, {{"DGs_Gs", -.20}});
    step(Change::width, {{"DGs_Gs", .12}});
    step(Change::mass, {{"xs_phi", 25.}});
    const double referenceV = p.at("a_0Vphi");
    step(Change::phiFormFactor, {{"a_0Vphi", 1.05 * referenceV}});
    step(Change::phiFormFactor, {{"a_0Vphi", referenceV}});
    step(Change::unchanged, {});
    std::cout << "Angular bins agree with fresh objects and have the expected integration counts.\n";
}
