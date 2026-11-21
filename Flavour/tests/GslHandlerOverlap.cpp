// Run from the repository root. For this test compile StandardModel.cpp with
// -Dgsl_integration_qags=testLEP2Qags, MVgamma.cpp with
// -Dgsl_integration_cquad=testMVgammaCquad, and MVll.cpp with
// -Dgsl_integration_cquad=testMVllCquad. MVll uses the real GSL quadrature.
#include "ScopedGslErrorHandlerOff.h"
#include "StandardModel.h"
#include "MVgamma.h"
#include "MVll.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {
unsigned errors = 0, probes = 0;
bool failNextQCDFIntegral = false;
int qcdfFailure = GSL_SUCCESS;
std::unique_ptr<ScopedGslErrorHandlerOff> firstScope;

void countError(const char*, const char*, int, int) { ++errors; }
void otherHandler(const char*, const char*, int, int) { ++errors; }

void probe()
{
    // End the first scope while the production integration still owns its scope.
    firstScope.reset();
    gsl_error("handler probe", __FILE__, __LINE__, GSL_EFAILED);
    ++probes;
}

template<class Integration> void checkOverlap(const char* channel, Integration integrate, unsigned expected)
{
    auto original = gsl_set_error_handler(countError);
    errors = probes = 0;
    firstScope.reset(new ScopedGslErrorHandlerOff);
    integrate();
    firstScope.reset();
    auto restored = gsl_set_error_handler(original);
    if (probes != expected || errors != 0 || restored != countError)
        throw std::runtime_error(std::string(channel) + " lost the process-wide handler");
}

void checkPersistentHandler()
{
    auto original = gsl_set_error_handler(countError);
    errors = 0;
    {
        ScopedGslErrorHandlerOff first;
        ScopedGslErrorHandlerOff second;
        ScopedGslErrorHandlerOff::setDefaultHandler(otherHandler);
        gsl_error("handler probe", __FILE__, __LINE__, GSL_EFAILED);
    }
    auto restored = gsl_set_error_handler(original);
    if (errors != 0 || restored != otherHandler)
        throw std::runtime_error("Persistent handler change was applied at the wrong time");
}

void init(StandardModel& model, MVgamma& radiative, MVll& decay)
{
    std::map<std::string,double> p;
    const std::string dir = "Configurations/Temporary/bsll_2026/SMEFT/";
    for (auto file : {"StandardModel.conf", "Flavour.conf", "FormFactors_BSZ_Lattice.conf"}) {
        std::ifstream cfg(dir + file);
        if (!cfg) throw std::runtime_error("Missing test inputs: " + dir + file);
        std::string line, kind, name;
        double value;
        while (std::getline(cfg, line)) {
            std::istringstream row(line);
            if (row >> kind >> name >> value && kind == "ModelParameter") p[name] = value;
        }
    }
    for (auto key : {"delSin2th_q", "delSin2th_b", "delsigma0H", "delR0l", "delR0c",
                     "delGammaWlv", "delGammaWqq"}) p.emplace(key, 0.);
    for (const auto& names : {radiative.initializeMVgammaParameters(), decay.initializeMVllParameters()}) {
        model.addParameters(names);
        for (const auto& name : names) p.emplace(name, 0.);
    }
    if (!model.Init(p)) throw std::runtime_error("Invalid test inputs");
}

void checkMVllUpdates(StandardModel& model, MVll& decay)
{
    auto original = gsl_set_error_handler(countError);
    errors = 0;
    auto setLambdaB = [&](double value) {
        if (!model.Update({{"lambdaB", value}})) throw std::runtime_error("Model update failed");
    };
    auto checkHandler = [] {
        if (gsl_set_error_handler(countError) != countError || errors != 0)
            throw std::runtime_error("MVll lost or invoked the caller's GSL handler");
    };
    auto bin = [&] { return decay.integrateSigma(0, 1., 2.); };
    setLambdaB(0.35);
    const double initial = bin();
    if (!(std::isfinite(initial) && initial > 0.)) throw std::runtime_error("Invalid initial bin");
    checkHandler();
    setLambdaB(0.30);
    const double reference = bin();
    if (!(std::isfinite(reference) && reference > 0.)) throw std::runtime_error("Invalid reference bin");
    checkHandler();

    setLambdaB(0.35);
    failNextQCDFIntegral = true;
    const double failed = bin();
    checkHandler();
    if (qcdfFailure == GSL_SUCCESS || !std::isnan(failed)
            || !model.getFlavour().getUpdateFlag(QCD::B_D, QCD::K_star, QCD::MU))
        throw std::runtime_error("MVll accepted a failed QCDF integration");
    const double retry = bin();
    if (!(std::isfinite(retry) && retry > 0. && std::abs(retry / reference - 1.) > 1.e-6))
        throw std::runtime_error("MVll did not recompute the bin on retry");

    setLambdaB(0.001); // Positive stress point for exp(-x) Ei(x), outside the fit range.
    const double stress = decay.getSigma(0, 1.5);
    if (!(std::isfinite(stress) && stress > 0.)) throw std::runtime_error("QCDF overflow at small lambdaB");
    checkHandler();
    setLambdaB(0.);
    if (!std::isnan(decay.integrateDelta(0, 1., 2.))) throw std::runtime_error("Accepted zero lambdaB");
    checkHandler();
    setLambdaB(0.30);
    if (!(std::abs(bin() / reference - 1.) < 1.e-10)) throw std::runtime_error("MVll failed to recover");
    checkHandler();
    gsl_set_error_handler(original);
}
}

extern "C" int testLEP2Qags(const gsl_function*, double, double, double, double, size_t,
                            gsl_integration_workspace*, double* result, double* err)
{
    probe();
    *result = *err = 0.;
    return GSL_SUCCESS;
}

extern "C" int testMVgammaCquad(const gsl_function*, double, double, double, double,
                                gsl_integration_cquad_workspace*, double* result, double* err, size_t*)
{
    probe();
    *result = *err = 0.;
    return GSL_SUCCESS;
}

extern "C" int testMVllCquad(const gsl_function* f, double a, double b, double epsabs, double epsrel,
                             gsl_integration_cquad_workspace* w, double* result, double* err, size_t* nevals)
{
    if (failNextQCDFIntegral && a == 0. && b == 1.) {
        failNextQCDFIntegral = false;
        qcdfFailure = gsl_integration_cquad(f, a, b, -1., 0., w, result, err, nevals);
        return qcdfFailure;
    }
    return gsl_integration_cquad(f, a, b, epsabs, epsrel, w, result, err, nevals);
}

int main()
{
    StandardModel model;
    if (!model.InitializeModel()) throw std::runtime_error("Could not initialize model");
    model.setFlag("NoApproximateGammaZ", true);
    model.setFlag("LEP2WeakBox", false);
    model.setFlag("FlagCsi", false);
    model.setFlag("UsezExpansion", false);
    auto& radiative = model.getFlavour().getMVgamma(QCD::B_D, QCD::K_star);
    auto& decay = model.getFlavour().getMVll(QCD::B_D, QCD::K_star, QCD::MU);
    init(model, radiative, decay);
    checkOverlap("LEP2", [&] { model.LEP2sigmaMu(200.*200.); }, 1);
    checkOverlap("MVgamma", [&] { radiative.updateParameters(); }, 4);
    checkPersistentHandler();
    checkMVllUpdates(model, decay);
    std::cout << "GSL handler overlap, MVll repeated updates and failure recovery passed.\n";
}
