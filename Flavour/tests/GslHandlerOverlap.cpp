// Run from the repository root. For this test compile StandardModel.cpp with
// -Dgsl_integration_qags=testLEP2Qags, MVgamma.cpp with
// -Dgsl_integration_cquad=testMVgammaCquad, and MVll.cpp with
// -Dgsl_integration_cquad=testMVllCquad. Compile MPll.cpp with
// -Dgsl_integration_cquad=testMPllCquad. Test both settings of SPLINE.
#include "ScopedGslErrorHandlerOff.h"
#include "StandardModel.h"
#include "MVgamma.h"
#include "MVll.h"
#include "MPll.h"
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
bool failNextMPllIntegral = false;
int mpllFailure = GSL_SUCCESS;
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

void init(StandardModel& model, MVgamma& radiative, MVll& decay, MPll& pseudoscalar)
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
    for (const auto& names : {radiative.initializeMVgammaParameters(), decay.initializeMVllParameters(),
                              pseudoscalar.initializeMPllParameters()}) {
        model.addParameters(names);
        for (const auto& name : names) p.emplace(name, 0.);
    }
    if (!model.Init(p)) throw std::runtime_error("Invalid test inputs");
}

void checkDirectAmplitudes(StandardModel& model, MVll& decay)
{
    const double first = decay.H_V_0(1.5, false).abs();
    if (!(std::isfinite(first) && first > 0.)) throw std::runtime_error("Invalid direct amplitude on first use");
    for (auto amplitude : {&MVll::H_V_0, &MVll::H_V_p, &MVll::H_V_m}) {
        // Each helicity must handle a pending update without a prior bin or getter.
        if (!model.Update({{"lambdaB", 0.35}})) throw std::runtime_error("Model update failed");
        failNextQCDFIntegral = true;
        qcdfFailure = GSL_SUCCESS;
        if (!std::isnan((decay.*amplitude)(1.5, false).real()) || qcdfFailure == GSL_SUCCESS)
            throw std::runtime_error("Direct amplitude skipped its pending QCDF update");
        auto retry = (decay.*amplitude)(1.5, false);
        if (!(std::isfinite(retry.abs()) && retry.abs() > 0.)
                || model.getFlavour().getUpdateFlag(QCD::B_D, QCD::K_star, QCD::MU))
            throw std::runtime_error("Direct amplitude failed to retry at the same model point");
    }
}

void checkMPllUpdates(StandardModel& model, MPll& decay)
{
    auto original = gsl_set_error_handler(countError);
    errors = 0;
    const double bin = decay.integrateSigma(0, 1., 2.);
    if (!(std::isfinite(bin) && bin > 0.) || gsl_set_error_handler(countError) != countError)
        throw std::runtime_error("Invalid initial MPll bin or handler");
    const double mb = model.getQuarks(QCD::BOTTOM).getMass();
    if (!model.Update({{"mbottom", mb + 0.05}})) throw std::runtime_error("Model update failed");
#if !SPLINE
    failNextMPllIntegral = true;
#endif
    decay.getSigma(0, 1.5);
    if (gsl_set_error_handler(countError) != countError || errors != 0)
        throw std::runtime_error("MPll preparation invoked or lost the caller's handler");
#if !SPLINE
    if (mpllFailure == GSL_SUCCESS || !model.getFlavour().getUpdateFlag(QCD::B_P, QCD::K_P, QCD::MU))
        throw std::runtime_error("MPll accepted a failed fit sample");
    const double retry = decay.getSigma(0, 1.5);
    if (!(std::isfinite(retry) && retry > 0.) || model.getFlavour().getUpdateFlag(QCD::B_P, QCD::K_P, QCD::MU))
        throw std::runtime_error("MPll failed to prepare the fit on retry");
#endif
    const double lambdaB = model.getMesons(QCD::B_P).getLambdaM();
    if (!model.Update({{"lambdaB", 0.001}, {"mbottom", mb + 0.02}})) throw std::runtime_error("Model update failed");
    const double stress = decay.getSigma(0, 1.5);
    if (!(std::isfinite(stress) && stress > 0.)) throw std::runtime_error("MPll overflow at small lambdaB");
    if (!model.Update({{"lambdaB", lambdaB}, {"mbottom", mb}})) throw std::runtime_error("Model update failed");
    if (!(std::abs(decay.integrateSigma(0, 1., 2.) / bin - 1.) < 1.e-10))
        throw std::runtime_error("MPll failed to recover the reference bin");
    if (gsl_set_error_handler(original) != countError || errors != 0)
        throw std::runtime_error("MPll failed to restore the caller's handler");
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

extern "C" int testMPllCquad(const gsl_function* f, double a, double b, double epsabs, double epsrel,
                             gsl_integration_cquad_workspace* w, double* result, double* err, size_t* nevals)
{
    if (failNextMPllIntegral && a == 0. && b == 1.) {
        failNextMPllIntegral = false;
        mpllFailure = gsl_integration_cquad(f, a, b, -1., 0., w, result, err, nevals);
        return mpllFailure;
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
    model.setFlag("MPll_FNALMILC_flag", true);
    auto& pseudoscalar = model.getFlavour().getMPll(QCD::B_P, QCD::K_P, QCD::MU);
    init(model, radiative, decay, pseudoscalar);
    checkOverlap("LEP2", [&] { model.LEP2sigmaMu(200.*200.); }, 1);
    checkOverlap("MVgamma", [&] { radiative.updateParameters(); }, 4);
    checkPersistentHandler();
    checkDirectAmplitudes(model, decay);
    checkMVllUpdates(model, decay);
    checkMPllUpdates(model, pseudoscalar);
    std::cout << "GSL handler overlap, direct amplitudes and repeated updates passed.\n";
}
