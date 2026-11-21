// Run from the repository root. For this test compile StandardModel.cpp with
// -Dgsl_integration_qags=testLEP2Qags and MVgamma.cpp with
// -Dgsl_integration_cquad=testMVgammaCquad. Only the quadrature is replaced.
#include "ScopedGslErrorHandlerOff.h"
#include "StandardModel.h"
#include "MVgamma.h"
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace {
unsigned errors = 0, probes = 0;
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

void init(StandardModel& model, MVgamma& radiative)
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
    const auto names = radiative.initializeMVgammaParameters();
    model.addParameters(names);
    for (const auto& name : names) p.emplace(name, 0.);
    if (!model.Init(p)) throw std::runtime_error("Invalid test inputs");
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

int main()
{
    StandardModel model;
    if (!model.InitializeModel()) throw std::runtime_error("Could not initialize model");
    model.setFlag("NoApproximateGammaZ", true);
    model.setFlag("LEP2WeakBox", false);
    model.setFlag("FlagCsi", false);
    model.setFlag("UsezExpansion", false);
    auto& radiative = model.getFlavour().getMVgamma(QCD::B_D, QCD::K_star);
    init(model, radiative);
    checkOverlap("LEP2", [&] { model.LEP2sigmaMu(200.*200.); }, 1);
    checkOverlap("MVgamma", [&] { radiative.updateParameters(); }, 4);
    checkPersistentHandler();
    std::cout << "LEP2 and MVgamma overlap and handler restoration passed.\n";
}
