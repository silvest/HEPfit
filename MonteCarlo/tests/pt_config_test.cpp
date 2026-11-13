#include "ParallelTemperingConfig.h"
#include <iostream>
#include <initializer_list>

static void check(bool condition) { if (!condition) throw std::runtime_error("PT config assertion failed"); }
static void reject(std::initializer_list<const char*> lines)
{
    bool failed = false;
    try {
        HEPfitPTConfig config;
        for (auto line : lines) config.ParseLine(line);
        config.Validate();
    } catch (const std::runtime_error&) { failed = true; }
    check(failed);
}
int main()
{
    HEPfitPTConfig config;
    check(!config.ParseLine("NChains 3"));
    check(config.ParseLine(" PTTemperatures 1:2:4,8 16 # ladder"));
    config.ParseLine("PTSwapInterval 7");
    config.ParseLine("PTDiagnostics false");
    config.ParseLine("ParallelTempering true");
    config.Validate();
    check(config.temperatures.size() == 5 && config.temperatures.back() == 16.);
    check(config.enabled && !config.diagnostics && config.swapInterval == 7);
    HEPfitPTConfig().Validate();
    reject({"ParallelTempering true"});
    reject({"ParallelTempering true", "PTTemperatures 1"});
    reject({"ParallelTempering true", "PTTemperatures 2 4"});
    reject({"ParallelTempering true", "PTTemperatures 1 4 2"});
    reject({"ParallelTempering true", "PTTemperatures 1 2 2"});
    reject({"PTTemperatures 1 2"});
    reject({"PTSwapInterval 10"});
    reject({"PTSwapInterval -1"});
    reject({"PTSwapInterval 0"});
    reject({"PTSwapInterval 1.5"});
    reject({"PTSwapInterval 9999999999999999999999999"});
    reject({"PTTemperatures 1 nan"});
    reject({"PTTemperatures 1 inf"});
    reject({"PTTemperatures 1 2junk"});
    reject({"PTTemperatures"});
    reject({"ParallelTempering yes"});
    reject({"PTDiagnostics true junk"});
    reject({"PTDiagnostics"});
    reject({"PTAdaptTemperatures yes"});
    reject({"PTAdaptTemperatures true"});
    reject({"ParallelTempering true", "PTTemperatures 1:2", "PTAdaptTemperatures true"});
    reject({"PTAdaptInterval 1"});
    reject({"PTAdaptInterval -1"});
    reject({"PTAdaptInterval 2147483648"});
    reject({"PTAdaptRounds 31"});
    reject({"PTAdaptRounds 1.5"});
    reject({"PTAdaptDamping 0"});
    reject({"PTAdaptDamping nan"});
    reject({"PTAdaptDamping 1.01"});
    config.ParseLine("PTAdaptTemperatures true");
    config.ParseLine("PTAdaptInterval 20");
    config.ParseLine("PTAdaptRounds 3");
    config.ParseLine("PTAdaptDamping .75");
    config.ValidateAdaptation(true, false, 5, 200, 150, 20, true);
    for (unsigned problem=0; problem<5; ++problem) {
        bool failed=false;
        try { config.ValidateAdaptation(problem!=0, problem==1, problem==2?0:5,
                                         problem==3?159:200, problem==4?190:150, 20, true); }
        catch (const std::runtime_error&) { failed=true; }
        check(failed);
    }
    std::cout << "PT configuration validation passed\n";
}
