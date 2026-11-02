#include "WarmupConfig.h"
#include <initializer_list>
#include <iostream>

static void check(bool condition)
{
    if (!condition) throw std::runtime_error("Warmup config assertion failed");
}

static void reject(std::initializer_list<const char*> lines, unsigned budget = 12,
                   unsigned factorized = 3, bool joint = true, bool load = false,
                   bool read = false, bool optimizer = false)
{
    bool failed = false;
    try {
        HEPfitWarmupConfig config;
        for (auto line : lines) config.ParseLine(line);
        config.ValidateRun(budget, factorized, joint, load, read, optimizer);
    } catch (const std::runtime_error&) { failed = true; }
    check(failed);
}

int main()
{
    HEPfitWarmupConfig defaults;
    defaults.Validate();
    check(!defaults.Enabled() && !defaults.staged && !defaults.smooth);
    check(defaults.settlingSweeps == 1000);
    check(!defaults.ParseLine("PrerunMaxIter 12"));
    HEPfitWarmupConfig config;
    check(config.ParseLine(" SettlingSweeps 7 # discard"));
    check(config.ParseLine("SmoothAdaptation true\r"));
    check(config.ParseLine("PreRunPolicy staged"));
    config.ValidateRun(12, 3, true, false, false, false);
    check(config.staged && config.smooth && config.settlingSweeps == 7);
    config.ParseLine("SmoothAdaptation false");
    check(config.staged && !config.smooth);
    config.ParseLine("SettlingSweeps 0");
    config.Validate();
    check(config.settlingSweeps == 0);
    HEPfitWarmupConfig independent;
    independent.ParseLine("SmoothAdaptation true");
    independent.ValidateRun(12, 3, true, false, false, false);
    check(!independent.staged && independent.smooth);
    reject({"PreRunPolicy unknown"});
    reject({"PreRunPolicy"});
    reject({"PreRunPolicy staged junk"});
    reject({"SmoothAdaptation yes"});
    reject({"SmoothAdaptation"});
    reject({"SmoothAdaptation true junk"});
    reject({"SettlingSweeps"});
    reject({"SettlingSweeps 2 junk"});
    reject({"SettlingSweeps 1"});
    reject({"PreRunPolicy staged", "SettlingSweeps -1"});
    reject({"PreRunPolicy staged", "SettlingSweeps 1.5"});
    reject({"PreRunPolicy staged", "SettlingSweeps 2x"});
    reject({"PreRunPolicy staged", "SettlingSweeps 2147483648"});
    reject({"PreRunPolicy staged", "SettlingSweeps 999999999999999999999999"});
    reject({"PreRunPolicy staged"}, 0);
    reject({"PreRunPolicy staged", "SettlingSweeps 1"}, 2147483647U);
    reject({"PreRunPolicy staged"}, 12, 12);
    reject({"SmoothAdaptation true"}, 12, 12);
    reject({"PreRunPolicy staged"}, 12, 3, true, true);
    reject({"SmoothAdaptation true"}, 12, 3, true, false, true);
    reject({"PreRunPolicy staged"}, 12, 3, true, false, false, true);
    std::cout << "Warmup defaults, independent controls and configuration validation passed\n";
}
