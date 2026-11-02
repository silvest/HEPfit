#ifndef HEPFIT_WARMUP_CONFIG_H
#define HEPFIT_WARMUP_CONFIG_H

#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

/// Order-independent warmup options, independent of the installed BAT API.
struct HEPfitWarmupConfig {
    bool staged = false;
    bool smooth = false;
    unsigned settlingSweeps = 1000;
    bool hasSettlingSweeps = false;

    bool ParseLine(std::string line)
    {
        line = line.substr(0, line.find('#'));
        std::istringstream input(line);
        std::string key, value, extra;
        input >> key;
        if (key != "PreRunPolicy" && key != "SettlingSweeps" && key != "SmoothAdaptation")
            return false;
        if (!(input >> value) || (input >> extra))
            throw std::runtime_error(key + " requires exactly one value.");
        if (key == "PreRunPolicy") {
            if (value != "legacy" && value != "staged")
                throw std::runtime_error("PreRunPolicy must be legacy or staged.");
            staged = value == "staged";
        } else if (key == "SmoothAdaptation") {
            if (value != "true" && value != "false")
                throw std::runtime_error("SmoothAdaptation must be true or false.");
            smooth = value == "true";
        } else {
            if (value.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("SettlingSweeps must be a non-negative integer no larger than INT_MAX.");
            unsigned long count;
            try { count = std::stoul(value); }
            catch (const std::exception&) { throw std::runtime_error("Invalid SettlingSweeps."); }
            if (count > static_cast<unsigned long>(std::numeric_limits<int>::max()))
                throw std::runtime_error("SettlingSweeps must be no larger than INT_MAX.");
            settlingSweeps = static_cast<unsigned>(count);
            hasSettlingSweeps = true;
        }
        return true;
    }

    bool Enabled() const { return staged || smooth; }

    void Validate() const
    {
        if (hasSettlingSweeps && !staged)
            throw std::runtime_error("SettlingSweeps requires PreRunPolicy staged.");
    }

    void ValidateRun(unsigned adaptationSweeps, unsigned factorizedSweeps,
                     bool multivariate, bool restoreRun, bool readPreRunData,
                     bool optimizerOnly) const
    {
        Validate();
        if (!Enabled()) return;
        if (restoreRun || readPreRunData)
            throw std::runtime_error("PreRunPolicy staged/SmoothAdaptation true cannot be combined with --loadRun or ReadPreRunData.");
        if (optimizerOnly)
            throw std::runtime_error("PreRunPolicy staged/SmoothAdaptation true require MCMC, not RunMinuitOnly or ComputeHessianOnly.");
        const unsigned limit = static_cast<unsigned>(std::numeric_limits<int>::max());
        if (!adaptationSweeps || adaptationSweeps > limit)
            throw std::runtime_error("Warmup requires a positive PrerunMaxIter no larger than INT_MAX.");
        if (staged && settlingSweeps > limit - adaptationSweeps)
            throw std::runtime_error("Staged adaptation plus settling must fit in INT_MAX sweeps.");
        if (multivariate && factorizedSweeps >= adaptationSweeps)
            throw std::runtime_error("Multivariate warmup requires NIterationsPreRunFactorized < PrerunMaxIter.");
    }
};
#endif
