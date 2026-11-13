#ifndef HEPFIT_PARALLEL_TEMPERING_CONFIG_H
#define HEPFIT_PARALLEL_TEMPERING_CONFIG_H

#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/// Order-independent PT options, shared by the runtime parser and regression tests.
struct HEPfitPTConfig {
    bool enabled = false;
    bool diagnostics = true;
    unsigned swapInterval = 10;
    std::vector<double> temperatures;
    bool hasSwapInterval = false;
    bool adaptTemperatures = false, hasAdaptOptions = false;
    unsigned adaptInterval = 1000, adaptRounds = 4;
    double adaptDamping = 0.5;

    /// Return false for a normal HEPfit option, true for a consumed PT option.
    bool ParseLine(std::string line)
    {
        line = line.substr(0, line.find('#'));
        std::istringstream input(line);
        std::string key, value, extra;
        input >> key;
        if (key != "ParallelTempering" && key != "PTTemperatures" &&
                key != "PTSwapInterval" && key != "PTDiagnostics" && key != "PTAdaptTemperatures" &&
                key != "PTAdaptInterval" && key != "PTAdaptRounds" && key != "PTAdaptDamping") return false;
        if (key == "PTTemperatures") {
            std::getline(input, line);
            for (char& c : line) if (c == ',' || c == ':') c = ' ';
            std::istringstream values(line);
            temperatures.clear();
            while (values >> value) {
                std::size_t used = 0;
                double temperature;
                try { temperature = std::stod(value, &used); }
                catch (const std::exception&) { throw std::runtime_error("Invalid PTTemperatures value: " + value); }
                if (used != value.size() || !std::isfinite(temperature))
                    throw std::runtime_error("Invalid PTTemperatures value: " + value);
                temperatures.push_back(temperature);
            }
            if (temperatures.empty()) throw std::runtime_error("PTTemperatures needs a temperature ladder.");
        } else {
            if (!(input >> value) || (input >> extra))
                throw std::runtime_error(key + " requires exactly one value.");
            if (key == "PTAdaptDamping") {
                std::size_t used = 0;
                try { adaptDamping = std::stod(value, &used); }
                catch (const std::exception&) { throw std::runtime_error("Invalid PTAdaptDamping."); }
                if (used != value.size() || !std::isfinite(adaptDamping) || adaptDamping <= 0. || adaptDamping > 1.)
                    throw std::runtime_error("PTAdaptDamping must be in (0,1].");
                hasAdaptOptions = true;
            } else if (key == "PTSwapInterval" || key == "PTAdaptInterval" || key == "PTAdaptRounds") {
                if (value.find_first_not_of("0123456789") != std::string::npos)
                    throw std::runtime_error(key + " must be a positive integer.");
                unsigned long interval;
                try { interval = std::stoul(value); }
                catch (const std::exception&) { throw std::runtime_error("Invalid " + key); }
                const unsigned long limit = key == "PTAdaptRounds" ? 30 :
                    (key == "PTAdaptInterval" ? std::numeric_limits<int>::max() : std::numeric_limits<unsigned>::max());
                if (!interval || interval > limit) throw std::runtime_error("Invalid " + key);
                if (key == "PTSwapInterval") { swapInterval = static_cast<unsigned>(interval); hasSwapInterval = true; }
                else { (key == "PTAdaptInterval" ? adaptInterval : adaptRounds) = static_cast<unsigned>(interval); hasAdaptOptions = true; }
            } else {
                if (value != "true" && value != "false")
                    throw std::runtime_error(key + " must be true or false.");
                if (key == "ParallelTempering") enabled = value == "true";
                else if (key == "PTAdaptTemperatures") adaptTemperatures = value == "true";
                else diagnostics = value == "true";
            }
        }
        return true;
    }

    void Validate() const
    {
        if (hasAdaptOptions && !adaptTemperatures)
            throw std::runtime_error("PTAdaptInterval/PTAdaptRounds/PTAdaptDamping require PTAdaptTemperatures true.");
        if (adaptTemperatures && (!enabled || temperatures.size() < 3))
            throw std::runtime_error("PTAdaptTemperatures requires ParallelTempering true and at least three temperatures.");
        if (!enabled) {
            if (!temperatures.empty() || hasSwapInterval)
                throw std::runtime_error("PTTemperatures/PTSwapInterval require ParallelTempering true.");
            return;
        }
        if (temperatures.size() < 2 || temperatures.front() != 1.)
            throw std::runtime_error("PTTemperatures must contain at least two temperatures, starting at 1.");
        for (std::size_t i = 0; i < temperatures.size(); ++i)
            if (!std::isfinite(temperatures[i]) || (i && temperatures[i] <= temperatures[i-1]))
                throw std::runtime_error("PTTemperatures must be finite and strictly increasing from 1.");
    }

    void ValidateAdaptation(bool staged, bool loadProposal, unsigned settling, unsigned budget,
                            unsigned factorized, unsigned check, bool joint) const
    {
        Validate();
        if (!adaptTemperatures) return;
        if (!staged || loadProposal || !settling)
            throw std::runtime_error("PTAdaptTemperatures requires PreRunPolicy staged, positive SettlingSweeps and no LoadProposal.");
        const unsigned long long total = static_cast<unsigned long long>(adaptInterval)*((1ULL << adaptRounds)-1);
        const unsigned long long last = total > (joint ? factorized : 0) ? total : (joint ? factorized : 0);
        if (!check || last+check > budget || adaptInterval < swapInterval)
            throw std::runtime_error("PTAdaptTemperatures batches must include exchanges and leave a full proposal-check window on the final ladder within PrerunMaxIter.");
    }
};
#endif
