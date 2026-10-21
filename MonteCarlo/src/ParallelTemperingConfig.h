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

    /// Return false for a normal HEPfit option, true for a consumed PT option.
    bool ParseLine(std::string line)
    {
        line = line.substr(0, line.find('#'));
        std::istringstream input(line);
        std::string key, value, extra;
        input >> key;
        if (key != "ParallelTempering" && key != "PTTemperatures" &&
                key != "PTSwapInterval" && key != "PTDiagnostics") return false;
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
            if (key == "PTSwapInterval") {
                if (value.find_first_not_of("0123456789") != std::string::npos)
                    throw std::runtime_error("PTSwapInterval must be a positive integer.");
                unsigned long interval;
                try { interval = std::stoul(value); }
                catch (const std::exception&) { throw std::runtime_error("Invalid PTSwapInterval."); }
                if (!interval || interval > std::numeric_limits<unsigned>::max())
                    throw std::runtime_error("Invalid PTSwapInterval.");
                swapInterval = static_cast<unsigned>(interval);
                hasSwapInterval = true;
            } else {
                if (value != "true" && value != "false")
                    throw std::runtime_error(key + " must be true or false.");
                if (key == "ParallelTempering") enabled = value == "true";
                else diagnostics = value == "true";
            }
        }
        return true;
    }

    void Validate() const
    {
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
};
#endif
