#ifndef HEPFIT_PROPOSAL_CONFIG_H
#define HEPFIT_PROPOSAL_CONFIG_H

#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

/// Full BAT proposal caches, distinct from legacy HEPfit prerun data.
struct HEPfitProposalConfig {
    std::string save;
    std::string load;
    unsigned warmup = 1000;
    bool hasWarmup = false;

    bool ParseLine(const std::string& line)
    {
        std::istringstream input(line);
        std::string key, value;
        input >> key;
        if (key != "SaveProposal" && key != "LoadProposal" && key != "ProposalWarmup")
            return false;
        input >> std::ws;
        if (input.peek() == '#' || !(input >> std::quoted(value)) || value.empty())
            throw std::runtime_error(key + " requires exactly one non-empty value.");
        input >> std::ws;
        if (!input.eof() && input.peek() != '#')
            throw std::runtime_error(key + " requires exactly one value (quote paths containing spaces).");
        if (key == "SaveProposal") save = value;
        else if (key == "LoadProposal") load = value;
        else {
            if (value.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("ProposalWarmup must be a non-negative integer no larger than INT_MAX.");
            unsigned long count;
            try { count = std::stoul(value); }
            catch (const std::exception&) { throw std::runtime_error("Invalid ProposalWarmup."); }
            if (count > static_cast<unsigned long>(std::numeric_limits<int>::max()))
                throw std::runtime_error("ProposalWarmup must be no larger than INT_MAX.");
            warmup = static_cast<unsigned>(count);
            hasWarmup = true;
        }
        return true;
    }

    bool Enabled() const { return !load.empty() || !save.empty(); }

    void ValidateRun(bool staged, bool smooth, bool restoreRun,
                     bool readPreRunData, bool optimizerOnly) const
    {
        if (hasWarmup && load.empty())
            throw std::runtime_error("ProposalWarmup requires LoadProposal.");
        if (!Enabled()) return;
        if (restoreRun || readPreRunData)
            throw std::runtime_error("SaveProposal/LoadProposal cannot be combined with --loadRun or ReadPreRunData.");
        if (optimizerOnly)
            throw std::runtime_error("SaveProposal/LoadProposal require MCMC, not RunMinuitOnly or ComputeHessianOnly.");
        if (!load.empty() && (staged || smooth))
            throw std::runtime_error("LoadProposal skips adaptation: remove PreRunPolicy staged and SmoothAdaptation true; use ProposalWarmup for discarded fixed-proposal sweeps.");
    }
};
#endif
