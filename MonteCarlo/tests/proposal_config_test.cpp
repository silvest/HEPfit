#include "ProposalConfig.h"
#include <initializer_list>
#include <iostream>

static void check(bool condition)
{
    if (!condition) throw std::runtime_error("Proposal config assertion failed");
}

static void reject(std::initializer_list<const char*> lines, bool staged = false,
                   bool smooth = false, bool loadRun = false, bool read = false,
                   bool optimizer = false)
{
    bool failed = false;
    try {
        HEPfitProposalConfig config;
        for (auto line : lines) config.ParseLine(line);
        config.ValidateRun(staged, smooth, loadRun, read, optimizer);
    } catch (const std::runtime_error&) { failed = true; }
    check(failed);
}

int main()
{
    HEPfitProposalConfig config;
    check(!config.Enabled() && config.warmup == 1000);
    config.ValidateRun(false, false, false, false, false);
    check(!config.ParseLine("PrerunMaxIter 12"));
    check(config.ParseLine(" ProposalWarmup 7 # discard"));
    check(config.ParseLine("SaveProposal \"path with spaces/saved#proposal.txt\" # comment"));
    check(config.ParseLine("LoadProposal input.txt\r"));
    config.ValidateRun(false, false, false, false, false);
    check(config.Enabled() && config.warmup == 7 && config.load == "input.txt");
    check(config.save == "path with spaces/saved#proposal.txt");
    config.ParseLine("ProposalWarmup 0");
    config.ValidateRun(false, false, false, false, false);
    check(config.warmup == 0);
    HEPfitProposalConfig save;
    save.ParseLine("SaveProposal tuned.txt");
    save.ValidateRun(true, true, false, false, false);
    reject({"SaveProposal"});
    reject({"LoadProposal"});
    reject({"LoadProposal # missing"});
    reject({"SaveProposal \"\""});
    reject({"LoadProposal \"unclosed"});
    reject({"SaveProposal a b"});
    reject({"ProposalWarmup"});
    reject({"ProposalWarmup 7"});
    reject({"LoadProposal in", "ProposalWarmup -1"});
    reject({"LoadProposal in", "ProposalWarmup 1.5"});
    reject({"LoadProposal in", "ProposalWarmup 2x"});
    reject({"LoadProposal in", "ProposalWarmup 2147483648"});
    reject({"LoadProposal in", "ProposalWarmup 9999999999999999999999"});
    reject({"LoadProposal in"}, true);
    reject({"LoadProposal in"}, false, true);
    reject({"LoadProposal in"}, false, false, true);
    reject({"LoadProposal in"}, false, false, false, true);
    reject({"SaveProposal out"}, false, false, true);
    reject({"SaveProposal out"}, false, false, false, true);
    reject({"LoadProposal in"}, false, false, false, false, true);
    reject({"SaveProposal out"}, false, false, false, false, true);
    std::cout << "Proposal cache paths, warmup counts and incompatible modes validated\n";
}
