# Parallel Tempering and MCMC diagnostics in HEPfit

HEPfit uses the Parallel Tempering (PT) implementation in the external `bat_MPI`
repository. Configure HEPfit against the upgraded library; its bundled legacy BAT
MPI patch does not contain PT. Ordinary MCMC remains supported with legacy BAT.
The statistical analyzer remains in `bat_MPI/diagnostics/analyze_mcmc.py`;
HEPfit-specific region definitions belong in your HEPfit analysis directory.

## Build

Build and install `bat_MPI` to a dedicated prefix, then build HEPfit against it:

```sh
cmake -S /path/to/bat_MPI -B /path/to/bat_MPI/build-pt \
  -DCMAKE_BUILD_TYPE=Release
cmake --build /path/to/bat_MPI/build-pt -j 4
cmake --install /path/to/bat_MPI/build-pt --prefix /path/to/bat-pt
cmake -S /path/to/HEPfit -B /path/to/HEPfit/build-pt \
  -DMPIBAT=ON -DBAT_INSTALL=OFF -DBAT_INSTALL_DIR=/path/to/bat-pt \
  -DRGESOLVER_CONFIG_DIR=/path/to/RGESolver/bin
cmake --build /path/to/HEPfit/build-pt --target analysis -j 4
```

HEPfit prefers the exported BAT CMake package when found in BAT_INSTALL_DIR.
It links BAT::BAT so ROOT and MPI dependencies follow the selected installation.
Legacy installations without a CMake package remain supported for ordinary runs.

Use your usual ROOT, GSL, Boost and compiler settings. The configuration must
report `BAT supports Parallel Tempering and transport diagnostics`. An explicit
`BAT_INSTALL_DIR` takes precedence over environment search paths. Keep
`BAT_INSTALL=OFF` to prevent the legacy BAT download from replacing the upgraded
library. The optional `BATmodels` library is not needed by the analysis executable.

## Configure and run

Start with `Analysis/MonteCarlo_PT.conf`. The new options are:

| Option | Meaning |
| --- | --- |
| `ParallelTempering true` | Enable likelihood tempering; default is false. |
| `PTTemperatures 1:2:4:8:16:32` | Strictly increasing finite temperatures starting at 1; at least two. Commas or whitespace also work. |
| `PTSwapInterval 10` | Attempt exchanges every 10 sampler sweeps; positive integer, default 10. |
| `PTDiagnostics true` | Export production swap and replica round-trip CSV files; default true. |

With PT, `NChains` is the number of independent ladders and must be at least two.
For example, four ladders with six temperatures have 24 replicas and four cold
posterior chains. Configuration option order does not matter. MPI process count
can differ from replica count, including one process. More workers do not add
independent ladders. `NChains 0` retains the legacy meaning of MPI process count;
use an explicit ladder count for reproducible setup.

```sh
mpirun -n 8 /path/to/HEPfit/build-pt/Analysis/analysis \
  Model.conf /path/to/HEPfit/Analysis/MonteCarlo_PT.conf --rootfile fit
```

Pass the configuration filename through a cluster job script. The colon
separator also avoids the comma interpretation of `qsub -v` if you generate
configuration values from environment variables.

The sampler targets `prior(theta) * likelihood(theta)^(1/T)` at temperature T.
HEPfit's model configuration still determines the prior support. Remove any
analysis-specific prior cuts that would exclude modes you want to explore;
PT cannot cross a region excluded by the prior. Each ladder has its own cold
chain. BAT's adaptation takes place during prerun and the production proposal
is fixed. Swaps and round-trip counters reported below exclude prerun.

PT forces cold-chain storage even with `WriteChain false`, because diagnostics
need ordered samples. All marginalized parameter, observable, correlation and
log-likelihood histograms use only cold chains. HEPfit's observable and physical
parameter trees also contain only cold chains, with `Chain`, `Iteration`, `Phase`
and `Replica` branches. `Chain` identifies the fixed cold slot; `Replica` tracks
the moving walker. The physical-parameter tree requires the existing
`WriteParametersChains true` option. Unrotated BAT sampling parameters are in
`<model>_mcmc`; HEPfit's physical parameters are in `<model>_Parameters`.
Observable values are in `<model>_Observables`. Phase 1 means production.
These additional metadata branches are also written for ordinary runs.

Legacy `--loadRun` and `ReadPreRunData` are rejected for PT: they do not restore
a tempered ensemble. `RunMinuitOnly` and `ComputeHessianOnly` cannot enable PT.

## Independent warmup controls (ordinary MCMC and PT)

Both patches are opt-in and use the implementation in `bat_MPI`. Reinstall BAT
and reconfigure/rebuild HEPfit against the same installation. CMake reports
`BAT staged warmup: 1; smooth adaptation: 1` when both APIs are present. Older BAT
remains supported with the default controls; requesting an unavailable feature
produces an error before sampling.

| Option | Meaning |
| --- | --- |
| `PreRunPolicy legacy` | Default: existing acceptance/R-value prerun stopping policy. |
| `PreRunPolicy staged` | Adapt for exactly `PrerunMaxIter` sweeps, freeze the proposal, discard a separate settling phase, then start production. |
| `SettlingSweeps 1000` | Number of discarded fixed-proposal sweeps, default 1000; specifying it requires `PreRunPolicy staged`. Zero is allowed. |
| `SmoothAdaptation false` | Default: existing proposal scale updates. |
| `SmoothAdaptation true` | Second patch: bounded, diminishing log-scale updates during adaptation, independently of the stopping policy. |

The staged policy ignores the legacy minimum-prerun and early-stop gates.
Acceptance estimates and legacy R values remain available as tuning information.
Covariance and scale adaptation stop before settling; settling samples are excluded
from chains, histograms and posterior diagnostics. `Iterations` still counts only
production sweeps. Completion of warmup does not certify posterior convergence.

Smooth adaptation targets local acceptance 0.25, with gain 2 and decay 0.6.
At update window k, the learning rate is `2/k^0.6`; the log standard deviation
step is that rate times `(acceptance - 0.25)`, clipped to [-0.5, 0.5]. Joint
proposals apply twice this step to their variance multiplier. Configured scale
limits still apply. Coordinate and joint phases have separate update counters,
and covariance-shape adaptation is unchanged. These engineering defaults do not
guarantee better mixing for a particular fit.

For a controlled comparison, make two copies of your existing configuration.
Use identical temperatures, chains, proposal type, adaptation/production budgets
and `NIterationsUpdateMax` in both. For example, these budget values are
illustrative and require assessment for your model:

```text
# Test A: first patch, existing scale tuning
PreRunPolicy staged
PrerunMaxIter 50000
SettlingSweeps 20000
Iterations 100000
SmoothAdaptation false
```

Test B uses the same settings with `SmoothAdaptation true`. Keep each keyword
once in each file. If using multivariate proposals, both configurations need
`NIterationsPreRunFactorized < PrerunMaxIter`, so adaptation actually reaches
the joint kernel. `MultivariateProposal false` selects coordinate proposals.
Smooth tuning can also be tested with `PreRunPolicy legacy`, without
`SettlingSweeps`, to retain the original stopping policy.

Run the two jobs in separate working directories to keep HEPfit's other outputs
separate, and use distinct ROOT prefixes or job tags:

```sh
# In the Test A working directory
mpirun -n 8 /path/to/analysis Model.conf MonteCarlo_A.conf --rootfile fit_A
# In the Test B working directory
mpirun -n 8 /path/to/analysis Model.conf MonteCarlo_B.conf --rootfile fit_B
```

With the new BAT API, each fresh MCMC run writes `<prefix>_warmup_status.csv`,
including any job tag. This file is updated during warmup and records the policy,
completion and stop reason, adaptation/settling sweep counts, acceptance outliers,
zero-acceptance kernels, scale method and its settings, and separate update counts.
`posterior_convergence_assessed` remains false: use production R-hat, ESS, MCSE,
and regional/transport diagnostics to assess the fit. Compare scientific precision
per time or likelihood cost, rather than acceptance alone.

PT streaming diagnostics distinguish `prerun` (adaptation), `warmup` (discarded
fixed settling) and `production`. All use the existing prefix; the current BAT
analyzer reads the status file into `summary.json` and prefers production when
available, otherwise settling, otherwise adaptation. `PTDiagnostics false` only
suppresses PT transport output; the warmup status is still written. Ordinary runs
also receive the status file, without PT tables.

The new modes reject `--loadRun`, `ReadPreRunData`, `RunMinuitOnly` and
`ComputeHessianOnly`, where adaptation would be skipped or legacy initialization
would conflict with these tests. `ReadPreRunData` is the older HEPfit position/scale
file. Use the full proposal-cache interface below to reuse a tuned kernel.

## Saving and reloading full proposals through MonteCarlo.conf

Reconfigure and rebuild HEPfit against an updated external `bat_MPI` installation.
CMake reports `BAT proposal cache: 1` when the API is present. These options work
in ordinary MCMC and PT, with coordinate or multivariate proposals:

| Option | Meaning |
| --- | --- |
| `SaveProposal tuned_proposal.txt` | Save the tuned proposal and starting positions after warmup, before production. Optional; no cache is written unless requested. |
| `LoadProposal tuned_proposal.txt` | Restore the proposal/positions, skip all adaptation, then perform discarded fixed-proposal warmup before production. |
| `ProposalWarmup 1000` | Discard this many full sweeps on reload; default 1000, zero allowed, maximum INT_MAX. Specifying it requires `LoadProposal`. |

Paths are relative to the process's working directory, rather than the
configuration file's directory. Absolute paths are also accepted. Quote paths
containing spaces, for example `SaveProposal "cache/tuned proposal.txt"`.
The parent directory must already exist. Inline comments are allowed after a
value. Option order does not matter.

For a fresh fit with staged warmup and smooth adaptation, add:

```text
PreRunPolicy staged
PrerunMaxIter 50000
SettlingSweeps 20000
SmoothAdaptation true
SaveProposal fit_proposal.txt
Iterations 100000
```

For a subsequent longer production fit, retain the same model/prior, chains,
temperatures and `MultivariateProposal` settings. Replace the fresh adaptation
options above with:

```text
LoadProposal fit_proposal.txt
ProposalWarmup 20000
Iterations 1000000
```

The numerical lengths are examples, not convergence requirements. On reload,
omit `PreRunPolicy staged`, `SettlingSweeps` and `SmoothAdaptation true`: the
loaded proposal is already fixed, so these settings are rejected. Explicit
`PreRunPolicy legacy` and `SmoothAdaptation false` are allowed. Prerun iteration
settings do not trigger fresh adaptation after a successful load. Use a new ROOT
prefix or job directory for the new production output.

The cache retains coordinate scales or the full multivariate covariance and
actual Cholesky transform, proposal distribution settings, replica temperatures,
parameter metadata and positions. BAT validates sampling parameter names/order,
bounds/fixed values, proposal type, replica counts and temperature ordering;
incompatible or malformed files produce errors. With PT, `NChains` must describe
the same ladder layout as the cache. MPI worker count may change. Preserve the
physical model, likelihood, priors and any correlated-parameter definitions:
the file is not a complete HEPfit model configuration or a model fingerprint.

Saved log probabilities, histograms, convergence flags and random-generator
state are not restored. The current likelihood/prior are evaluated at the
loaded positions. This is proposal reuse with a new production run, not an exact
continuation from the previous production endpoint. The cache is saved after
adaptation/settling; it is not refreshed with the final production position.

`SaveProposal` can also accompany `LoadProposal` to export the unchanged kernel
with positions reached after the new discarded warmup. BAT replaces the output
cache atomically after a successful write; use distinct paths for concurrent
writers. The warmup status records `proposal_loaded=true`, zero adaptation sweeps,
the actual discarded sweep count and `stop_reason=loaded_proposal`. Scale-method
settings on reload describe the current run; the cache does not record which
adaptation method originally trained the proposal. PT exports have `warmup` and
`production` phases, and no adaptive `prerun` phase. All saved chains and
histograms still exclude discarded warmup.

Both cache options reject `--loadRun`, `ReadPreRunData`, `RunMinuitOnly` and
`ComputeHessianOnly`, avoiding mixed restoration paths or silently unused cache
requests. Legacy `WritePreRunData` output is a separate format and cannot be
passed to `LoadProposal`. Older BAT builds retain their existing behavior when
cache options are absent and reject cache requests explicitly when unsupported.

## Analyze the output

For `--rootfile fit` without a job tag, HEPfit produces `fit.root` plus:

- `fit_pt_swaps.csv`: aggregate acceptance per adjacent temperature pair.
- `fit_pt_swaps_by_ladder.csv`: attempts and acceptance for each ladder/pair.
- `fit_pt_roundtrips.csv`: cold–hot–cold round trips and endpoint visits per walker.

The prefix includes any job tag; pass that same prefix to the analyzer.
`PTDiagnostics false` suppresses the CSV export, but cold-chain storage remains
active. Python is needed only for offline diagnostics:

```sh
python3 -m venv diagnostics-venv
. diagnostics-venv/bin/activate
python -m pip install -r /path/to/bat_MPI/diagnostics/requirements.txt
python /path/to/bat_MPI/diagnostics/analyze_mcmc.py fit.root \
  --output diagnostics_parameters --transport-prefix fit_ --fail-on-warning
python /path/to/bat_MPI/diagnostics/analyze_mcmc.py fit.root \
  --tree StandardModel_Observables --var Mw \
  --output diagnostics_observables --transport-prefix fit_
```

Replace the model/tree/variable names with those in your output. Select varying
columns with repeated `--var` options when using HEPfit custom trees, which can
contain fixed physical parameters or constant predictions without BAT's fixed
parameter metadata. For physical parameters use `--tree <model>_Parameters`.
The same analyzer works on ordinary runs with saved chains, omitting
`--transport-prefix`. Existing old observable trees without `Phase` cannot be
analyzed directly; use their BAT `<model>_mcmc` tree.

The analyzer provides rank-normalized split/folded R-hat, bulk/tail ESS, mean
MCSE, trace/rank plots, per-ladder swap statistics, and round trips. Explicit
region definitions via `--regions regions.json` add mode weights by chain,
block and cumulative time, indicator ESS and mode-weight MCSE. Copy and adapt
the schema from `bat_MPI/diagnostics/regions.example.json`. Use `--angular name`
for phase variables that do not follow the `_arg` naming convention.

Assess both cold-posterior convergence and transport. Good swap acceptance
alone does not establish mixing between modes. Check round trips, agreement of
mode weights across independent ladders/time blocks, R-hat, ESS and MCSE for
parameters and scientific predictions. Reported numerical thresholds are
screening criteria, not proof that all modes have been found. Constant or
unvisited mode indicators cannot provide a reliable zero MCSE. The analyzer's
`--fail-on-warning` exits 2 for diagnostic warnings/incomplete diagnostics,
1 for invalid input and 0 when requested checks pass.

## Regression tests

BAT owns the generic sampler and statistical estimator tests. HEPfit adds
configuration validation and an MPI integration test using a one-parameter
Standard Model EW fit. Enable them with:

```sh
cmake -S /path/to/HEPfit -B /path/to/HEPfit/build-pt \
  -DHEPFIT_BUILD_PT_TESTS=ON \
  -DHEPFIT_DIAGNOSTICS_PYTHON=/path/to/diagnostics-venv/bin/python \
  -DHEPFIT_BAT_SOURCE_DIR=/path/to/bat_MPI
cmake --build /path/to/HEPfit/build-pt --target analysis HEPfitPTConfigTest HEPfitWarmupConfigTest HEPfitProposalConfigTest -j 4
ctest --test-dir /path/to/HEPfit/build-pt --output-on-failure
```

The integration test requires PT-enabled BAT and MPIBAT. It checks MPI and
single-process PT, ordinary MPI MCMC, cold-only histograms/trees, synchronized
replica identities and likelihoods, transport exports, offline parameter and
observable analysis, and termination on invalid PT configuration. Configure
`MPIEXEC_PREFLAGS` when your MPI installation needs launcher-specific flags.

With both new BAT APIs present, a second MPI integration test covers staged
warmup with legacy/smooth scales and coordinate/joint proposals, for ordinary and
PT runs, plus single-process PT, zero settling, smooth adaptation with legacy
stopping, and status output with PT exports disabled. It checks exact sweep budgets,
production-only ROOT entries and histogram counts, synchronized likelihoods and
replica metadata, independent phase exports, status ingestion by the analyzer,
and termination of MPI workers on invalid warmup options. Configuration tests
also run independently of the BAT feature set. Generic kernel-freezing, update
mathematics and posterior-moment regression tests remain in `bat_MPI`.

When the cache API is present, a proposal integration test saves/reloads real
HEPfit EW fits for ordinary/PT sampling and both proposal types. It compares
every saved kernel field on reload, checks zero/nonzero discarded warmup,
changing worker counts, production-only trees/histograms, phase/status exports,
and clean MPI termination for invalid configuration, missing/corrupt caches,
incompatible bounds, temperatures, replica counts and proposal types. The
standalone proposal configuration test also runs with legacy BAT installations.


## Periodic all-temperature diagnostics

With a newly rebuilt bat_MPI exposing `SetPTDiagnosticsOutput`, `PTDiagnostics true`
also enables periodic sampler output: exact counter windows every 1000 sweeps and
state traces every 100 sweeps. Files use the existing output prefix. Counters cover
local proposal acceptance per temperature/kernel, block swap rates, and decimated
physical replica trajectories, separately for prerun and production. HEPfit currently
provides no model-specific region classifier, so regional tables are empty; a model
can implement BAT's generic partition/classifier hooks. GoldenModes definitions
remain in its own repository. `PTDiagnostics false` disables these files too.

Rebuild HEPfit after upgrading BAT; matching headers and shared libraries are required.
CMake detects the new API separately, so older PT-enabled BAT builds retain their
existing end-of-run transport diagnostics. For new output, the BAT analyzer supports:

```sh
python /path/to/bat_MPI/diagnostics/analyze_mcmc.py --pt-only \
  --transport-prefix /path/to/output_ --output /path/to/transport-report
```

This also works before a job finishes: it reads complete committed counter windows,
marks unfinished phases, and ignores uncommitted trailing output. It does not check
cold posterior R-hat/ESS. Run the ordinary ROOT analysis for those checks. Details of
all files and diagnostic interpretation are in `bat_MPI/diagnostics/README.txt`.
