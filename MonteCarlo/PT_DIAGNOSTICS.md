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
cmake --build /path/to/HEPfit/build-pt --target analysis HEPfitPTConfigTest -j 4
ctest --test-dir /path/to/HEPfit/build-pt --output-on-failure
```

The integration test requires PT-enabled BAT and MPIBAT. It checks MPI and
single-process PT, ordinary MPI MCMC, cold-only histograms/trees, synchronized
replica identities and likelihoods, transport exports, offline parameter and
observable analysis, and termination on invalid PT configuration. Configure
`MPIEXEC_PREFLAGS` when your MPI installation needs launcher-specific flags.


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
