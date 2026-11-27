# CMFV fit

Constrained Minimal Flavour Violation: new physics enters only as real shifts of the gauge-independent master
functions, `Ftt` (dS), `CMFV_dX`, `CMFV_dY`, `CMFV_dZ` and `CMFV_dC7` (model `CMFV`, see
`FlavourWilsonCoefficient/src/CMFV.h`). The likelihood is the UTfit Winter27/Current one plus the rare decays of
`Observables_CMFV.conf`:

| Observable | Measurement | Likelihood |
|---|---|---|
| BR(B -> Xs gamma), E > 1.6 GeV | 3.49 +- 0.19 (1e-4), HFLAV | Gaussian, times CMFV_r_bsg |
| BR(B -> Xs l l), [1, 6] GeV^2 | BaBar 1.60 +- 0.46, Belle 1.49 +- 0.63 (1e-6) | one Gaussian (their product), times CMFV_r_bsll |
| BR(Bs -> mu mu), BR(B0 -> mu mu) | 3.40 -0.23 +0.25 (1e-9), 0.85 -0.32 +0.34 (1e-10) | split Gaussians, times CMFV_r_Bmumu |
| BR(K+ -> pi+ nu nu) | 9.6 -1.8 +1.9 (1e-11), NA62 | split Gaussian |
| BR(KL -> pi0 nu nu) | < 2.2e-9 at 90% CL, KOTO | Gaussian at 0, sigma = limit / 1.645 |
| BR(B+ -> K+ nu nu) | 2.3 +- 0.7 (1e-5), Belle II | Gaussian (short distance) |
| BR(B0 -> K*0 nu nu) | < 1.8e-5 at 90% CL, Belle | Gaussian at 0, sigma = limit / 1.645 |

`CMFV_r_bsg`, `CMFV_r_bsll`, `CMFV_r_Bmumu` are multiplicative theory-uncertainty factors with Gaussian priors of
4.2%, 7.5% and 1.5% (`CMFV.conf`), applied through theory names of the form `<parameter>*<observable>`. The B_q ->
mu mu values are an average of the LHCb Run 1-3, CMS and ATLAS results, which replaces the UTfit Bs -> mu mu input.
eps'/eps is not used: the CMFV shifts do not enter the Delta S = 1 matching of K -> pi pi yet.

## Run directories

```sh
Configurations/Temporary/CMFV/setup_cmfv.sh <dir> [--zdom] [--smlike] [--predictions]
cd <dir> && mpirun -n <N> analysis StandardModel.conf MonteCarlo_PT.conf
```

| Option | Fit |
|---|---|
| (none) | all five shifts free on dS [-2, 4], dX [-6, 3], dY [-3, 1.5], dZ [-4, 4], dC7 [-0.6, 1.5] |
| `--smlike` | only the SM-like solution (no flipped sign of X, Y, C7): dX > -2, dY > -0.95, dC7 < 0.5 |
| `--zdom` | Z0-penguin dominance, dX = dY = dZ = dC (flag CMFV_ZPenguinDominance), dC = CMFV_dX in [-3, 1.5] |
| `--zdom --smlike` | Z0-penguin dominance, SM-like solution: dC in [-0.95, 1.5], dC7 < 0.5 |
| `--predictions` | rare decays as noMCMC predictions only (for a likelihood computed outside HEPfit) |

The free fit has separated sign-flipped solutions of X, Y and C7 (C7 together with C9), hence parallel tempering:
`MonteCarlo_PT.conf` is a starting point (4 ladders x 6 temperatures, pilot length). Check the mixing between
solutions with `bat_MPI/diagnostics/analyze_mcmc.py --regions`, with regions split at dX = -2, dY = -0.95 and
dC7 = 0.5. Along the Delta M_s constraint (dS against B_Bs, f_Bs and the CKM factors) an ordinary run of the
Delta S sector needed ~1e6 iterations to converge.

Cost: one likelihood evaluation takes about 0.6 s (lxplus, one core), almost all of it BR(B -> Xs l l); the other
rare decays take 1-3 ms each.

ROOT prints "Illegal data type" for UTfit parameters such as BK(1/2)1 when writing the chains: their names are not
valid branch names, so they are missing from the parameter tree; the fit is not affected.

## Validation

On random points, the rare-decay log-likelihood computed here equals the one of an independent implementation
(MFVfit, used for the nested-sampling fits) up to the Gaussian normalisation constant, with and without the
theory factors (spread 5e-9, from the 10-digit combined B -> Xs l l measurement). With CMFV_ZPenguinDominance the
predictions are identical to those with dY = dZ = dX set explicitly. The CMFV responses BR(shift)/BR(SM) agree with the validated earlier version of the
model to 0.07-0.3%.
