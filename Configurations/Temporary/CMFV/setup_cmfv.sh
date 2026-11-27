#!/bin/bash
# Assemble a run directory for the CMFV fit: the UTfit Winter27/Current configuration of this repository with the
# model CMFV, plus the files of this directory.
#   Configurations/Temporary/CMFV/setup_cmfv.sh <output dir> [--zdom] [--smlike] [--predictions]
#   (default) all five shifts dS (Ftt), dX, dY, dZ, dC7 free, rare-decay measurements in the likelihood
#   --zdom        Z0-penguin dominance, dX = dY = dZ = dC (flag CMFV_ZPenguinDominance; dC is CMFV_dX, on the
#                 overlap [-3, 1.5] of the dX, dY, dZ ranges; CMFV_dY and CMFV_dZ are fixed and ignored)
#   --smlike      only the SM-like solution, no flipped sign: dX > -2, dY > -0.95, dC7 < 0.5; with --zdom, dC > -0.95
#   --predictions rare decays as noMCMC predictions only (likelihood computed outside HEPfit)
# then: cd <output dir> && mpirun -n <N> analysis StandardModel.conf MonteCarlo_PT.conf
set -e
HERE=$(cd $(dirname $0) && pwd); HF=$(cd $HERE/../../.. && pwd)
OUT=$1; shift || true
[ -n "$OUT" ] || { sed -n 2,10p $0; exit 1; }
ZDOM=0; SMLIKE=0; PRED=0
for a in "$@"; do
  case $a in
    --zdom) ZDOM=1 ;;
    --smlike) SMLIKE=1 ;;
    --predictions) PRED=1 ;;
    *) echo "unknown option $a"; exit 1 ;;
  esac
done
mkdir -p $OUT && cp -r $HF/Configurations/Temporary/UTFit/Winter27/Current/. $OUT/ && cd $OUT
rm -f generateConf.sh run_*.sh submit_job*.sh NPDF2.conf NPWC*.conf UTfit_NP.conf
cp $HERE/*.conf .
# B -> K(*) form factors: the dedicated B -> K nu nu configuration when available (cubic DM expansion)
FF=$HF/Configurations/Temporary/BKnunu/SM/FormFactors_DM.conf
[ -f $FF ] && cp $FF FormFactors_DM.conf
sed -i '1s/^StandardModel\s*$/CMFV/' StandardModel.conf
grep -qE '^\s*ModelFlag\s+FlagMtPole\s+true' *.conf || sed -i '1a ModelFlag       FlagMtPole              true' StandardModel.conf
# cubic DM form-factor parameters: only if MPll needs them and the configuration does not define them
if grep -q b_3_fplus $HF/Flavour/src/MPll.cpp && ! grep -qhE '^\s*ModelParameter\s+b_3_fplus\s' $(ls *.conf | grep -v cubic_placeholder); then
  echo "NOTE: using FormFactors_DM_cubic_placeholder.conf (b_3 = 0)"
else rm -f FormFactors_DM_cubic_placeholder.conf; fi
# extra parameters already defined by UTfit
for p in $(grep -oE '^ModelParameter\s+\S+' Flavour_extra.conf | awk '{print $2}'); do
  grep -q -E "^ModelParameter\s+$p\s" SM_parameters.conf Flavour.conf UTfit.conf UTfit_SM.conf && sed -i -E "/^ModelParameter\s+$p\s/d" Flavour_extra.conf
done
if [ $PRED = 1 ]; then
  mv Observables_CMFV_predictions.conf Observables_CMFV.conf
  sed -i -E '/^ModelParameter\s+CMFV_r_/d' CMFV.conf
else rm -f Observables_CMFV_predictions.conf; fi
# flat prior on [lo, hi]: "ModelParameter name ave 0. errf"
setflat() { awk -v lo=$2 -v hi=$3 'BEGIN{printf "%g %g\n", (lo + hi) / 2, (hi - lo) / 2}' | { read a f; sed -i -E "s/^(ModelParameter\s+$1\s+).*/\1$a     0.   $f/" CMFV.conf; }; }
if [ $ZDOM = 1 ]; then
  setflat CMFV_dX -3 1.5; setflat CMFV_dY 0 0; setflat CMFV_dZ 0 0
  echo "ModelFlag       CMFV_ZPenguinDominance  true" >> CMFV.conf
  if [ $SMLIKE = 1 ]; then setflat CMFV_dX -0.95 1.5; fi
elif [ $SMLIKE = 1 ]; then
  setflat CMFV_dX -2 3; setflat CMFV_dY -0.95 1.5
fi
if [ $SMLIKE = 1 ]; then setflat CMFV_dC7 -0.6 0.5; fi
for f in Flags_bsll.conf CMFV.conf Flavour_extra.conf FormFactors_DM.conf FormFactors_DM_cubic_placeholder.conf h_lambda.conf Observables_CMFV.conf; do
  [ -f $f ] && { grep -q "IncludeFile $f" StandardModel.conf || echo "IncludeFile $f" >> StandardModel.conf; }
done
# eps'/eps: the CMFV shifts do not enter the Delta S = 1 matching of K -> pi pi yet; B_s -> mu mu: our average of the
# LHCb Run 1-3, CMS and ATLAS results enters through Observables_CMFV.conf instead of the UTfit one
sed -i -E 's/^(\s*Observable\s+(EpsilonP_O_Epsilon|BRbar_Bsmumu)\s.*)\bMCMC(\s+)weight/\1noMCMC\3weight/' UTfit_SM.conf
grep -E '^\s*Observable\s+(EpsilonP_O_Epsilon|BRbar_Bsmumu)\s' UTfit_SM.conf | grep -q -w MCMC && { echo "failed to disable eps'/eps, Bs->mumu"; exit 1; }
echo "CMFV configuration written to $OUT"
