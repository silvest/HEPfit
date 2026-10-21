#!/usr/bin/env python3
"""Exercise HEPfit's MPI worker protocol, cold-only output and BAT diagnostics."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import uproot


def run(command, cwd, expect_success=True):
    result = subprocess.run(command, cwd=cwd, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=75)
    if expect_success != (result.returncode == 0):
        raise AssertionError(f"Command failed expectation: {command}\n{result.stdout[-14000:]}")
    return result


def model_config():
    # One varying parameter, a likelihood and a prediction through real HEPfit EW code.
    values = dict(lambda_=0.225787, A=0.826, rhob=0.148, etab=0.348,
                  Mz=91.1875, AlsMz=0.1185, GF=1.1663787e-5,
                  ale=7.2973525698e-3, dAle5Mz=0.0275, mHl=125.09,
                  delMw=0, delSin2th_l=0, delSin2th_q=0, delSin2th_b=0,
                  delGammaZ=0, delsigma0H=0, delR0l=0, delR0c=0, delR0b=0,
                  delGammaWlv=0, delGammaWqq=0,
                  mneutrino_1=0, mneutrino_2=0, mneutrino_3=0,
                  melectron=5.109989e-4, mmu=0.10565837, mtau=1.77682,
                  muw=80.385, mup=0.0023, mdown=0.0048, mcharm=1.3,
                  mstrange=0.0938, mtop=163.0, mbottom=4.18,
                  muc=1.3, mub=4.8, mut=164.1)
    lines = ['StandardModel']
    for name, value in values.items():
        name = name.rstrip('_')
        lines.append(f'ModelParameter {name} {value} {2.0 if name == "mtop" else 0} 0')
    lines += ['Observable Mw Mw Mw 79 82 MCMC weight 80.36 0.04 0',
              'Observable MwPrediction Mw MwPrediction 79 82 noMCMC noweight']
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--expect-streaming', action='store_true')
    parser.add_argument('--analysis', required=True)
    parser.add_argument('--bat-source', required=True)
    parser.add_argument('--mpiexec', default='mpirun')
    parser.add_argument('--mpi-flags', default='')
    args = parser.parse_args()
    analysis = str(Path(args.analysis).resolve())
    diagnostic = Path(args.bat_source).resolve() / 'diagnostics/analyze_mcmc.py'
    flags = args.mpi_flags.split(';') if args.mpi_flags else []
    flags = [flag for flag in flags if flag]
    common = '''NChains 2
PrerunMaxIter 100
NIterationsUpdateMax 50
NIterationsPreRunFactorized 0
Iterations 400
Seed 9876
WriteChain true
WriteLogLikelihoodChain true
WriteParametersChains true
FindModeWithMinuit false
CalculateNormalization false
PrintAllMarginalized false
PrintCorrelationMatrix false
PrintKnowledgeUpdatePlots false
PrintParameterPlot false
PrintTrianglePlot false
MultivariateProposal true
HistogramBufferSize 0
'''
    with tempfile.TemporaryDirectory(prefix='hepfit-pt-') as temporary:
        base = Path(temporary)
        for label, processes, pt in [('pt_partial', 3, True), ('pt_serial', 1, True), ('ordinary', 3, False)]:
            work = base / label
            work.mkdir()
            (work / 'model.conf').write_text(model_config())
            # PT before NChains tests that ladder setup is independent of option order.
            options = ('ParallelTempering true\nPTTemperatures 1:2:4\nPTSwapInterval 5\n' if pt else '')
            # No trailing newline exercises the runtime parser, not only the helper.
            export_transport = pt and processes > 1
            settings = common
            if pt and processes == 1:
                settings = settings.replace('WriteChain true', 'WriteChain false')
                options += 'PTDiagnostics false\n'
            (work / 'mc.conf').write_text(options + settings.rstrip('\n'))
            command = [args.mpiexec, '-n', str(processes), *flags, analysis, 'model.conf', 'mc.conf']
            run(command, work)
            with uproot.open(work / 'MCout.root') as root:
                trees = {key: root[f'StandardModel_{key}'].arrays(library='np')
                         for key in ('mcmc', 'Observables', 'Parameters')}
                for key, arrays in trees.items():
                    assert len(arrays['Chain']) == 800, (label, key, len(arrays['Chain']))
                    assert set(arrays['Chain']) == {0, 1}
                    assert set(arrays['Phase']) == {1}
                    np.testing.assert_array_equal(arrays['Chain'], trees['mcmc']['Chain'])
                    np.testing.assert_array_equal(arrays['Iteration'], trees['mcmc']['Iteration'])
                    np.testing.assert_array_equal(arrays['Replica'], trees['mcmc']['Replica'])
                np.testing.assert_allclose(trees['Observables']['LogLikelihood'], trees['mcmc']['LogLikelihood'])
                np.testing.assert_allclose(trees['Parameters']['mtop'], trees['mcmc']['mtop'])
                np.testing.assert_allclose(trees['Observables']['Mw'], trees['Observables']['MwPrediction'])
                histogram_count = 0
                for key in root.keys():
                    obj = root[key]
                    if obj.classname.startswith('TH1') and any(n in key for n in ('Mw', 'LogLikelihood')):
                        histogram_count += 1
                        assert obj.member('fEntries') == 800, (label, key, obj.member('fEntries'))
                assert histogram_count >= 3
            if export_transport:
                with (work / 'MCout_pt_swaps_by_ladder.csv').open() as handle:
                    swaps = list(csv.DictReader(handle))
                assert len(swaps) == 4
                assert all(int(row['attempts']) == 40 for row in swaps)
                assert sum(int(row['accepted']) for row in swaps) > 0
                assert {int(row['ladder']) for row in swaps} == {0, 1}
                assert (work / 'MCout_pt_roundtrips.csv').is_file()
                assert (work / 'MCout_pt_swaps.csv').is_file()
                if args.expect_streaming:
                    for phase, sweeps in [('prerun',100),('production',400)]:
                        with (work/f'MCout_pt_progress_{phase}.csv').open() as handle:
                            progress=list(csv.DictReader(handle))
                        assert len(progress)==1 and progress[0]['completed']=='1'
                        assert int(progress[0]['last_sweep'])==sweeps
                    with (work/'MCout_pt_local_moves.csv').open() as handle:
                        counters=list(csv.DictReader(handle))
                    assert {float(row['temperature']) for row in counters}=={1.,2.,4.}
                    assert {int(row['slot']) for row in counters}==set(range(6))
                    assert {row['phase'] for row in counters}=={'prerun','production'}
                    with (work/'MCout_pt_partitions.csv').open() as handle:
                        assert not list(csv.DictReader(handle))  # Generic HEPfit has no region cuts.
                    run([sys.executable,str(diagnostic),'--pt-only','--transport-prefix','MCout_',
                         '--output','transport_report','--no-plots'],work)
            else:
                assert not list(work.glob('*pt_*.csv'))
            # Statistical convergence is not promised by this short integration fit.
            # Successful analysis/output is checked; warnings are expected and retained.
            for tree, variable in [('mcmc', 'mtop'), ('Observables', 'MwPrediction'), ('Parameters', 'mtop')]:
                run([sys.executable, str(diagnostic), 'MCout.root', '--tree', f'StandardModel_{tree}',
                     '--var', variable, '--output', f'diagnostics_{tree}', '--no-plots',
                     *(['--transport-prefix', 'MCout_'] if export_transport else [])], work)
            print(f'{label}: cold trees, histograms, replica identities and diagnostics passed')
        # Rejected configuration must terminate all MPI workers without hanging.
        for options, message in [('ParallelTempering true\nPTTemperatures 1 1', 'strictly increasing'),
                                 ('ParallelTempering true\nPTTemperatures 1 2\nNChains 1', 'independent ladders')]:
            (work / 'mc.conf').write_text(common + options)
            result = run([args.mpiexec, '-n', '3', *flags, analysis, 'model.conf', 'mc.conf'], work, False)
            assert message in result.stdout, result.stdout
        print('Invalid PT configurations terminate MPI workers correctly')


if __name__ == '__main__':
    main()
