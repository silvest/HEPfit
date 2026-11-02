#!/usr/bin/env python3
"""Check HEPfit's two independent warmup controls through real MPI fits."""
import argparse
import csv
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
import uproot

from test_pt_integration import model_config, run


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis', required=True)
    parser.add_argument('--bat-source', required=True)
    parser.add_argument('--mpiexec', default='mpirun')
    parser.add_argument('--mpi-flags', default='')
    parser.add_argument('--expect-streaming', action='store_true')
    args = parser.parse_args()
    flags = [flag for flag in args.mpi_flags.split(';') if flag]
    command = [args.mpiexec, '-n', '3', *flags, str(Path(args.analysis).resolve()),
               'model.conf', 'mc.conf', '--rootfile', 'fit', '--job_tag', '_test']
    diagnostic = Path(args.bat_source).resolve() / 'diagnostics/analyze_mcmc.py'
    common = '''NChains 2
PrerunMaxIter 12
NIterationsUpdateMax 5
NIterationsPreRunFactorized 3
Iterations 9
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
HistogramBufferSize 0
'''
    pt_options = 'ParallelTempering true\nPTTemperatures 1:2\nPTSwapInterval 2\n'
    with tempfile.TemporaryDirectory(prefix='hepfit-warmup-') as temporary:
        base = Path(temporary)
        cases = [(pt, joint, smooth, 7, 3)
                 for pt in (False, True) for joint in (False, True) for smooth in (False, True)]
        cases += [(True, True, True, 0, 1)]
        for number, (pt, joint, smooth, settling, processes) in enumerate(cases):
            label = f'{number}_pt{pt}_joint{joint}_smooth{smooth}'
            work = base / label
            work.mkdir()
            (work / 'model.conf').write_text(model_config())
            # Options precede the older settings: validation must use their final values.
            options = (f'SettlingSweeps {settling}\nSmoothAdaptation {str(smooth).lower()}\n'
                       f'PreRunPolicy staged\nMultivariateProposal {str(joint).lower()}\n')
            (work / 'mc.conf').write_text(options + (pt_options if pt else '') + common.rstrip('\n'))
            case_command = command.copy()
            case_command[2] = str(processes)
            run(case_command, work)
            prefix = work / 'fit_test_'
            status = rows(work / 'fit_test_warmup_status.csv')[0]
            assert status['policy'] == 'staged' and status['stage'] == 'complete', status
            assert status['adaptation_sweeps'] == '12' and int(status['settling_sweeps']) == settling
            assert status['adaptation_complete'] == status['settling_complete'] == '1'
            assert status['posterior_convergence_assessed'] == status['proposal_loaded'] == '0'
            assert status['stop_reason'] == 'adaptation_budget_exhausted'
            assert status['scale_adaptation'] == ('smooth' if smooth else 'legacy')
            assert float(status['scale_target']) == .25
            assert float(status['scale_gain']) == 2 and float(status['scale_decay']) == .6
            updates = [int(status[name]) for name in ('factorized_scale_updates', 'joint_scale_updates')]
            expected_updates = ([1, 2] if joint else [3, 0]) if smooth else [0, 0]
            assert updates == expected_updates, updates
            with uproot.open(work / 'fit_test.root') as root:
                trees = {key: root[f'StandardModel_{key}'].arrays(library='np')
                         for key in ('mcmc', 'Observables', 'Parameters')}
                for key, arrays in trees.items():
                    # Only nine production draws from each cold slot; no settling entries.
                    assert len(arrays['Chain']) == 18, (label, key, len(arrays['Chain']))
                    assert set(arrays['Chain']) == {0, 1} and set(arrays['Phase']) == {1}
                    for metadata in ('Chain', 'Iteration', 'Replica'):
                        np.testing.assert_array_equal(arrays[metadata], trees['mcmc'][metadata])
                np.testing.assert_allclose(trees['Parameters']['mtop'], trees['mcmc']['mtop'])
                np.testing.assert_allclose(trees['Observables']['LogLikelihood'], trees['mcmc']['LogLikelihood'])
                histograms = [root[key] for key in root.keys() if root[key].classname.startswith('TH1')
                              and any(name in key for name in ('Mw', 'LogLikelihood'))]
                assert len(histograms) >= 3
                assert all(hist.member('fEntries') == 18 for hist in histograms)
            if pt and args.expect_streaming:
                phases = [('prerun', 12), ('production', 9)]
                if settling:
                    phases.insert(1, ('warmup', settling))
                for phase, sweeps in phases:
                    progress = rows(work / f'fit_test_pt_progress_{phase}.csv')
                    assert progress[-1]['completed'] == '1' and int(progress[-1]['last_sweep']) == sweeps
                run([sys.executable, str(diagnostic), '--pt-only', '--transport-prefix', str(prefix),
                     '--output', 'transport_report', '--no-plots'], work)
                summary = json.loads((work / 'transport_report/summary.json').read_text())
                recorded = summary['streaming_pt']['runs'][0]['warmup_status']
                assert recorded['adaptation_sweeps'] == 12 and recorded['settling_sweeps'] == settling
                assert recorded['scale_adaptation'] == status['scale_adaptation']
            if not pt:
                assert not list(work.glob('*pt_*.csv'))
            print(f'{label}: exact budgets, separate controls, production-only output and status passed')

        # The scale flag works independently of staged stopping and PT CSV output.
        work = base / 'legacy_smooth_no_transport'
        work.mkdir()
        (work / 'model.conf').write_text(model_config())
        (work / 'mc.conf').write_text(common + pt_options +
            'PreRunPolicy legacy\nSmoothAdaptation true\nMultivariateProposal false\nPTDiagnostics false\n')
        run(command, work)
        status = rows(work / 'fit_test_warmup_status.csv')[0]
        assert status['policy'] == 'legacy' and status['scale_adaptation'] == 'smooth'
        assert status['stage'] == 'complete' and status['settling_sweeps'] == '0'
        assert int(status['factorized_scale_updates']) > 0 and status['joint_scale_updates'] == '0'
        assert not list(work.glob('*pt_*.csv'))
        print('Smooth adaptation with legacy stopping and disabled PT exports passed')

        # Runtime errors must abort the MPI worker loop promptly, before sampler collectives.
        for options, message in [
            ('PreRunPolicy staged\nSettlingSweeps -1', 'SettlingSweeps'),
            ('SettlingSweeps 7', 'requires PreRunPolicy staged'),
            ('PreRunPolicy other', 'legacy or staged'),
            ('SmoothAdaptation yes', 'true or false'),
            ('PreRunPolicy staged\nReadPreRunData unused', 'ReadPreRunData'),
            ('SmoothAdaptation true\nReadPreRunData unused', 'ReadPreRunData'),
            ('PreRunPolicy staged\nRunMinuitOnly true', 'require MCMC'),
            ('SmoothAdaptation true\nComputeHessianOnly true', 'require MCMC'),
            ('PreRunPolicy staged\nMultivariateProposal true\nNIterationsPreRunFactorized 12',
             'NIterationsPreRunFactorized < PrerunMaxIter'),
            ('SmoothAdaptation true\nMultivariateProposal true\nNIterationsPreRunFactorized 12',
             'NIterationsPreRunFactorized < PrerunMaxIter'),
            ('PreRunPolicy staged\nPrerunMaxIter 2147483647\nSettlingSweeps 1', 'INT_MAX'),
        ]:
            (work / 'mc.conf').write_text(common + options)
            result = run(command, work, False)
            assert message in result.stdout, result.stdout
        (work / 'mc.conf').write_text(common + 'PreRunPolicy staged\n')
        result = run([*command, '--loadRun', '--inFile', 'unused'], work, False)
        assert '--loadRun or ReadPreRunData' in result.stdout, result.stdout
        print('Invalid warmup configurations terminate MPI workers correctly')


if __name__ == '__main__':
    main()
