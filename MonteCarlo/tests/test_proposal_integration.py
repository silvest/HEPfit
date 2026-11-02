#!/usr/bin/env python3
"""Round-trip HEPfit proposal caches through MonteCarlo.conf and real MPI fits."""
import argparse
import csv
import json
from pathlib import Path
import shlex
import sys
import tempfile

import numpy as np
import uproot

from test_pt_integration import model_config, run


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def kernel(path):
    """Independent cache reader: compare every proposal field, excluding positions."""
    tokens = iter(shlex.split(path.read_text()))
    assert next(tokens) == 'BAT_MPI_PROPOSAL' and next(tokens) == '1'
    n, dimension, chains, joint = [int(next(tokens)) for _ in range(4)]
    dof, updates = next(tokens), next(tokens)
    parameters = [tuple(next(tokens) for _ in range(5)) for _ in range(n)]
    replicas = []
    for _ in range(chains):
        temperature, size = next(tokens), int(next(tokens))
        scales = [next(tokens) for _ in range(size)]
        matrices = [next(tokens) for _ in range(2 * dimension * dimension)] if joint else []
        positions = [next(tokens) for _ in range(n)]
        assert len(positions) == n
        replicas.append((temperature, scales, matrices))
    assert next(tokens) == 'END' and next(tokens, None) is None
    return (n, dimension, chains, joint, dof, updates, parameters, replicas)


def output_checks(work, prefix):
    with uproot.open(work / (prefix + '.root')) as root:
        trees = {key: root[f'StandardModel_{key}'].arrays(library='np')
                 for key in ('mcmc', 'Observables', 'Parameters')}
        for key, arrays in trees.items():
            assert len(arrays['Chain']) == 18, (prefix, key, len(arrays['Chain']))
            assert set(arrays['Phase']) == {1} and set(arrays['Chain']) == {0, 1}
            for metadata in ('Chain', 'Iteration', 'Replica'):
                np.testing.assert_array_equal(arrays[metadata], trees['mcmc'][metadata])
        np.testing.assert_allclose(trees['Parameters']['mtop'], trees['mcmc']['mtop'])
        np.testing.assert_allclose(trees['Observables']['LogLikelihood'], trees['mcmc']['LogLikelihood'])
        histograms = [root[k] for k in root.keys() if root[k].classname.startswith('TH1')
                      and any(n in k for n in ('Mw', 'LogLikelihood'))]
        assert len(histograms) >= 3 and all(h.member('fEntries') == 18 for h in histograms)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis', required=True)
    parser.add_argument('--bat-source', required=True)
    parser.add_argument('--mpiexec', default='mpirun')
    parser.add_argument('--mpi-flags', default='')
    parser.add_argument('--expect-streaming', action='store_true')
    args = parser.parse_args()
    flags = [f for f in args.mpi_flags.split(';') if f]
    command = [args.mpiexec, '-n', '3', *flags, str(Path(args.analysis).resolve()), 'model.conf', 'mc.conf']
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
    with tempfile.TemporaryDirectory(prefix='hepfit-proposal-') as temporary:
        base = Path(temporary)
        for pt in (False, True):
            for joint in (False, True):
                work = base / f'pt{pt}_joint{joint}'
                work.mkdir()
                (work / 'model.conf').write_text(model_config())
                settings = common + f'MultivariateProposal {str(joint).lower()}\n' + (pt_options if pt else '')
                original = work / 'tuned proposal.txt'
                # Cache requests precede the remaining options: parsing is order independent.
                (work / 'mc.conf').write_text('SaveProposal "tuned proposal.txt"\n' + settings +
                    'PreRunPolicy staged\nSettlingSweeps 7\nSmoothAdaptation true\n')
                run([*command, '--rootfile', 'fresh'], work)
                expected = kernel(original)
                assert expected[2] == (4 if pt else 2) and bool(expected[3]) == joint
                status = rows(work / 'fresh_warmup_status.csv')[0]
                assert status['proposal_loaded'] == '0' and status['adaptation_sweeps'] == '12'
                output_checks(work, 'fresh')
                for sweeps in (0, 7):
                    prefix = f'loaded{sweeps}'
                    # Leave a factorized budget that would prevent fresh joint adaptation:
                    # a reload must bypass both adaptation and the coordinate-to-joint switch.
                    loaded = settings.replace('PrerunMaxIter 12', 'PrerunMaxIter 1')
                    (work / 'mc.conf').write_text(f'ProposalWarmup {sweeps}\nLoadProposal "tuned proposal.txt"\n'
                        f'SaveProposal reexport{sweeps}.txt\n' + loaded.rstrip('\n'))
                    case_command = command.copy()
                    if sweeps == 0: case_command[2] = '1'  # Worker count may differ from the saved replica layout.
                    run([*case_command, '--rootfile', prefix], work)
                    assert kernel(work / f'reexport{sweeps}.txt') == expected
                    status = rows(work / f'{prefix}_warmup_status.csv')[0]
                    assert status['proposal_loaded'] == '1' and status['adaptation_sweeps'] == '0'
                    assert int(status['settling_sweeps']) == sweeps and status['stop_reason'] == 'loaded_proposal'
                    assert status['stage'] == 'complete' and status['settling_complete'] == '1'
                    assert status['factorized_scale_updates'] == status['joint_scale_updates'] == '0'
                    output_checks(work, prefix)
                    if pt and args.expect_streaming:
                        assert not (work / f'{prefix}_pt_progress_prerun.csv').exists()
                        production = rows(work / f'{prefix}_pt_progress_production.csv')[-1]
                        assert production['completed'] == '1' and production['last_sweep'] == '9'
                        if sweeps:
                            warmup = rows(work / f'{prefix}_pt_progress_warmup.csv')[-1]
                            assert warmup['completed'] == '1' and int(warmup['last_sweep']) == sweeps
                            run([sys.executable, str(diagnostic), '--pt-only', '--transport-prefix', prefix + '_',
                                 '--output', prefix + '_report', '--no-plots'], work)
                            report = json.loads((work / f'{prefix}_report/summary.json').read_text())
                            assert report['streaming_pt']['runs'][0]['warmup_status']['proposal_loaded']
                    if not pt: assert not list(work.glob('*pt_*.csv'))
                print(f'PT={pt}, joint={joint}: saved/restored kernel, MPI/serial reload and discarded warmup passed')

        # Saving also works with the original prerun stopping and scale updates.
        (work / 'mc.conf').write_text(settings.replace(pt_options, '') + 'SaveProposal legacy-proposal.txt\n')
        run([*command, '--rootfile', 'legacy'], work)
        assert kernel(work / 'legacy-proposal.txt')[2] == 2
        status = rows(work / 'legacy_warmup_status.csv')[0]
        assert status['policy'] == status['scale_adaptation'] == 'legacy'
        assert status['proposal_loaded'] == '0' and int(status['adaptation_sweeps']) > 0
        output_checks(work, 'legacy')
        print('Ordinary legacy warmup exports a full proposal cache')

        # Current work directory contains a valid PT joint cache; all errors must terminate workers.
        for options, message in [
            ('ProposalWarmup 1', 'requires LoadProposal'),
            ('LoadProposal "tuned proposal.txt"\nProposalWarmup -1', 'ProposalWarmup'),
            ('LoadProposal "tuned proposal.txt"\nPreRunPolicy staged', 'skips adaptation'),
            ('LoadProposal "tuned proposal.txt"\nSmoothAdaptation true', 'skips adaptation'),
            ('LoadProposal "tuned proposal.txt"\nReadPreRunData unused', 'ReadPreRunData'),
            ('SaveProposal unused.txt\nRunMinuitOnly true', 'require MCMC'),
            ('LoadProposal "tuned proposal.txt"\nComputeHessianOnly true', 'require MCMC'),
            ('LoadProposal missing.txt', 'missing file or unsupported format/version'),
            ('LoadProposal "tuned proposal.txt"\nNChains 3', 'parameter/replica counts'),
            ('LoadProposal "tuned proposal.txt"\nPTTemperatures 1:3', 'temperature ladder'),
            ('LoadProposal "tuned proposal.txt"\nMultivariateProposal false', 'proposal type'),
            ('SaveProposal missing-directory/tuned.txt', 'ERROR'),
        ]:
            (work / 'mc.conf').write_text(settings + options + '\n')
            result = run(command, work, False)
            assert message in result.stdout, result.stdout
        # Bounds, malformed caches and legacy ROOT loading must also be rejected.
        (work / 'model.conf').write_text(model_config().replace('mtop 163.0 2.0 0', 'mtop 163.0 3.0 0'))
        (work / 'mc.conf').write_text(settings + 'LoadProposal "tuned proposal.txt"\n')
        result = run(command, work, False)
        assert 'bounds or fixed value differ' in result.stdout, result.stdout
        (work / 'model.conf').write_text(model_config())
        (work / 'broken.txt').write_text('BAT_MPI_PROPOSAL 1\n')
        (work / 'mc.conf').write_text(settings + 'LoadProposal broken.txt\n')
        result = run(command, work, False)
        assert 'invalid header' in result.stdout, result.stdout
        (work / 'mc.conf').write_text(settings + 'LoadProposal "tuned proposal.txt"\n')
        result = run([*command, '--loadRun', '--inFile', 'unused'], work, False)
        assert '--loadRun or ReadPreRunData' in result.stdout, result.stdout
        print('Invalid proposal configurations and cache mismatches terminate MPI workers correctly')


if __name__ == '__main__':
    main()
