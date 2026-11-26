#!/usr/bin/env python3
"""Sample a HEPfit parameter with two separated flat intervals in one MPI fit."""
import argparse
from pathlib import Path
import tempfile

import numpy as np
import uproot

from test_pt_integration import model_config, run


def check_samples(mass):
    left = (mass >= 161) & (mass <= 163.5)
    right = (mass >= 165) & (mass <= 167)
    assert np.all(left | right), mass[~(left | right)]
    assert np.any(left) and np.any(right), 'Both intervals should be visited'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--analysis', required=True)
    parser.add_argument('--mpiexec', required=True)
    parser.add_argument('--mpi-flags', default='')
    args = parser.parse_args()
    analysis = str(Path(args.analysis).resolve())
    flags = [flag for flag in args.mpi_flags.split(';') if flag]

    with tempfile.TemporaryDirectory(prefix='hepfit-disjoint-') as directory:
        work = Path(directory)
        # mtop is the only varied parameter. The envelope midpoint lies in the gap.
        (work / 'model.conf').write_text(model_config(mtop_prior='flat_intervals 161 163.5 165 167'))
        settings = '''NChains 2
PrerunMaxIter 100
NIterationsUpdateMax 50
Iterations 400
Seed 9876
MCMCInitialPosition Center
WriteChain true
FindModeWithMinuit false
CalculateNormalization false
PrintAllMarginalized false
PrintCorrelationMatrix false
PrintKnowledgeUpdatePlots false
PrintParameterPlot false
PrintTrianglePlot false
'''
        (work / 'mc.conf').write_text('ParallelTempering true\nPTTemperatures 1:2:4\nPTSwapInterval 5\n' + settings)
        command = [args.mpiexec, '-n', '3', *flags, analysis, 'model.conf', 'mc.conf']
        warning = 'Parallel tempering is disabled'
        result = run(command, work)
        assert warning not in result.stdout, result.stdout
        with uproot.open(work / 'MCout.root') as output:
            samples = output['StandardModel_mcmc'].arrays(library='np')
        assert set(samples['Phase']) == {1}
        np.testing.assert_allclose(samples['LogPrior'], -np.log(4.5))
        for chain in range(2):
            mass = samples['mtop'][samples['Chain'] == chain]
            assert len(mass) == 400
            check_samples(mass)

        result = run([analysis, 'model.conf', '--noMC', '--it=50', '--output_folder', 'events'], work)
        assert warning not in result.stdout, result.stdout
        events = np.loadtxt(work / 'GeneratedEvents/events/Parameters/mtop.txt')
        assert len(events) == 51 and events[0] == 163
        check_samples(events[1:])

        model = work / 'model.conf'
        model.write_text(model.read_text().replace('ModelParameter mHl 125.09 0 0',
            'ModelParameter mHl 125.09 flat_intervals 124 126 127 128'))
        for option in ('ParallelTempering false\n', ''):
            (work / 'mc.conf').write_text(option + settings)
            result = run(command, work)
            warnings = [line for line in result.stdout.splitlines() if warning in line]
            assert len(warnings) == 1, result.stdout
            assert 'mtop' in warnings[0] and 'mHl' in warnings[0], warnings[0]
            assert (work / 'log.txt').read_text().count(warning) == 1

        (work / 'model.conf').write_text(model_config())
        result = run(command, work)
        assert warning not in result.stdout, result.stdout


if __name__ == '__main__':
    main()
