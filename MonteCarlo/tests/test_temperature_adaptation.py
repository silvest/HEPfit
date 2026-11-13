#!/usr/bin/env python3
"""Exercise config wiring, grid-aware proposal caches and production diagnostics."""
import argparse
import json
from pathlib import Path
import sys
import tempfile
from test_pt_integration import model_config, run
from test_proposal_integration import rows, kernel, output_checks

p = argparse.ArgumentParser()
p.add_argument('--analysis', required=True)
p.add_argument('--bat-source', required=True)
p.add_argument('--mpiexec', required=True)
p.add_argument('--mpi-flags', default='')
a = p.parse_args()
command = [a.mpiexec, '-n', '3', *[f for f in a.mpi_flags.split(';') if f],
           str(Path(a.analysis).resolve()), 'model.conf', 'mc.conf']
common = '''NChains 2
PrerunMaxIter 25
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
ParallelTempering true
PTTemperatures 1:1.01:4
PTSwapInterval 2
'''
tuning = '''PTAdaptTemperatures true
PTAdaptInterval 4
PTAdaptRounds 2
PTAdaptDamping .5
PreRunPolicy staged
SettlingSweeps 7
SmoothAdaptation true
'''
with tempfile.TemporaryDirectory(prefix='hepfit-temperatures-') as directory:
    work = Path(directory)
    (work / 'model.conf').write_text(model_config())
    (work / 'regions.json').write_text(json.dumps({'regions': [
        {'name': 'low_mtop', 'conditions': [{'parameter': 'mtop', 'op': '<', 'value': 163}]},
        {'name': 'high_mtop', 'conditions': [{'parameter': 'mtop', 'op': '>=', 'value': 163}]},
    ]}))
    for joint in (False, True):
        (work / 'mc.conf').write_text(tuning + common + f'MultivariateProposal {str(joint).lower()}\nSaveProposal tuned.txt\n')
        run([*command, '--rootfile', 'fresh'], work)
        s = rows(work / 'fresh_warmup_status.csv')[0]
        assert s['temperature_updates'] == '2' and s['last_temperature_update'] == '12'
        assert s['temperatures_frozen'] == '1' and s['adaptation_sweeps'] == '25'
        history = rows(work / 'fresh_pt_temperature_adaptation.csv')
        assert {int(r['sweep']) for r in history} == {0, 4, 12}
        expected = kernel(work / 'tuned.txt')
        assert expected[-1] is not None
        final_grid = [float(t) for t in expected[-1][1]]
        assert final_grid[0] == 1 and final_grid[-1] == 4 and final_grid[1] != 1.01
        output_checks(work, 'fresh')
        for r in rows(work / 'fresh_pt_checkpoints.csv'):
            if r['phase'] == 'prerun':
                assert all(not (int(r['first_sweep']) <= t < int(r['last_sweep'])) for t in (4, 12))
        for r in rows(work / 'fresh_pt_local_moves.csv'):
            if r['phase'] in ('warmup', 'production'):
                assert float(r['temperature']) == final_grid[int(r['slot']) // 2]
        run([sys.executable, str(Path(a.bat_source) / 'diagnostics/analyze_mcmc.py'), 'fresh.root',
             '--tree', 'StandardModel_mcmc', '--regions', 'regions.json',
             '--transport-prefix', 'fresh_', '--output', 'report', '--no-plots'], work)
        report = json.loads((work / 'report/summary.json').read_text())
        assert len(rows(work / 'report/diagnostics.csv')) > 0
        weights = rows(work / 'report/region_weights.csv')
        assert len(weights) == 2 and abs(sum(float(r['weight']) for r in weights)-1) < 1.e-12
        assert report['streaming_pt']['runs'][0]['warmup_status']['temperatures_frozen'] is True
        assert report['streaming_pt']['runs'][0]['temperature_adaptation']['updates'] == 2
        # The initial intermediates can differ, but layout and endpoints must match.
        (work / 'mc.conf').write_text(common.replace('1:1.01:4', '1:2:4') +
            f'MultivariateProposal {str(joint).lower()}\nLoadProposal tuned.txt\nProposalWarmup 7\nSaveProposal restored.txt\n')
        serial = command.copy(); serial[2] = '1'
        run([*serial, '--rootfile', 'loaded'], work)
        assert kernel(work / 'restored.txt') == expected
        output_checks(work, 'loaded')
        s = rows(work / 'loaded_warmup_status.csv')[0]
        assert s['temperature_updates'] == '0' and s['proposal_loaded'] == '1'
    for options, message in [
        ('PTAdaptTemperatures true', 'PreRunPolicy staged'),
        ('PTAdaptInterval 4', 'require PTAdaptTemperatures'),
        ('PTAdaptTemperatures true\nPreRunPolicy staged\nPTAdaptInterval 1000', 'proposal-check window'),
        (tuning + 'SettlingSweeps 0', 'positive SettlingSweeps'),
        (tuning + 'LoadProposal tuned.txt', 'skips adaptation'),
        (tuning + 'PTTemperatures 1:4', 'at least three'),
    ]:
        (work / 'mc.conf').write_text(common + options + '\n')
        result = run(command, work, False)
        assert message in result.stdout, result.stdout
print('PASS: HEPfit temperature tuning, production grid, diagnostics, optimized cache reload and invalid configurations')
