#!/usr/bin/env python3
"""Verify one-shot construction and use of its native NSCF potential."""
import array
import hashlib
import math
import os
from pathlib import Path
import subprocess
import sys

build = Path(sys.argv[1]).resolve()
source = Path(sys.argv[2]).resolve()
script = Path(sys.argv[3]).resolve()
work = build / 'test' / 'defectCoulomb'
suffix = os.environ.get('JDFTX_SUFFIX', '')
env = dict(os.environ, SRCDIR=str(source))
engine = build / ('ConstructDefectPotential' + suffix)
jdftx = build / ('jdftx' + suffix)


def run(command, cwd=work, succeeds=True, environment=None):
    result = subprocess.run(list(map(str, command)), cwd=cwd, env=env if environment is None else environment,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=90)
    assert (result.returncode == 0) == succeeds, result.stdout + result.stderr
    return result


def data(path):
    values = array.array('d')
    values.frombytes(Path(path).read_bytes())
    if sys.byteorder != 'little':
        values.byteswap()
    return values


def close_files(a, b, tolerance=1e-10):
    a, b = data(a), data(b)
    assert len(a) == len(b)
    error = max(abs(x-y) for x, y in zip(a, b))
    assert error <= tolerance, error



def close_potentials(a, b, reference):
    # Native Vscloc stores the density gradient, dV times the physical potential.
    lines = Path(reference).read_bytes().partition(b'DATA\n')[0].decode().splitlines()
    grid = list(map(int, lines[1].split()))
    r = list(map(float, lines[2:11]))
    volume = abs(r[0]*(r[4]*r[8]-r[5]*r[7]) - r[1]*(r[3]*r[8]-r[5]*r[6]) + r[2]*(r[3]*r[7]-r[4]*r[6]))
    dv = volume / math.prod(grid)
    a, b = data(a), data(b)
    assert len(a) == len(b) == math.prod(grid)
    errors = [(x-y)/dv for x, y in zip(a, b)]
    maximum = max(map(abs, errors))
    rms = math.sqrt(sum(x*x for x in errors)/len(errors))
    # GGA derivatives amplify different FFT-plan rounding near the vacuum density cutoff.
    assert maximum < 1e-6 and rms < 1e-8, (maximum, rms)


def command(clean, defect, prefix):
    return [sys.executable, script, '--clean-input', source / (clean + '.in'),
            '--defect-input', source / (defect + '.in'),
            '--clean-density', str(work / clean) + '.$VAR',
            '--defect-density', str(work / defect) + '.$VAR',
            '--output-prefix', work / prefix, '--engine', engine]


assert work.exists(), 'Run the defectCoulomb fixture first.'
source_paths = [work / name for name in ('clean.n', 'defect.n', 'spinClean.n_up', 'spinClean.n_dn', 'hexDefect.deltaRho')]
hashes = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}

for clean, defect, prefix in [('clean', 'defect', 'potentialTest'),
                              ('clean', 'clean', 'potentialRecovery'),
                              ('spinClean', 'spinClean', 'potentialSpin'),
                              ('hexClean', 'hexDefect', 'potentialHex')]:
    run(command(clean, defect, prefix) + ['--overwrite'])
    log = (work / (prefix + '.defect.out')).read_text()
    assert 'Skipped wave function initialization' in log
    assert 'Constructed Vscloc once' in log
    assert 'SCF: Iter' not in log and 'BandDavidson: Iter' not in log
    if defect in ('defect', 'hexDefect'):
        close_files(work / (prefix + '.Vcorr'), work / (defect + '.Vcorr'))
        close_files(work / (prefix + '.deltaRho'), work / (defect + '.deltaRho'))
    else:
        assert max(map(abs, data(work / (prefix + '.Vcorr')))) < 1e-10
    # Compare the entire potential against the standard fixed-density calculation
    # with the same Vcorr supplied as an external potential.
    template = (work / (prefix + '.nscf.in')).read_text()
    physical = '\n'.join(line for line in template.splitlines()
                         if not line.startswith(('fix-electron-potential', 'dump-name', 'dump End')))
    direct_prefix = work / (prefix + '.direct')
    direct = physical + f'\nfix-electron-density {work / defect}.$VAR\n'
    direct += f'Vexternal {work / prefix}.Vcorr\n'
    direct += f'dump-name {direct_prefix}.$VAR\ndump End Vscloc ElecDensity\n'
    direct_path = work / (prefix + '.direct.in')
    direct_path.write_text(direct)
    run([jdftx, '-d', '-i', direct_path, '-o', work / (prefix + '.direct.out')], cwd=source)
    channels = ('Vscloc_up', 'Vscloc_dn') if clean == 'spinClean' else ('Vscloc',)
    for channel in channels:
        close_potentials(str(work / prefix) + '.' + channel, str(direct_prefix) + '.' + channel,
                         str(work / prefix) + '.defectReference')
    density_channels = ('n_up', 'n_dn') if clean == 'spinClean' else ('n',)
    for channel in density_channels:
        close_files(str(work / defect) + '.' + channel, str(direct_prefix) + '.' + channel, 0)
    print('Passed density-only construction and complete potential comparison:', prefix)

# Exercise an actual NSCF solve through fix-electron-potential using the generated input.
nscf = work / 'potentialTest.nscf.in'
run([jdftx, '-d', '-i', nscf, '-o', work / 'potentialTest.nscf.out'], cwd=source)
log = (work / 'potentialTest.nscf.out').read_text()
assert 'Done!' in log and 'SCF: Iter' not in log
assert (work / 'potentialTest.nscf.eigenvals').is_file()
print('Passed NSCF orbital solution at the generated fixed potential')

# A newly built utility can lack the library available to the original SCF binary.
# Exercise the explicit search directory without inheriting the environment setting.
if env.get('JDFTX_PSEUDO_DIR'):
    without_pseudo = {key: value for key, value in env.items() if key != 'JDFTX_PSEUDO_DIR'}
    run(command('clean', 'clean', 'potentialPseudoPath') +
        ['--pseudo-dir', env['JDFTX_PSEUDO_DIR'], '--overwrite'], environment=without_pseudo)
    template = (work / 'potentialPseudoPath.nscf.in').read_text()
    species = [line.split(maxsplit=1)[1] for line in template.splitlines() if line.startswith('ion-species ')]
    assert species and all(Path(path).is_absolute() and Path(path).is_file() for path in species)
    run([jdftx, '-d', '-i', work / 'potentialPseudoPath.nscf.in',
         '-o', work / 'potentialPseudoPath.nscf.out'], cwd=source, environment=without_pseudo)
    print('Passed explicit pseudopotential search and NSCF without the original search environment')

# Verify collision handling and protection of source files.
run(command('clean', 'defect', 'potentialTest'), succeeds=False)
for path, digest in hashes.items():
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
print('Passed overwrite guard and preservation of input densities')
