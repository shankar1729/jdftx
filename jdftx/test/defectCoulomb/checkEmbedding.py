#!/usr/bin/env python3
"""Embedded SCF, point-ion derivatives, one-shot potentials, and center validation."""
import array
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

build, source, script = map(lambda p: Path(p).resolve(), sys.argv[1:4])
work = build / 'test/defectEmbedding'
work.mkdir(parents=True, exist_ok=True)
suffix = os.environ.get('JDFTX_SUFFIX', '')
env = dict(os.environ, SRCDIR=str(source))
launcher = shlex.split(env.get('JDFTX_LAUNCH', ''))
jdftx = build / ('jdftx' + suffix)
engine = build / ('ConstructDefectPotential' + suffix)
derivatives = build / 'aux' / ('TestDefectCoulomb' + suffix)


def run(command, cwd=work, succeeds=True):
    result = subprocess.run(list(map(str, command)), cwd=cwd, env=env,
                            capture_output=True, text=True, timeout=180)
    assert (result.returncode == 0) == succeeds, result.stdout + result.stderr
    return result


def native(stem, dry=False, executable=jdftx, succeeds=True):
    return run(launcher + [executable, '-n' if dry else '-d', '-i', work / (stem + '.in'),
                          '-o', work / (stem + '.out')], succeeds=succeeds)


def values(path):
    x = array.array('d')
    x.frombytes(Path(path).read_bytes())
    if sys.byteorder != 'little':
        x.byteswap()
    assert all(map(math.isfinite, x))
    return x


def component(stem, name):
    for ending in ('.Ecomponents', '.ecomponents'):
        path = work / (stem + ending)
        if path.exists():
            return float(re.search(r'^\s*' + name + r'\s*=\s*(\S+)', path.read_text(), re.M)[1])
    raise AssertionError('Missing energy components: ' + stem)


def oneshot(clean, defect, output, center=False):
    command = [sys.executable, script, '--engine', engine,
               '--clean-input', work / (clean + '.in'), '--defect-input', work / (defect + '.in'),
               '--clean-density', work / (clean + '.n'), '--defect-density', work / (defect + '.n'),
               '--output-prefix', work / output, '--overwrite']
    if launcher:
        command += ['--launcher', shlex.join(launcher)]
    if center:
        command += ['--defect-center', '1.11', '0.33', '1.40']
    run(command)


for tag, common in [('ortho', 'common.in'), ('hex', 'hexCommon.in')]:
    physical = f'include {source / common}\ncoords-type Cartesian\ncoulomb-truncation-embed 0.63 -0.48 1.17\n'
    clean, recovery, defect = [tag + x for x in ('Clean', 'Recovery', 'Defect')]
    (work / (clean + '.in')).write_text(physical + f'ion He 0 0 0 0\ndump-name {clean}.$VAR\n'
        'dump End State ElecDensity DefectReference Ecomponents\n')
    native(clean)
    (work / (recovery + '.in')).write_text(physical + f'ion He 0 0 0 0\ninitial-state {clean}.$VAR\n'
        f'defect-coulomb {clean}.defectReference\ndump-only\ndump-name {recovery}.$VAR\n'
        'dump End DefectCoulomb Ecomponents\n')
    native(recovery)
    assert abs(component(recovery, 'EdefectCoulomb')) < 1e-12
    assert abs(component(recovery, 'Etot') - component(clean, 'Etot')) < 1e-9
    assert max(map(abs, values(work / (recovery + '.Vcorr')))) < 1e-10
    (work / (defect + '.in')).write_text(physical + 'ion He 0 0 0 0\nion He 0 0 3.2 0\n'
        f'defect-coulomb {clean}.defectReference\ndefect-coulomb-center 1.11 0.33 1.40\n'
        f'dump-name {defect}.$VAR\ndump End ElecDensity Vscloc DefectCoulomb Ecomponents\n')
    native(defect)
    assert 'SCF: Converged' in (work / (defect + '.out')).read_text()
    reference = (work / (clean + '.defectReference')).read_bytes()
    assert reference.startswith(b'JDFTX_DEFECT_REFERENCE 2\n')
    header = reference.partition(b'DATA\n')[0].decode().splitlines()
    nr = math.prod(map(int, header[1].split()))
    volume = 16 * 16 * 24 * (math.sqrt(3)/2 if tag == 'hex' else 1)
    dv = volume/nr
    rho = values(work / (defect + '.deltaRhoEffective'))
    v = values(work / (defect + '.Vcorr'))
    assert abs(sum(rho)*dv) < 1e-8
    assert abs(.5*sum(x*y for x, y in zip(rho, v))*dv-component(defect, 'EdefectCoulomb')) < 1e-9
    assert max(map(abs, v)) > 1e-6
    oneshot(clean, clean, tag + 'OneRecovery', center=True)
    assert max(map(abs, values(work / (tag + 'OneRecovery.Vcorr')))) < 1e-10
    # Preserve a center declared in the original corrected-SCF input without
    # requiring the user to repeat it as a script option.
    oneshot(clean, defect, tag + 'OneDefect')
    reconstructed_corr = values(work / (tag + 'OneDefect.Vcorr'))
    assert max(abs(x-y) for x, y in zip(reconstructed_corr, v)) < 1e-9
    reconstructed = values(work / (tag + 'OneDefect.Vscloc'))
    # Independently reconstruct the complete potential through ordinary fixed-density
    # JDFTx and its embedded point-ion path, with Vcorr supplied as Vexternal.
    recipe = (work / (tag + 'OneDefect.nscf.in')).read_text()
    base = '\n'.join(line for line in recipe.splitlines()
                     if not line.startswith(('fix-electron-potential', 'dump-name', 'dump End')))
    direct = tag + 'Direct'
    (work / (direct + '.in')).write_text(base + f'\nfix-electron-density {work / defect}.$VAR\n'
        f'Vexternal {work / tag}OneDefect.Vcorr\ndump-name {direct}.$VAR\ndump End Vscloc\n')
    native(direct)
    expected = values(work / (direct + '.Vscloc'))
    errors = [(x-y)/dv for x, y in zip(reconstructed, expected)]
    assert max(map(abs, errors)) < 1e-6 and math.sqrt(sum(x*x for x in errors)/nr) < 1e-8
    assert derivatives.is_file(), 'Build TestDefectCoulomb' + suffix + ' for embedded derivative checks.'
    derivative = tag + 'Derivatives'
    (work / (derivative + '.in')).write_text(physical + f'ion He 0 0 0 0\ninitial-state {clean}.$VAR\n'
        f'defect-coulomb {clean}.defectReference\ndefect-coulomb-center 1.11 0.33 1.40\n')
    native(derivative, executable=derivatives)
    assert 'All defect Coulomb derivative tests passed' in (work / (derivative + '.out')).read_text()
    assert 'Embedded PointChargeRight equivalence' in (work / (derivative + '.out')).read_text()
    default_derivative = tag + 'DefaultCenterDerivatives'
    (work / (default_derivative + '.in')).write_text(physical + f'ion He 0 0 0 0\ninitial-state {clean}.$VAR\n'
        f'defect-coulomb {clean}.defectReference\n')
    native(default_derivative, executable=derivatives)
    assert 'All defect Coulomb derivative tests passed' in (work / (default_derivative + '.out')).read_text()
    print('Passed embedded SCF, clean recovery, one-shot potential, and density/ionic derivatives:', tag, flush=True)

    # Reference metadata prevents silently switching the Hamiltonian during reconstruction.
    for name, changed, message in [
        ('Mode', physical.replace('coulomb-truncation-embed 0.63 -0.48 1.17\n', ''), 'embedding mode mismatch'),
        ('Center', physical.replace('0.63 -0.48 1.17', '0.63 -0.48 2.17'), 'embedding center mismatch'),
        ('Margin', physical + 'coulomb-truncation-ion-margin 4\n', 'embedding ionic margin mismatch'),
    ]:
        stem = tag + 'Wrong' + name
        (work / (stem + '.in')).write_text(changed + f'ion He 0 0 0 0\ndefect-coulomb {clean}.defectReference\n')
        native(stem, dry=True, succeeds=False)
        assert message in (work / (stem + '.out')).read_text()
    print('Passed embedding mode, center and ionic-margin compatibility checks:', tag, flush=True)
