"""Integration checks: a converged clean restart and neutral localized addition."""
import array
import math
import re
from pathlib import Path


def component(run, name):
    text = Path(run + '.Ecomponents').read_text() if Path(run + '.Ecomponents').exists() else Path(run + '.ecomponents').read_text()
    return float(re.search(r'^\s*' + name + r'\s*=\s*(\S+)', text, re.M)[1])


def values(path):
    data = array.array('d')
    data.frombytes(Path(path).read_bytes())
    import sys
    if sys.byteorder != 'little':
        data.byteswap()
    assert all(math.isfinite(x) for x in data)
    return data


def check(value, expected, tolerance, label):
    assert math.isfinite(value)
    print(value, expected, tolerance, label)


print(23)
check(component('recovery', 'EdefectCoulomb'), 0, 1e-12, 'Clean recovery correction energy')
check(component('recovery', 'Etot') - component('clean', 'Etot'), 0, 1e-9, 'Clean recovery total energy')
check(max(map(abs, values('recovery.Vcorr'))), 0, 1e-10, 'Clean recovery correction potential')
check(max(map(abs, values('recovery.deltaRho'))), 0, 1e-10, 'Clean recovery density difference')
rho = values('defect.deltaRho')
v = values('defect.Vcorr')
volume = 16 * 16 * 24
check(sum(rho) * volume / len(rho), 0, 1e-8, 'Defect charge neutrality')
check(0.5 * sum(x * y for x, y in zip(rho, v)) * volume / len(rho)
      - component('defect', 'EdefectCoulomb'), 0, 1e-9, 'Correction energy matches dumped fields')
check(float(max(map(abs, v)) > 1e-6), 1, 0, 'Localized addition produces nonzero potential')
check(float('SCF: Converged' in Path('defect.out').read_text()), 1, 0, 'Corrected SCF converges')

check(component('spinRecovery', 'EdefectCoulomb'), 0, 1e-12, 'Spin reference recovery correction energy')
check(max(map(abs, values('spinRecovery.Vcorr'))), 0, 1e-10, 'Spin reference recovery potential')
check(component('widthRecovery', 'EdefectCoulomb'), 0, 1e-12, 'Ionic plotting width does not affect correction')
check(component('widthRecovery', 'Etot') - component('clean', 'Etot'), 0, 1e-9, 'Ionic plotting width does not affect energy')

# The backend must reject incompatible/corrupt data before any SCF step.
import os
import subprocess
exe = Path.cwd().parents[1] / ('jdftx' + os.environ.get('JDFTX_SUFFIX', ''))
header, marker, payload = Path('clean.defectReference').read_bytes().partition(b'DATA\n')
lines = header.decode().splitlines()

def rejected(name, reference, expected):
    Path(name + '.in').write_text(
        'include ' + os.environ['SRCDIR'] + '/common.in\n'
        'ion He 0 0 0 0\ndefect-coulomb ' + reference + '\n')
    result = subprocess.run([str(exe), '-n', '-i', name + '.in', '-o', name + '.out'],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
    log = Path(name + '.out').read_text()
    check(float(result.returncode != 0 and expected in log), 1, 0, name + ' rejected')

rejected('missingReference', 'missing.defectReference', 'Cannot read defect Coulomb reference')
Path('truncated.defectReference').write_bytes(header + marker + payload[:-8])
rejected('truncatedReference', 'truncated.defectReference', 'density payload length')
for name, index, replacement, message in [
    ('wrongGrid', 1, '61 60 90', 'FFT grid mismatch'),
    ('wrongLattice', 2, '17', 'lattice mismatch'),
    ('wrongPseudo', 13, 'He 2 0000000000000000', 'host pseudopotential mismatch'),
    ('wrongSpin', 11, '15 60 2 1 0 0.001', 'slab direction/spin mismatch'),
]:
    changed = lines.copy()
    changed[index] = replacement
    Path(name + '.defectReference').write_bytes(('\n'.join(changed) + '\n').encode() + marker + payload)
    rejected(name, name + '.defectReference', message)
import struct
charged = bytearray(payload)
struct.pack_into('<d', charged, 0, struct.unpack_from('<d', charged, 0)[0] + 1)
Path('charged.defectReference').write_bytes(header + marker + charged)
rejected('chargedReference', 'charged.defectReference', 'non-neutral clean reference')

check(component('hexRecovery', 'EdefectCoulomb'), 0, 1e-12, 'Hexagonal clean recovery energy')
check(max(map(abs, values('hexRecovery.Vcorr'))), 0, 1e-10, 'Hexagonal clean recovery potential')
rho = values('hexDefect.deltaRho')
v = values('hexDefect.Vcorr')
volume = 16 * 16 * 24 * math.sqrt(3) / 2
check(0.5 * sum(x * y for x, y in zip(rho, v)) * volume / len(rho)
      - component('hexDefect', 'EdefectCoulomb'), 0, 1e-9, 'Hexagonal correction energy consistency')
check(float('SCF: Converged' in Path('hexDefect.out').read_text()), 1, 0,
      'Hexagonal potential mixing with smearing converges')
