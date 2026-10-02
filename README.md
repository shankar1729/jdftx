# JDFTx with mixed Coulomb truncation

This version of [JDFTx](https://github.com/shankar1729/jdftx) adds defect-only
mixed Coulomb truncation for localized, neutral defects on periodic slabs.
The clean host retains slab electrostatics, while the defect-induced total
charge is treated with an isolated Wigner–Seitz Coulomb kernel.

The implementation includes:

- A variational correction to the electronic SCF potential and energy, using
  a fixed clean-system reference.
- Support for Coulomb embedding and a separate isolated-defect embedding center.
- `defectDensityToVscloc`, a one-shot utility that constructs a native `Vscloc`
  and an NSCF input from clean and defect electronic densities and geometries.
- CPU/CUDA implementations and regression tests for SCF, potential construction,
  embedding, and finite-difference derivatives.

See [the method and usage guide](README-defect-coulomb.md) for input commands,
reference generation, charge conventions, embedding settings, limitations, and
validation. The original [method outline](Defect_Only_Coulomb_Truncation_O_Ag111.pdf)
is also included.

The difference charge includes both electrons and the ionic Coulomb source.
This matters when adding or removing atoms: the electronic density difference
alone need not be neutral. The current SCF method supports fixed ions and cell,
neutral vacuum calculations, and disabled spatial symmetries. Localization of
the difference charge must be checked for the chosen cell and embedding center.
A one-shot potential does not make its supplied density self-consistent.

## Build

Install the normal JDFTx dependencies (C++ compiler, CMake, GSL, FFTW,
BLAS/LAPACK, and MPI for the default configuration), then build from this root:

```bash
cmake -S jdftx -B build
cmake --build build --parallel --target jdftx ConstructDefectPotential TestDefectCoulomb
```

Use `build/jdftx` for the modified SCF calculation and
`jdftx/scripts/defectDensityToVscloc --help` for the one-shot interface.
Provide the required pseudopotentials through JDFTx's normal search paths or
`JDFTX_PSEUDO_DIR`. Python 3 is required for the one-shot wrapper and regression
scripts. CUDA build details and input examples are in the usage guide.

## Tests

With the targets above built and the test pseudopotentials available:

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
  ctest --test-dir build -R '^defect(Coulomb|Potential|Embedding)$' --output-on-failure
```

Build products and calculation output are not part of the source repository.

## Upstream and license

This repository preserves upstream history, based on JDFTx commit `52548a96`.
JDFTx and these modifications are distributed under the GNU General Public
License, version 3 or later; see [COPYING](jdftx/COPYING) and the source notices.
Retain JDFTx attribution and cite the software and methods used in calculations.
