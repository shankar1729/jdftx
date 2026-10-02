# Defect-only Coulomb truncation

This source tree implements the fixed-ion SCF method in
`Defect_Only_Coulomb_Truncation_O_Ag111.pdf`. At every electronic density/potential
evaluation it adds (for unembedded truncation)

```
deltaRho = (n + rhoIon) - (nClean + rhoIonClean)
Vcorr = (Kisolated - Kslab) deltaRho
EdefectCoulomb = 0.5 * integral(deltaRho * Vcorr)
```

The clean reference remains fixed. The ordinary slab calculation still treats the
periodic host; the isolated 3D Wigner–Seitz kernel acts on the defect-induced
charge. Exchange-correlation and short-range/nonlocal pseudopotentials keep their
ordinary dependence on the full electronic density. The correction enters the
SCF potential before orbital solution and mixing, and its energy enters the
reported total energy. Both density and potential mixing are supported.

JDFTx's charge convention is positive electronic number density and negative
ionic charge. Both parts must be included: subtracting only the electronic
densities after adding O leaves a charged difference. Core densities used only
for exchange-correlation are excluded. The ionic source is the **bare** charge
that generates the long-range part of `Vlocps`, independently of `ion-width`.
For unembedded truncation, both source and correction potential are projected onto the real FFT grid
with `I`/`J`; this removes unrepresentable point-ion Nyquist components and makes
the energy, real-space potential, and ionic derivative consistent on that grid.

## Clean reference

First relax and converge the clean slab with the lattice, pseudopotentials,
cutoffs, spin treatment, smearing, and FFT grid intended for the defect run.
Use slab truncation, with or without embedding. Preserve the Coulomb settings
that generated the density. In the final clean calculation add:

```text
coulomb-interaction Slab 001
symmetries none
# Keep your lattice, species, ions, k points, cutoffs and smearing here.
dump-name clean.$VAR
dump End State ElecDensity DefectReference
```

`clean.defectReference` contains the converged spin-summed electronic density,
the ionic Coulomb source at the clean geometry, and metadata. There is no need to
reconstruct clean ions in the defect calculation, or to supply an atomic-density
approximation. The stored reference can represent either the relaxed pristine
surface or a separately converged frozen-host geometry.

To export an already converged calculation with saved wavefunctions, use its
original input/settings and ionic positions, plus:

```text
initial-state clean.$VAR
dump-only
dump-name clean.$VAR
dump End DefectReference
```

Remove an explicit `wavefunction` command when using `initial-state`. A raw
`clean.n` file alone is insufficient: it lacks the reference ionic contribution
and the cell/grid metadata. The new reference dump packages the needed inputs.

## Defect SCF

Use the same cell and grid, and add the defect to your ionic input. Keep the
reference and current geometries in the same coordinate origin. For O/Ag(111),
keep Ag pseudopotentials identical to the clean run and add the O species/ion.

```text
coulomb-interaction Slab 001
symmetries none
electronic-scf
defect-coulomb clean.defectReference

dump-name defect.$VAR
dump End State ElecDensity Vscloc Ecomponents Dvac DefectCoulomb
```

The optional second argument of `defect-coulomb` is an absolute neutrality
tolerance in electrons, default `1e-5`. Each density evaluation checks the
integral of `deltaRho`; a failure aborts without silently inserting a charge
background. Each SCF report includes `Qdefect` and `Ecorr`.

At startup the reader checks the file version and payload length, cell, grid,
cutoffs, slab direction, spin treatment, smearing, finite density values,
reference neutrality, and the valence charge and content fingerprint of each
host pseudopotential. Host atoms may move between the separately prepared clean
and defect geometries; the difference includes their charge displacement.
Added defect species are permitted. Use the same k-point sampling and functional
in both runs; these are workflow requirements, not checked by the reference
header. Verify electronic convergence before exporting the reference.

## Outputs and supported scope

`dump End DefectCoulomb` writes `defect.deltaRho` and `defect.Vcorr` as standard
JDFTx real-space binary grids. `Dvac` and `Dtot` include the correction as well.
Use `potential-subtraction no` in both clean and defect runs when comparing the
unsmoothed electrostatic potentials directly. The energy component is named
`EdefectCoulomb`.

This implementation supports neutral, fixed-electron-number, vacuum SCF
calculations with fixed ions and lattice, slab truncation with optional embedding, and
`symmetries none`. It rejects ionic/lattice optimization, molecular dynamics,
stress, vibrations, variational perturbation calculations, fluids, external
charges/fields, and fixed chemical potential. The explicit ionic correction is
included in force evaluation and has a finite-difference regression check;
geometry optimization remains disabled pending surface-system force validation.
Spin-unpolarized, collinear and noncollinear densities are summed over their
diagonal charge components; regression checks cover unpolarized and collinear
cases. The existing exchange kernel is unchanged.

The density difference must be localized within the support allowed by isolated
WS truncation. In the unembedded scheme, charge separations must fit inside the
truncation cell (the usual half-cell confinement requirement). Plot
`deltaRho` and `Vcorr` and check their lateral behavior. The code checks
neutrality, but does not certify localization. See the
[JDFTx Coulomb-interaction documentation](https://jdftx.org/CommandCoulombInteraction.html)
for the truncation confinement requirement.

The regression system is a small neutral He host with a localized He addition;
it validates the implementation, rather than the physical approximation for
O/Ag(111). The PDF's
O/Ag force checks, lateral-size comparisons, and local Wannier-Hamiltonian
localization study remain necessary before interpreting production results.

## Build and reproduce the checks

A CPU binary is built in `build/jdftx`:

```bash
cmake -S jdftx -B build -DEnableCUDA=OFF -DEnableMPI=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target jdftx TestDefectCoulomb -j 12
export JDFTX_PSEUDO_DIR=/path/to/jdftx/pseudopotentials
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=2
ctest --test-dir build -R '^defectCoulomb$' --output-on-failure
```

The integration test has 23 checks: clean recovery, a nonzero converged defect
correction, charge neutrality, real-space energy/potential consistency,
collinear spin recovery, a hexagonal cell with potential mixing and smearing, independence from ionic plotting width, and rejection
of missing, truncated, incompatible, or charged reference files. Like other
JDFTx tests, successful `.out` files are reused; clear the test run's outputs
when rerunning calculations after source changes.

The additional native test checks density energy derivatives, insertion into
the KS potential and total energy, the frozen reference, preservation of the
ionic source, and the explicit ionic derivative through the pseudopotential
force path:

```bash
cd build/test/defectCoulomb
export SRCDIR=../../../jdftx/test/defectCoulomb
../../aux/TestDefectCoulomb -i "$SRCDIR/recovery.in" -o derivatives.out
```

For MPI, set `JDFTX_LAUNCH='mpirun -n 2'` before rerunning the integration test.
GPU builds use the same field operators; enable CUDA and build `jdftx_gpu` and
`TestDefectCoulomb_gpu`. For the RTX 4090s in this workspace, the GPU build is
configured in `build-gpu` with CUDA architecture 89.

The main implementation is in `jdftx/electronic/DefectCoulomb.{h,cpp}`, with input
handling in `jdftx/commands/defect_coulomb.cpp` and integration into `ElecVars`,
`IonInfo`, `SCF`, and `Dump`.

Validation completed in this workspace: all 23 integration checks passed on
CPU (including two MPI processes) and CUDA. Native density, KS-potential, energy,
and ionic finite-difference checks passed on CPU/MPI and CUDA; orthorhombic and
hexagonal cells were checked. Existing `openShell`, `graphene`, and `checkpoint`
regressions passed. These checks preceded the embedded tests described below.

## One-shot potential construction for NSCF

`jdftx/scripts/defectDensityToVscloc` takes the saved clean and defect electronic
densities and constructs the corrected **full local KS potential** once. It uses
JDFTx itself to evaluate slab Hartree, local pseudopotentials, XC (including
partial core corrections), and the defect correction. It preserves JDFTx's
internal `Vscloc` normalization and spin channel naming, so the result can be
read directly with `fix-electron-potential`.

Run from the source-tree root:

```bash
./jdftx/scripts/defectDensityToVscloc \
  --clean-input /path/to/clean.in \
  --defect-input /path/to/defect.in \
  --clean-density /path/to/clean.n \
  --defect-density /path/to/defect.n \
  --output-prefix /path/to/nscf/defect-mixed
```

The two input files supply the lattice, grid, ionic geometries,
pseudopotentials, spin treatment and functional. Their ionic positions must
match the saved densities, including any prior relaxation. The script expands
includes and environment variables through JDFTx's input reader, and runs each
stage from that original input file's directory. It ignores the original
SCF/minimization, wavefunction/restart and dump controls; no orbitals are
allocated or solved, and no density iteration occurs. It retains the original
FFT grid, then disables symmetry averaging of the corrected potential.

For wildcard species such as `SG15/$ID_ONCV_PBE.upf`, the script uses the usual
JDFTx search paths. If `JDFTX_PSEUDO_DIR` is unset, it also looks for the library
beside `jdftx` or `jdftx_gpu` on your `PATH` (including the installed library
under `share/jdftx/pseudopotentials`). Use `--pseudo-dir /path/to/pseudopotentials`
to set the additional search directory explicitly. The selected directory is
printed, and the generated NSCF input records absolute paths to the actual
pseudopotential files used for potential construction.

For collinear spin, supply quoted filename patterns instead:

```bash
--clean-density '/path/to/clean.$VAR' \
--defect-density '/path/to/defect.$VAR'
```

JDFTx reads `n_up` and `n_dn`; noncollinear calculations additionally need `n_re`
and `n_im`. A raw `.n_up` filename is also accepted and converted to the pattern.
The raw files must be standard JDFTx `ElecDensity` dumps on matching grids.
Their byte counts are checked, but raw densities contain no cell metadata:
the supplied input files must actually describe the calculations that produced
them. The clean ionic charge is reconstructed at the clean geometry and stored
in the generated reference, rather than inferred from the defect geometry.

Outputs with this prefix are:

- `.Vscloc` (or `.Vscloc_up`, `.Vscloc_dn`, and noncollinear off-diagonal channels)
- `.Vcorr` and `.deltaRho` for diagnostics (also `.deltaRhoEffective` with embedding)
- `.defectReference`, containing the clean density and reconstructed clean ions
- `.clean.in`/`.defect.in` recipes and `.clean.out`/`.defect.out` logs
- `.nscf.in`, a fixed-potential NSCF input with the exact FFT grid

Run the generated NSCF input from the original defect-input directory. It
already contains:

```text
fix-electron-potential /path/to/nscf/defect-mixed.$VAR
```

Adjust k points and `elec-n-bands` there for the desired NSCF calculation. Keep
the cell, density grid, spin treatment, ionic geometry and pseudopotentials
consistent. Do not add `defect-coulomb` or `electronic-scf` to this NSCF input;
the saved potential already includes the correction. The ordinary nonlocal
pseudopotential operator is provided by the retained species settings.

This path accepts internal LDA/GGA functionals, neutral vacuum slabs
and slab truncation with optional embedding. It rejects hybrids, meta-GGAs, orbital-dependent
functionals and DFT+U, which require data beyond two electronic densities.
The same neutrality and reference compatibility checks apply as in SCF mode.
Use `--charge-tolerance` to change the default `1e-5` electron tolerance, and
`--overwrite` to reuse an output prefix. Source densities are preserved.

Add `--gpu` to use CUDA, or `--launcher 'mpirun -n 2'` to launch the backend with
MPI. `--engine` selects a specific `ConstructDefectPotential` executable. Both
CPU and GPU backends have been built in this workspace; to rebuild:

```bash
cmake -S jdftx -B build
cmake --build build --target ConstructDefectPotential
# For an already CUDA-configured build:
cmake -S jdftx -B build-gpu
cmake --build build-gpu --target ConstructDefectPotential_gpu
```

The potential is evaluated at the supplied defect density. If that density came
from a conventional periodic calculation, the one-shot correction does **not**
make it self-consistent under mixed electrostatics; use the corrected SCF path
when that electronic response is required.

`ctest --test-dir build -R '^defectPotential$' --output-on-failure` runs the
required density fixture and the end-to-end one-shot regression. It compares
the complete potential with an ordinary fixed-density JDFTx calculation having
the identical `Vcorr` added externally, covers spin and hexagonal/smearing cases,
checks source-density preservation and overwrite protection, and performs an
actual NSCF solve using the generated potential.

The one-shot regression passed on CPU and CUDA in this workspace, including
actual fixed-potential NSCF solves. The two-process MPI script launch was also
validated.

## Embedded mixed Coulomb truncation

Keep `coulomb-truncation-embed` in the clean and defect inputs. Both the normal
embedded slab Hamiltonian and the mixed correction are supported in SCF and in
the one-shot utility. The isolated operator doubles its auxiliary Coulomb grid
in all three directions; the slab operator doubles only its truncated direction.
The wavefunction cell and density grid stay the same.

Version-2 clean references record the embedding mode, center in lattice
coordinates, and ionic margin, and preserve the complete bare ionic Fourier
source. These settings must agree between reference generation and the defect
run. Centers differing by integer lattice vectors are equivalent. Version-1
references remain readable for unembedded calculations; regenerate them for
embedded calculations. Raw density files contain no record of the Coulomb
settings: supply the physical settings that actually generated the densities.
Removing embedding from a reconstruction input changes its underlying potential.

With embedding, let `S` be JDFTx's Gaussian range-separation filter (its width is
taken from the slab Coulomb operator), and `D=Kisolated-Kslab` on smooth sources.
The common ionic short-range kernel cancels between the two operators. The
numerical correction is

```text
deltaN = n - nClean
deltaIon = rhoIonBare - rhoIonBareClean
deltaRhoEffective = deltaN + S deltaIon
Vcorr = D deltaRhoEffective
EdefectCoulomb = 0.5 * integral(deltaRhoEffective * Vcorr)
ionicGradientPotential = S Vcorr
```

This is equivalent to applying the difference of JDFTx's `PointChargeRight`
operators to the ionic difference, and smooth operators to the electron
difference. The ionic derivative uses the adjoint Gaussian filter so the energy
and forces remain consistent. It does not change the local/nonlocal
pseudopotentials or add another ionic short-range contribution. Range separation
is independent of the `ion-width` used for fluid/plotting charges.

`dump End DefectCoulomb` additionally writes `deltaRhoEffective` when embedding
is enabled. Use this field with `Vcorr` for the real-space energy check;
`deltaRho` retains the unsmoothed total difference for charge diagnostics.

The isolated operator inherits all three components of the slab embedding
center. Its lateral coordinates matter: place the isolated center near the
localized defect so the difference does not straddle a cut. You can set a
separate isolated center without changing the slab center:

```text
coords-type Cartesian
coulomb-interaction Slab 001
coulomb-truncation-embed 0 0 18.105968516010567
defect-coulomb clean.defectReference
defect-coulomb-center 2.771410394108804 8.001439334819281 18.105968516010567
```

The example uses Cartesian coordinates in bohr. Lattice coordinates are also
supported. `defect-coulomb-center` requires embedding. The one-shot utility
preserves this center if it is present in the defect input. Use
`--defect-center C0 C1 C2` to override it explicitly; its
coordinates follow the defect input's `coords-type`. The generated NSCF input
retains the slab embedding command; the saved potential already includes the
correction.

The confinement condition now refers to the doubled Coulomb box. Charge must
remain localized away from the cuts of the original cell centered on the
specified embedding center. Check lateral localization of the defect-induced
charge, particularly for a small metallic surface supercell. Embedding does not
establish localization by itself.

The `defectEmbedding` regression covers orthorhombic and hexagonal embedded SCF,
clean reference recovery, nonzero corrections, one-shot reconstruction against
ordinary fixed-density JDFTx, finite-difference density and ionic derivatives,
equivalence to `PointChargeRight`, and rejection of inconsistent embedding
modes, centers, and ionic margins. Build `TestDefectCoulomb` (or its GPU target)
before running:

```bash
ctest --test-dir build -R '^defectEmbedding$' --output-on-failure
```
