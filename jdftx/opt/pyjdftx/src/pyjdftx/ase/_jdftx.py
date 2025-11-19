from typing import Optional, Union
from functools import cache

import numpy as np
from ase import Atoms
from ase.units import Bohr, Hartree, eV
from ase.calculators.calculator import Calculator, all_changes, kptdensity2monkhorstpack

from .. import JDFTxWrapper


@cache
def valid_commands() -> set[str]:
    """Commands that can be specified in the jdftx input-file format while
    initializing the JDFTx calculator. This excludes `reserved_commands`
    that should be specified using the ASE standard syntax instead."""
    result = set(cmd.lower() for cmd in JDFTxWrapper.getCommands())
    return result.difference(reserved_commands.keys())


reserved_commands: dict[str, str] = {
    "ion-species": "pseudopotentials",
    "elec-ex-corr": "xc",
    "elec-smearing": "smearing",
    "kpoint-folding": "kpts",
    "kpoint": "kpts",
    "elec-initial-charge": "charge",
    "elec-n-bands": "nbands",
}  #: unallowed jdftx input-file commands and the ASE parameters that replace them

pseudopotential_sets: dict[str, str] = {
    "SG15": "SG15/$ID_ONCV_PBE.upf",
    "GBRV": "GBRV/$ID_pbe.uspp",
    "GBRV-pbe": "GBRV/$ID_pbe.uspp",
    "GBRV-lda": "GBRV/$ID_lda.uspp",
    "GBRV-pbesol": "GBRV/$ID_pbesol.uspp",
}  #: mapping of recognized pseudopotential sets to filenames

def default_cutoff(pseudopotentials: str) -> str:
    """Default plane-wave cutoff(s) for recognized pseudopotential sets."""
    if pseudopotentials.startswith("GBRV"):
        return "20 100"
    elif pseudopotentials.startswith("SG15"):
        return "30"
    else:
        return "20"  # default in JDFTx

xc_map: dict[str, str]  = {
    "LDA": "lda",
    "PBE": "gga-PBE",
}  #: maps from ASE standard functional names to JDFTx names


class JDFTx(Calculator):
    """Calculator using JDFTx (https://jdftx.org) through libjdftx."""
    
    implemented_properties = ["energy", "forces", "stress"]
    
    default_parameters = dict(
        xc="PBE",
        kpts=3.5,
        smearing=None,
        charge=0,
        nbands=0,
        pseudopotentials="GBRV",
        commands=None,
    )

    def __init__(self, **kwargs) -> None:
        """
        In addition to standard ASE Calculator options,
        the following keyword arguments are recognized:
        
        xc
            Exchange-correlation functional, using either the ASE standard
            options LDA and PBE, or any combination supported by JDFTx.
        kpts
            Specify as a tuple of 3 integers for a Monkhorst-Pack grid,
            with an additional key 'gamma' to force Gamma-centered mesh.
            Alternately, specify a single k-point density in Angstroms.
        smearing
            Specify as None or (type, width_in_eV), where type may be Fermi-Dirac,
            Gaussian, Methfessel-Paxton or Cold. Note that the width refers
            to the kT for Fermi-Dirac and the Gaussian sigma for all others.
            (This differs from JDFTx input, where sigma = 2 kT for the others.)
            Methfessel-Paxton has an extra argument for order, which must be 1.
        charge
            Total (initial) charge in system. Note that positive implies removing
            electrons, opposite to the JDFTx electron-is-positive convention.
        nbands
            Number of electron bands, with default zero requesting minimum number
            of bands to accommodate the electrons (including margin for smearing
            if present). Can be of the form '0 1.1' to request 1.1x the default
            number of bands, exactly as in elec-n-bands.
        pseudopotentials
            Either the name of a recognized pseudopotential set, including
            SG15, GBRV, GBRV-pbe, GBRV-pbesol, GBRV-lda, or a single string
            with the $ID wildcard to be replaced by element symbols.
            (Individual filenames per element are deliberately not supported
            because these may cause issues with atom ordering. If not using
            these sets, create a directory with the relevant pseudopotentials
            or links to them, so that you can use the wildcard syntax.)
            Recommended cutoffs are automatically set for the recognized
            sets; specify elec-cutoff explicitly in `commands` for others.
        commands
            Recognized JDFTx commands and corresponding arguments specified
            as a dictionary or a list of pairs of command names and arguments
            (all arguments together as a single string). Note that the settings
            listed above that are supported using keyword arguments directly
            to Calculator (following the ASE standard interface) must be
            specified there, and cannot be set within commands.
            See `valid_commands` and `reserved_commands` for more info.
        """
        # Note that all parameters are processed by Calculator
        # and merged with default_parameters specified above.
        # The parameters are all handled then in initialize() below.
        super().__init__(**kwargs)

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        if {"numbers", "pbc"}.intersection(system_changes):
            self.initialize(self.atoms)
        elif {"cell", "positions"}.intersection(system_changes):
            R = self.atoms.cell[:].T  # current lattice vectors
            R_prev = self.atoms_calculated.cell[:].T  # previous lattice vectors
            fractional = self.atoms.get_positions() @ np.linalg.inv(R.T)
            fractional_prev = self.atoms_calculated.get_positions() @ np.linalg.inv(R_prev.T)
            self.jdftx_wrapper.move(fractional - fractional_prev, (R - R_prev) / Bohr)
            self.atoms_calculated = self.atoms.copy()
        
        self.results["energy"] = self.jdftx_wrapper.getEnergy() * Hartree
        self.results["forces"] = self.jdftx_wrapper.getForces() * (Hartree/Bohr)
        self.results["stress"] = self.jdftx_wrapper.getStress() * (Hartree/Bohr**3)

    def initialize(self, atoms: Atoms) -> None:
        # Geometry:
        RT = atoms.get_cell() / Bohr
        lattice_arg = ""
        for i in range(3):
            for j in range(3):
                lattice_arg += f"{RT[j, i]} "
        commands = [("lattice", lattice_arg)]
        commands.append(("coords-type", "Cartesian"))
        positions = atoms.get_positions() / Bohr  # Cartesian positions in bohrs
        names = atoms.get_chemical_symbols()
        for name, (x, y, z) in zip(names, positions):
            commands.append(("ion", f"{name} {x:.15f} {y:.15f} {z:.15f} 1"))

        # Exchange-correlation:
        xc = self.parameters.xc
        if xc in xc_map:
            xc = xc_map[xc]
        commands.append(("elec-ex-corr", xc))

        # K-point mesh:
        kpts = self.parameters.kpts
        force_gamma = False
        if isinstance(kpts, (float, int)):
            folding = kptdensity2monkhorstpack(atoms, kpts)
        elif isinstance(kpts, tuple):
            folding = kpts[:3]
            if len(kpts) == 4 and kpts[3] == "gamma":
                force_gamma = True
        else:
            raise ValueError("kpts must be a density, grid size (+ optional 'gamma')")
        offsets = [(0 if (force_gamma or (fold % 2)) else 0.5) for fold in folding]
        commands.append(("kpoint", f"{offsets[0]} {offsets[1]} {offsets[2]} 1"))
        commands.append(("kpoint-folding", f"{folding[0]} {folding[1]} {folding[2]}"))

        # Smearing:
        smearing = self.parameters.smearing
        if smearing is not None:
            smear_key = smearing[0].lower()
            smear_map = {
                "fermi-dirac": "Fermi",
                "gaussian": "Gaussian",
                "methfessel-paxton": "MP1",
                "cold": "Cold",
            }
            if smear_key not in smear_map:
                raise KeyError(
                    f"Unrecognized smearing {smear_key}; "
                    f"must be one of {', '.join(smear_map.keys())}"
                )
            smear_type = smear_map[smear_key]
            if smear_type == "MP1":
                assert smearing[2] == 1  # only first-order M-P smearing supported
            kT = smearing[1] * (1.0 if smear_type == "Fermi" else 0.5) * eV/Hartree
            commands.append(("elec-smearing", f"{smear_type} {kT}"))

        # Settings that translate (nearly) verbatim:
        commands.append(("elec-initial-charge", f"{-self.parameters.charge}"))
        commands.append(("elec-n-bands", f"{self.parameters.nbands}"))

        # Pseudopotentials:
        pseudopotentials = self.parameters.pseudopotentials
        if pseudopotentials in pseudopotential_sets:
            pseudopotentials = pseudopotential_sets[pseudopotentials]
        commands.append(("ion-species", pseudopotentials))

        # Validate extra commands specified in JDFTx form:
        extra_commands = self.parameters.commands
        if extra_commands is None:
            extra_commands = []
        if isinstance(extra_commands, dict):
            extra_commands = extra_commands.items()
        if not isinstance(commands, list):
            raise ValueError("commands should be a dict or list of command-argument pairs")
        commands_in = set(cmd.lower() for cmd, _ in extra_commands)
        invalid_commands = commands_in.difference(valid_commands())
        if invalid_commands:
            message = f"Invalid command(s): {', '.join(invalid_commands)}."
            for cmd in invalid_commands.intersection(reserved_commands.keys()):
                message += f" Use parameter {reserved_commands[cmd]} instead of {cmd}."
            raise KeyError(message)
        commands.extend(extra_commands)

        # Pseudopotential-set-dependent cutoff:
        if "elec-cutoff" not in commands_in:
            commands.append(("elec-cutoff", default_cutoff(pseudopotentials)))

        self.jdftx_wrapper = JDFTxWrapper(commands, True)
        self.atoms_calculated = atoms.copy()  # atoms for which results are current
    
    @staticmethod
    def help(command: str) -> None:
        """Retrieve documentation of JDFTx commands."""
        reserved = (command in reserved_commands)
        valid = (command in valid_commands())
        if reserved or valid:
            print(JDFTxWrapper.getCommandDoc(command))
        if reserved:
            parameter = reserved_commands[command]
            print(
                f"Note: {command} is valid in JDFTx input files, but it isn't allowed\n"
                f"in the JDFTx calculator. Use parameter {parameter} instead.\n"
            )
