from typing import Optional

import numpy as np
from ase import Atoms
from ase.units import Bohr, Hartree
from ase.calculators.calculator import Calculator, all_changes

from .. import pyjdftx


class JDFTx(Calculator):
    """Calculator using JDFTx (https://jdftx.org) through libjdftx."""
    implemented_properties = ["energy", "forces", "stress"]

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        if ("numbers" in system_changes) or ("pbc" in system_changes):
            self.initialize(self.atoms)
        
        if ("cell" in system_changes) or ("positions" in system_changes):
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
        # Prepare geometry commands:
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
        commands.extend(self.parameters["commands"])
        
        self.jdftx_wrapper = pyjdftx.JDFTxWrapper(commands, True)
        self.atoms_calculated = atoms.copy()  # atoms for which results are current

        
