from typing import Optional, Union

import numpy as np
from ase import Atoms
from ase.units import Bohr, Hartree
from ase.calculators.calculator import Calculator, all_changes

from .. import JDFTxWrapper


class JDFTx(Calculator):
    """Calculator using JDFTx (https://jdftx.org) through libjdftx."""
    implemented_properties = ["energy", "forces", "stress"]

    def __init__(
        self,
        *,
        commands: Union[None, dict[str, str], list[tuple[str, str]]] = None,
        **kwargs
    ) -> None:
        # Validate extra commands specified:
        if commands is None:
            commands = []
        if isinstance(commands, dict):
            commands = commands.items()
        if not isinstance(commands, list):
            raise ValueError("commands should be a dict or list of command-argument pairs")
        commands_in = set(cmd.lower() for cmd, _ in commands)
        valid_commands = set(cmd.lower() for cmd in JDFTxWrapper.getCommands())
        invalid_commands = commands_in.difference(valid_commands)
        if invalid_commands:
            raise KeyError(f"Invalid JDFTx command(s): {', '.join(invalid_commands)}")
    
        super().__init__(commands=commands, **kwargs)

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
        
        self.jdftx_wrapper = JDFTxWrapper(commands, True)
        self.atoms_calculated = atoms.copy()  # atoms for which results are current

        
