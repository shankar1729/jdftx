This is the ASE calculator class for JDFTx

Copy jdftx.py into your ase calculators folder (/path-to-ASE/ase/calculators/)
You can then import it like any other ASE calculator
    from ase.calculators.jdftx import JDFTx

JDFTx interface needs to know where your executable and pseudopotentials are located.
This can be done at runtime by specifying paths to the class constructor, 
or by setting environment variables in your system.

A JDFTx calculator object uses variables executable and pseudoDir 
to find the executable and the pseudopotential directory.
By default, they are set to the envionment variables $JDFTx and $JDFTx_pseudo.
These can be set in the constructor or changed later.

The easiest way to use the interface is to set these environment variables.
On Ubuntu, the easiest way to accomplish this is to have the lines
    export JDFTx=/path/to/jdftx
    export JDFTx_pseudo=/path/to/pseudopotential/directory
in your ~/.bashrc or ~/.profile files.

On clusters with resource management programs (e.g. slurm, torque...),
you can override the executable with a command that submits a job to the queue.
For example, on slurm, the following would work:
    calculator = JDFTx(executable='srun /path/to/jdftx')
Take note that when you do that, the command that submits the job (in this case, srun)
should not release the shell until the job is completed.
For example, in slurm, srun would work, but not sbatch.

/*-------------------------------------------------------------------
Copyright 2012 Deniz Gunceler

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/
