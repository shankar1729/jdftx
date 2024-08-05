#!/usr/bin/env python

# Atomistic Simulation Environment (ASE) calculator interface for JDFTx
# See http://jdftx.org for JDFTx and https://wiki.fysik.dtu.dk/ase/ for ASE
# Authors: Deniz Gunceler, Ravishankar Sundararaman

from __future__ import print_function #For Python2 compatibility

import os, scipy, subprocess, tempfile, re
from ase.calculators.calculator import Calculator
from ase.units import Bohr, Hartree

#Run shell command and return output as a string:
def shell(cmd):
	return subprocess.check_output(cmd, shell=True)

#Return var, replacing it with environment variable varName if it is None
def replaceVariable(var, varName):
	if (var is None) and (varName in os.environ):
		return os.environ[varName]
	else:
		return var

class JDFTx(Calculator):

	def __init__(self, executable=None, pseudoDir=None, pseudoSet='GBRV', commands=None, ignoreStress=False):
		#Valid pseudopotential sets (mapping to path and suffix):
		pseudoSetMap = {
			'SG15' : 'SG15/$ID_ONCV_PBE.upf',
			'GBRV' : 'GBRV/$ID_pbe.uspp',
			'GBRV-pbe' : 'GBRV/$ID_pbe.uspp',
			'GBRV-lda' : 'GBRV/$ID_lda.uspp',
			'GBRV-pbesol' : 'GBRV/$ID_pbesol.uspp'
		}

		#Get default values from environment:
		self.executable = replaceVariable(executable, 'JDFTx')      #Path to the jdftx executable (cpu or gpu)
		self.pseudoDir = replaceVariable(pseudoDir, 'JDFTx_pseudo') #Path to the pseudopotentials folder
		self.ignoreStress = ignoreStress
		
		if (self.executable is None):
			raise Exception('Specify path to jdftx in argument \'executable\' or in environment variable \'JDFTx\'.')
		if (self.pseudoDir is None) and (not (pseudoSet in pseudoSetMap)):
			raise Exception('Specify pseudopotential path in argument \'pseudoDir\' or in environment variable \'JDFTx_pseudo\', or specify a valid \'pseudoSet\'.')
		
		if pseudoSet in pseudoSetMap:
			self.pseudoSetCmd = 'ion-species ' + pseudoSetMap[pseudoSet]
		else:
			self.pseudoSetCmd = ''
		
		# Gets the input file template
		self.acceptableCommands = set(['electronic-SCF'])
		template = str(shell('%s -t' % (self.executable)))
		for match in re.findall(r"# (\S+) ", template):
			self.acceptableCommands.add(match)

		self.input = [
			('dump-name', '$VAR'),
			('initial-state', '$VAR')
		]

		# Parse commands, which can be a dict or a list of 2-tuples.
		if isinstance(commands, dict):
			commands = commands.items()
		elif commands is None:
			commands = []
		for cmd, v in commands:
			self.addCommand(cmd, v)

		# Accepted pseudopotential formats
		self.pseudopotentials = ['fhi', 'uspp', 'upf']

		# Current results
		self.E = None
		self.Forces = None

		# History
		self.lastAtoms = None
		self.lastInput = None

		# k-points
		self.kpoints = None

		# Dumps
		self.dumps = []
		self.addDump("End", "State")
		self.addDump("End", "Forces")
		self.addDump("End", "Ecomponents")
		
		#Run directory:
		self.runDir = tempfile.mkdtemp()
		print('Set up JDFTx calculator with run files in \'' + self.runDir + '\'')

	########### Interface Functions ###########

	def addCommand(self, cmd, v):
		if(not self.validCommand(cmd)):
			raise IOError('%s is not a valid JDFTx command!\nLook at the input file template (jdftx -t) for a list of commands.' % (cmd))
		self.input.append((cmd, v))

	def addDump(self, when, what):
		self.dumps.append((when, what))

	def addKPoint(self, pt, w=1):
		b1, b2, b3 = pt
		if self.kpoints is None:
			self.kpoints = []
		self.kpoints.append((b1, b2, b3, w))
	
	def clean(self):
		shell('rm -rf ' + self.runDir)

	def calculation_required(self, atoms, quantities):
		if((self.E is None) or (self.Forces is None)):
			return True
		if((self.lastAtoms != atoms) or (self.input != self.lastInput)):
			return True
		return False

	def get_forces(self, atoms):
		if(self.calculation_required(atoms, None)):
			self.update(atoms)
		return self.Forces

	def get_potential_energy(self, atoms, force_consistent=False):
		if(self.calculation_required(atoms, None)):
			self.update(atoms)
		return self.E

	def get_stress(self, atoms):
		if self.ignoreStress:
			return scipy.zeros((3,3))
		else:
			raise NotImplementedError('Stress calculation not implemented in JDFTx interface: set ignoreStress=True to ignore.')

	################### I/O ###################

	def __readEnergy(self, filename):
		Efinal = None
		for line in open(filename):
			tokens = line.split()
			if len(tokens)==3:
				Efinal = float(tokens[2])
		if Efinal is None:
			raise IOError('Error: Energy not found.')
		return Efinal * Hartree #Return energy from final line (Etot, F or G)

	def __readForces(self, filename):
		idxMap = {}
		symbolList = self.lastAtoms.get_chemical_symbols()
		for i, symbol in enumerate(symbolList):
			if symbol not in idxMap:
				idxMap[symbol] = []
			idxMap[symbol].append(i)
		forces = [0]*len(symbolList)
		for line in open(filename):
			if line.startswith('force '):
				tokens = line.split()
				idx = idxMap[tokens[1]].pop(0) # tokens[1] is chemical symbol
				forces[idx] = [float(word) for word in tokens[2:5]] # tokens[2:5]: force components

		if(len(forces) == 0):
			raise IOError('Error: Forces not found.')
		return (Hartree / Bohr) * scipy.array(forces)

	############## Running JDFTx ##############

	def update(self, atoms):
		self.runJDFTx(self.constructInput(atoms))

	def runJDFTx(self, inputfile):
		""" Runs a JDFTx calculation """
		#Write input file:
		fp = open(self.runDir+'/in', 'w')
		fp.write(inputfile)
		fp.close()
		#Run jdftx:
		shell('cd %s && %s -i in -o out' % (self.runDir, self.executable))
		self.E = self.__readEnergy('%s/Ecomponents' % (self.runDir))
		self.Forces = self.__readForces('%s/force' % (self.runDir))

	def constructInput(self, atoms):
		""" Constructs a JDFTx input string using the input atoms and the input file arguments (kwargs) in self.input """
		inputfile = ''

		# Add lattice info
		R = atoms.get_cell() / Bohr
		inputfile += 'lattice \\\n'
		for i in range(3):
			for j in range(3):
				inputfile += '%f  ' % (R[j, i])
			if(i != 2):
				inputfile += '\\'
			inputfile += '\n'

		# Construct most of the input file
		inputfile += '\n'
		for cmd, v in self.input:
			inputfile += '%s %s\n' % (cmd, str(v))

		# Add ion info
		atomPos = [x / Bohr for x in list(atoms.get_positions())]  # Also convert to bohr
		atomNames = atoms.get_chemical_symbols()   # Get element names in a list
		inputfile += '\ncoords-type cartesian\n'
		for i in range(len(atomPos)):
			inputfile += 'ion %s %f %f %f \t 1\n' % (atomNames[i], atomPos[i][0], atomPos[i][1], atomPos[i][2])
		del i

		# Add k-points
		if self.kpoints:
			for pt in self.kpoints:
				inputfile += 'kpoint %.8f %.8f %.8f %.14f\n' % pt

		#Add pseudopotentials
		inputfile += '\n'
		if not (self.pseudoDir is None):
			added = []  # List of pseudopotential that have already been added
			for atom in atomNames:
				if(sum([x == atom for x in added]) == 0.):  # Add ion-species command if not already added
					for filetype in self.pseudopotentials:
						try:
							shell('ls %s | grep %s.%s' % (self.pseudoDir, atom, filetype))
							inputfile += 'ion-species %s/%s.%s\n' % (self.pseudoDir, atom, filetype)
							added.append(atom)
							break
						except:
							pass
		inputfile += self.pseudoSetCmd + '\n' #Pseudopotential sets

		# Add truncation info (periodic vs isolated)
		inputfile += '\ncoulomb-interaction '
		pbc = list(atoms.get_pbc())
		if(sum(pbc) == 3):
			inputfile += 'periodic\n'
		elif(sum(pbc) == 0):
			inputfile += 'isolated\n'
		elif(sum(pbc) == 1):
			inputfile += 'wire %i%i%i\n' % (pbc[0], pbc[1], pbc[2])
		elif(sum(pbc) == 2):
			inputfile += 'slab %i%i%i\n' % (not pbc[0], not pbc[1], not pbc[2])
		#--- add truncation center:
		if(sum(pbc) < 3):
			center = scipy.mean(scipy.array(atomPos), axis=0)
			inputfile += 'coulomb-truncation-embed %g %g %g\n' % tuple(center.tolist())
		
		#Add dump commands
		inputfile += "".join(["dump %s %s\n" % (when, what) for when, what in self.dumps])

		# Cache this calculation to history
		self.lastAtoms = atoms.copy()
		self.lastInput = list(self.input)
		return inputfile

	############## JDFTx command structure ##############

	def validCommand(self, command):
		""" Checks whether the input string is a valid jdftx command \nby comparing to the input template (jdft -t)"""
		if(type(command) != str):
			raise IOError('Please enter a string as the name of the command!\n')
		return command in self.acceptableCommands

	def help(self, command=None):
		""" Use this function to get help about a JDFTx command """
		if(command is None):
			print('This is the help function for JDFTx-ASE interface. \
				  \nPlease use the command variable to get information on a specific command. \
				  \nVisit jdftx.sourceforge.net for more info.')
		elif(self.validCommand(command)):
			raise NotImplementedError('Template help is not yet implemented')
		else:
			raise IOError('%s is not a valid command' % (command))
