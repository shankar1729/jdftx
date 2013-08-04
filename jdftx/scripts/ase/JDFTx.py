#!/usr/bin/env python

# Atomistic Simulation Environment (ASE) calculator interface for JDFTx
# JDFTx is located at http://jdftx.sourceforge.net
# ASE is located at https://wiki.fysik.dtu.dk/ase/
#
# Author: Deniz Gunceler, June 2013
#
# If you use this interface (or JDFTx) without proper credit,
# gnomes will come to haunt you!

from os import system as shell            # Put system calls
from commands import getoutput as shellO  # Put system calls and retrieve output
import copy, random, re, scipy

from ase.calculators.interface import Calculator
from ase.units import Bohr, Hartree


class JDFTx(Calculator):

    def __init__(self, executable='$JDFTx', pseudoDir='$JDFTx_pseudo', commands={}):

        self.executable = shellO('echo %s' % (executable))  # Path to the jdftx executable (cpu or gpu)
        self.pseudoDir = shellO('echo %s' % (pseudoDir))    # Path to the pseudopotentials folder

        # Gets the input file template
        self.template = shellO('%s -t' % (executable))

        # Check commands for consistency
        for item in commands.items():
            if(not self.validCommand(item[0])):
                raise IOError('%s is not a valid JDFTx command!\nLook at the input file template (jdftx -t) for a list of commands.' % (item[0]))
        self.input = copy.deepcopy(commands)

        # Accepted pseudopotential formats
        self.pseudopotentials = ['fhi', 'uspp']

        # Current results
        self.E = None
        self.Forces = None

        # History
        self.lastAtoms = None
        self.lastInput = None

    ########### Interface Functions ###########

    def calculation_required(self, atoms, quantities):

        if((self.E == None) or (self.Forces == None)):
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

        raise NotImplementedError('Stress interface with JDFTx are not yet implemented')

    ################### I/O ###################

    def __readEnergy(self, filename):

        readin = shellO('cat %s | grep Etot' % (filename)).split()
        if(len(readin) == 0):  # If could not find Etot, search for F (used in fermi surface smearing of metals)
            readin = shellO('cat %s | grep F' % (filename)).split()
        if(len(readin) == 0):  # If still couldn't find total energy, there's a problem with the file
            raise IOError('Could not find Etot or F in %s!' % (filename))

        return float(readin[2]) * Hartree

    def __readForces(self, filename):

        readin = shellO('cat %s' % (filename)).splitlines()

        forces = []
        for line in readin:
            if(line.find('force ') > -1):
                word = line.split()
                forces.append(scipy.array([float(word[2]), float(word[3]), float(word[4])]))
        del line

        if(len(forces) == 0):
            raise IOError('Error: Forces not found.')

        return (Hartree / Bohr) * scipy.array(forces)

    ############## Running JDFTx ##############

    def update(self, atoms):

        self.runJDFTx(self.constructInput(atoms))

    def runJDFTx(self, inputfile):
        """ Runs a JDFTx calculation """

        # Make a temp directory
        directory = 'temp.%s' % (int(round(100000 * random.random())))
        shell('mkdir %s' % (directory))

        shell('cd %s && echo \'%s\' | %s -o temp.out' % (directory, inputfile, self.executable))

        self.E = self.__readEnergy('%s/temp.Ecomponents' % (directory))
        self.Forces = self.__readForces('%s/temp.force' % (directory))

        # Delete the temp directory
        shell('rm -rf %s' % (directory))

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
        for item in self.input.items():
            inputfile += '%s %s\n' % (item[0], str(item[1]))

        # Add ion info
        atomPos = [x / Bohr for x in list(atoms.get_positions())]  # Also convert to bohr
        atomNames = atoms.get_chemical_symbols()   # Get element names in a list
        inputfile += '\ncoords-type cartesian\n'
        for i in range(len(atomPos)):
            inputfile += 'ion %s %f %f %f \t 1\n' % (atomNames[i], atomPos[i][0], atomPos[i][1], atomPos[i][2])
        del i

        # Add pseudopotentials
        inputfile += '\n'
        added = []  # List of pseudopotential that have already been added
        for atom in atomNames:
            if(sum([x == atom for x in added]) == 0.):  # Add ion-species command if not already added
                for filetype in self.pseudopotentials:
                    searchPseudo = shellO('ls %s | grep %s.%s' % (self.pseudoDir, atom, filetype))
                    if(len(searchPseudo.splitlines()) > 0):
                        inputfile += 'ion-species %s/%s.%s\n' % (self.pseudoDir, atom, filetype)
                        added.append(atom)
                        break

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

        # Add dump commands
        inputfile += '\ndump-name temp.$VAR'
        inputfile += '\ndump End Forces Ecomponents'

        # Cache this calculation to history
        self.lastAtoms = copy.deepcopy(atoms)
        self.lastInput = copy.deepcopy(self.input)

        return inputfile

    ############## JDFTx command structure ##############

    def validCommand(self, command):
        """ Checks whether the input string is a valid jdftx command \nby comparing to the input template (jdft -t)"""

        if(type(command) != str):
            raise IOError('Please enter a string as the name of the command!\n')

        return (self.template.find('# %s ' % (command)) > -1)

    def help(self, command=None):

        """ Use this function to get help about a JDFTx command """

        if(command == None):
            print('This is the help function for JDFTx-ASE interface. \
                  \nPlease use the command variable to get information on a specific command. \
                  \nVisit jdftx.sourceforge.net for more info.')
        elif(self.validCommand(command)):
            raise NotImplementedError, 'Template help is not yet implemented'
        else:
            raise IOError('%s is not a valid command' % (command))
