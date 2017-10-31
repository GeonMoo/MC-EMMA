# Copyright (C) 2008 CSC - Scientific Computing Ltd.
"""This module defines an ASE interface to GULP.

The user should also set the environmental flag $GULP_SCRIPT pointing
to a python script looking something like::

   import os
   exitcode = os.system('gulp')

Alternatively, user can set the environmental flag $GULP_COMMAND pointing
to the command use the launch gulp e.g. 'gulp' or 'mpirun -n 16 gulp'

http://projects.ivec.org/gulp/

-Matthew Dyer msd30@liverpool.ac.uk
"""

import os
import sys
import re
from general import Calculator
from os.path import join, isfile, islink

import numpy as np

import ase

class Gulp(Calculator):
    def __init__(self, restart=None, output_template='gulp', track_output=False,
                 keywords=['sing',],options=None,filename='ase-gulp',
                 shells=None,setups=None):
        """Initialize calculator

           keywords: Keywords with no parameters which appear on first line of 
                     GULP input file (Default is sing)
           options: Options with parameters. This should be a list of options 
                    with the attached parameters
           filename: GULP will use the files filename.gin and filename.gout
                     (Default is ase-gulp.gin and ase-gulp.out)
           shells: A dictionary object with atom labels as keys and charges as 
                   values for each atomic species to have a shell attached
           setups: A list to define special labels for specific atoms
                   It must contain pairs of lists of indices and a string to 
                   be appended to the atomic symbol for that atom
                   (e.g. [[[0,1,2],"_A"],[[3,4,5],"B"]])
        """
        self.name = 'Gulp'
        self.keywords = keywords
        self.options = options
        self.filename = filename
        self.shells = shells
        self.setups = setups

        self.restart = restart
        self.track_output = track_output
        self.output_template = output_template
        if restart:
            self.restart_load()
            return

        self.atoms = None
        self.positions = None
        self.run_counts = 0

    def update(self, atoms):
        if self.calculation_required(atoms, ['energy']):
            if (self.atoms is None or
                self.atoms.positions.shape != atoms.positions.shape):
                # Completely new calculation just reusing the same
                # calculator, so delete any old GULP files found.
                self.clean()
            self.calculate(atoms)

    def initialize(self, atoms):
        """Initialize a GULP calculation
        """

        self.all_symbols = atoms.get_chemical_symbols()
        self.natoms = len(atoms)
        atomtypes = atoms.get_chemical_symbols()

        # Check through setups for special labels
        atoms_labels = self.all_symbols
        # Each entry of setups must contain a list of indices followed by a 
        # string
        if self.setups:
            for pair in self.setups:
                for index in pair[0]:
                    atom_labels[index]=atom_labels[index]+pair[1]
        self.atoms_labels = atoms_labels
        self.converged = None
        self.setups_changed = None

    def calculate(self, atoms):
        """Generate necessary files in the working directory and run GULP.

        The method first write GULP input files, then calls the method
        which executes GULP. When the GULP run is finished energy, forces,
        etc. are read from the GULP output.
        """

        # Write input
        from ase.io.gulp import write_gulp
        self.initialize(atoms)
        write_gulp(self.filename+".gin", atoms, 
                   keywords=self.keywords, options=self.options, 
                   shells=self.shells, symbols = self.atoms_labels)

        # Execute GULP
        self.run()
        # Read output if relaxation has happened
        if (self.keywords.count('opti')>0 or self.keywords.count('optimise')>0 or 
           self.keywords.count('grad')>0 or self.keywords.count('gradient')>0 or 
           self.keywords.count('fit')>0):
	   from ase.io.gulp import read_gulp_out as rgo
           tmp_atoms=rgo(self.filename+'.gout')
           atoms.positions = tmp_atoms.positions
           atoms.cell = tmp_atoms.cell
        self.converged = self.read_convergence()
        self.set_results(atoms)

    def set_results(self, atoms):
        self.read(atoms)
        self.old_keywords = self.keywords
        self.old_options = self.options
        self.old_shells = self.shells
        self.old_setups = self.setups
        self.old_filename = self.filename
        self.atoms = atoms.copy()
        self.version = self.read_version()
        self.niter = self.read_number_of_iterations()

    def run(self):
        """Method which explicitely runs GULP."""

        if self.track_output:
            self.out = self.output_template+str(self.run_counts)+'.out'
            self.run_counts += 1
        else:
            self.out = self.output_template+'.out'
        if os.environ.has_key('GULP_COMMAND'):
            gulp = os.environ['GULP_COMMAND']
            exitcode = os.system('%s %s > %s' % (gulp, self.filename, self.out))
        elif os.environ.has_key('GULP_SCRIPT'):
            gulp = os.environ['GULP_SCRIPT']
            locals={}
            execfile(vasp, {}, locals)
            exitcode = locals['exitcode']
        else:
            raise RuntimeError('Please set either GULP_COMMAND or GULP_SCRIPT environment variable')
        if exitcode != 0:
            raise RuntimeError('Gulp exited with exit code: %d.  ' % exitcode)

    def clean(self):
        """Method which cleans up after a calculation.

        The default files generated by GULP will be deleted IF this
        method is called.

        """
        files = [self.filename+".gin",self.filename+".gout"]
        for f in files:
            try:
                os.remove(f)
            except OSError:
                pass

    def calculation_required(self, atoms, quantities):
        if (self.positions is None or
            (self.atoms != atoms) or
            (self.keywords != self.old_keywords) or
            (self.options != self.old_options) or
            (self.shells != self.old_shells) or
            (self.setups != self.old_setups) or
            (self.filename != self.old_filename) or
            not self.converged):
            return True
        return False

    def set_atoms(self, atoms):
        if (atoms != self.atoms):
            self.converged = None
        self.atoms = atoms.copy()

    def get_atoms(self):
        atoms = self.atoms.copy()
        atoms.set_calculator(self)
        return atoms

    def get_version(self):
        self.update(self.atoms)
        return self.version

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        return self.energy

    def get_number_of_iterations(self):
        self.update(self.atoms)
        return self.niter

    def read_number_of_iterations(self):
        niter = None
        for line in open(self.filename+'.gout','r'):
            if line.find('Cycle:') != -1: # find the last iteration number
                niter = int(line.split()[1])
        return niter

    def get_forces(self, atoms):
        self.update(atoms)
        return self.forces

    def get_stress(self, atoms):
        self.update(atoms)
        return self.stress

    def get_dipole_moment(self, atoms):
        """Returns total dipole moment of the system."""
        self.update(atoms)
        return self.dipole

  # Methods for reading information from output files:
    def read(self, atoms):
       self.positions = atoms.get_positions()
       self.energy = self.read_energy()
       self.forces = self.read_forces(atoms)
       self.stress = self.read_stress(atoms)
       self.atoms = atoms.copy()
       try:
           self.nbands = self.read_nbands()
       except NotImplementedError:
           do_nothing = True
       except AttributeError:
           do_nothing = True
       try:
           self.dipole = self.read_dipole()
       except NotImplementedError:
           do_nothing = True
       return

    def read_energy(self):
        energy=0
        for line in open(self.filename+'.gout', 'r'):
            # Free energy
            if line.find('Primitive unit cell')!=-1 and line.find('eV')!=-1:
		try:
			energy = float(line.split()[-2]) 
		except ValueError: 
			energy =line.split()[-2]
            elif line.find('Total lattice energy')!=-1 and line.find('eV')!=-1:
		try:
			energy = float(line.split()[-2]) 
		except ValueError: 
			energy =line.split()[-2]
        return energy

    def read_forces(self, atoms):
        """Method that reads forces from output file.
        """
        file = open(self.filename+'.gout','r')
        data = file.readlines()
        file.close()
        forces=None
        for n,line in enumerate(data):
            if 'Total number atoms/shells' in line:
                ncores_shells=int(data[n].split()[4])
            if 'Final internal derivatives :' in line:
                forces = []
                for iatom in range(ncores_shells):
                    temp    = data[n+6+iatom].split()
                    #Only include cores
                    if temp[2]=='c':
			try:
	                        forces += [[float(temp[3]),float(temp[4]),float(temp[5])]]
			except(TypeError,ValueError):
				forces = None
        return np.array(forces)

    def read_stress(self,atoms):
        stress=[]
        for line in open(self.filename+'.gout','r'):
            if line.find('dE/de1(xx)') != -1:	
		try:
	                stress.append(float(line.split()[4]))
		except (TypeError,ValueError):
			pass
            if line.find('dE/de2(yy)') != -1:
		try:
	                stress.append(float(line.split()[4]))
		except (TypeError,ValueError):
			pass
            if line.find('dE/de3(zz)') != -1:
		try:
	               stress.append(float(line.split()[4]))
		except(TypeError,ValueError):
			pass
            if line.find('dE/de4(yz)') != -1:
		try:
	                stress.append(float(line.split()[4]))
		except(TypeError,ValueError):	
                        pass
            if line.find('dE/de5(xz)') != -1:
		try:
	                stress.append(float(line.split()[4]))
		except(TypeError,ValueError):
			pass
            if line.find('dE/de6(xy)') != -1:
		try:
	                stress.append(float(line.split()[4]))
		except:
			pass
        #Convert from eV/strain to eV/Angstrom^3
        stress = -np.array(stress)/atoms.get_volume()
        return stress

    def read_dipole(self):
        #dipolemoment=np.zeros([1,3])
        #for line in open('OUTCAR', 'r'):
        #    if line.rfind('dipolmoment') > -1:
        #        dipolemoment=np.array([float(f) for f in line.split()[1:4]])
        #return dipolemoment
        raise NotImplementedError

    def read_convergence(self):
        """Method that checks whether a calculation has converged."""
        converged = None
        # First check electronic convergence
        for line in open(self.filename+'.gout', 'r'):
            if line.rfind('Components of energy') > -1:
                converged = True
        # Then if ibrion > 0, check whether ionic relaxation condition been fulfilled
        if (self.keywords.count('opti')>0 or self.keywords.count('optimise')>0 or 
           self.keywords.count('grad')>0 or self.keywords.count('gradient')>0 or 
           self.keywords.count('fit')>0):
            if not self.read_relaxed():
                converged = False
            else:
                converged = True
        return converged

    def read_relaxed(self):
        for line in open(self.filename+'.gout', 'r'):
            if line.rfind('**** Optimisation achieved ****') > -1:
                return True
        return False

    def read_version(self):
        version = None
        for line in open(self.filename+'.gout'):
            if line.find('* Version =') != -1: # find the first occurence
                version = line[len('* Version ='):].split()[0]
                break
        return version
