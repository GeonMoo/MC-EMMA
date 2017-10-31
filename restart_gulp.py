import sys
import os
from ase.all import *
from ase.io.gulp import *

#target = "temp01.res"
#restart= "restart.res"
#head = "head2.txt"

def restart_gulp(target='', restart='', head=''):
	converged = False
	t = False
	if os.path.isfile("ase-gulp.gout"):
		for line in open("ase-gulp.gout",'r'):
			if line.rfind("**** Maximum number of function calls has been reached ****") > -1: # check to make sure that the previous part of the run completed
				t=True
	else:
		t = False
	if t == True:		
		if os.path.isfile(target):
			a = open(target,'r') #1. extract structure, inc. shells from previous run
			b = a.readlines()
		else:
			converged = False
			return converged,energy,atoms
		index = []
		for i in range(len(b)):
			if "# Options" in b[i]: # these are the markers in gulp output files above/below the structure
				index.append(i)
			if "totalenergy" in b[i]:
				index.append(i)
				
		new1 = open(head,'r').readlines() #open input options for gulp restart
		new2 = b[index[0]:index[1]] # select just the structure positions
		new3 = open(restart,'w') # make restart file
		
		for i in range(len(new1)):
			new3.write(new1[i]) # write header to restart file
			
		new3.write("\ndump temp.res\n") # add in option for dump file
		
		for i in range(len(new2)):
			new3.write(new2[i]) # write structure to restart file
		new3.close()	
		
		gulp=os.environ['GULP_COMMAND'] # re-run gulp
		inpu=restart
		command=str(gulp + " < " + inpu + " > ase-gulp.gout")
		os.system(command)
		
		converged = False # check for convergance from restart file
		for line in open("ase-gulp.gout",'r'):
			if line.rfind('**** Optimisation achieved ****') > -1:
				converged=True
		
		energy = 0. # read energy from restart file
		for line in open("ase-gulp.gout",'r'):
			if line.rfind('Total lattice energy')!=-1 and line.find('eV')!=-1:
				try:
					energy = float(line.split()[-2])
					atoms = read_gulp("temp.res") # return relaxed atoms object
				except ValueError:
					atoms = read_gulp(target)
					converged = False
					energy = 1.0e20
		
			
		return converged,energy,atoms
				
	else:		# if previous run did not finish it's cyclces, do nothing.
		converged=False
		energy = 0
		atoms = read_gulp(target)
		return converged,energy,atoms
	
