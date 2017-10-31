ase-gulp

These files provide the python package ase-gulp, which integrates the classical
force-field programme Gulp (http://projects.ivec.org/gulp/) with the Atomic 
Simulation Environment (https://wiki.fysik.dtu.dk/ase/).

Following installation you will be able to read and write Gulp files from within
ase and you should be able to run Gulp as one of ase internal calculators.

INSTALLATION
------------

The main difficulty with installation is that the ase/io/__init__.py file found 
in your ase distribution needs to be altered to include references to the gulp
file types.

This can be done by copying your exisiting ase/io/__init__.py file from the ase
installation directory into the ase/io directory within ase-gulp

cp <ase_dir>/ase/io/__init__.py ase/io

Then change directory to the ase-gulp ase/io directory

cd ase/io

And run the patch file provided with the ase-gulp package

patch < ../../__init__.patch

This should run without errors. If there are errors, please contact the authors.

After this, installation should be a simple case of changing directory back into
the ase-gulp root directory and running setup.py as usual.

cd ../../
python setup.py install

Or if you need to be superuser to install new python packages:

sudo python setup.py install

That should be done. Please note that you will need to do this after each update
of ase since the ase/io/__init__.py file needs to be updated each time.

USAGE
-----

Once ase-gulp is installed, structures should be readable using the usual ase 
read() command from files ending in the .res, .gin and .gout suffixes.
Alternatively one can specify type="gulp" as a keywords within the read() 
command. 

Similarly, files can be written using the ase write() command. There are some
additonal optional gulp specific keywords:

write_gulp(filename, images, keywords=['sing',], options=None, shells = None, 
           symbols = None)

keywords: Keywords with no parameters which appear on first line of 
          GULP input file (Default is sing)
options: Options with parameters. This should be a list of options 
         with all their parameters.
         e.g. options = ["library my_library.lib","dump every 1 structure.res"]
shells: A dictionary object with atom labels as keys and charges as 
        values for each atomic species to have a shell attached
        e.g. shells = {"O":-2.8,"Ba":0.2}
symbols: A list of symbols for each atom if they are not wanted to be
         the same as in the original atoms object

The gulp calculator object should work in the same way as other ase
calculators, inparticular it is similar to the VASP calculator.

To import and set up the calculator use

from ase.calculators.gulp import Gulp

calculator=Gulp()

The same keywords are available as for writing files, with the additon of the
filename keyword
filename: GULP will use the files filename.gin and filename.gout
         (Default is ase-gulp.gin and ase-gulp.out)



