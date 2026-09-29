The primary use of this script is as an alternative, Open, back-end to
Materials Studio's 'Calculate Elastic Constants' feature.

It was designed with CASTEP in mind, and only deals with the set of strains
that Materials Studio generates. As well as CASTEP's native formats, it also
accepts input in CML format, and can generate CML output. It also has the
option to generate an at-a-glance view, which allows the user to verify the
success, or otherwise, of all of the various linear fits.

It should be fairly easy to convert/generalize it for other strains and
simulation codes.

How-To Guide
------------

Below is a list of steps to follow to use this set of scripts to obtain the elastic constants for a structure, starting from an optimised geometry and using a set of shears.

**First:** Starting with a `<seedname>.cell` and `<seedname>.param` file for your structure, perform the CASTEP `GeometryOptimisation` task. Note that, to obtain *good* elastic constants, you will need to use tighter tolerances on energies, forces, and atomic displacements than you would in a usual geometry optimsation. See, *e.g.* the [CASTEP webiste](https://www.tcm.phy.cam.ac.uk/castep/documentation/WebHelp/content/modules/castep/tskcastepsetelecquality.htm) for a list of suggested tolerances. Make sure that you have the options `WRITE_CELL_STRUCTURE: TRUE` and `CALCULATE_STRESS: TRUE`, as these will be important later.

**Second:** Move the optimised cell to the file `<seedname>.cell`. (If you like, save the original geometry to something like `<seedname>.cell.old` first.)

**Third:** Run the command
```
python3 /path/to/castep_elastic_constants/generate_strain.py <seedname> <options>
```
to generate the strained structures. Note that options like `numsteps = N` and `maxstrain = s` will control how many strains are generated and how large they will be for each shearing pattern, respectively. This will generate a whole series of `<seedname>_cij__i__j.param` and `<seedname>_cij__i__j.cell` files with the cells and atomic coordinates in the `.cell` files sheared appropriately. Note also that the cell constraints will be updated, with `FIX_ALL_CELL true`.

**Fourth:** Loop over all of the sheared structures and perform fixed cell geometry optimisations on all of these using a script which looks something like this (if using slurm):
```
#!/bin/bash
#SBATCH --nodes=N
#SBATCH --ntasks-per-node=n
#SBATCH --mem-per-cpu=m
#SBATCH --time=24:00:00

module restore CASTEP

cd $SLURM_SUBMIT_DIR

for i in `seq 1 6`
do
  for j in `seq 1 6`
  do
    srun -n <n*N> castep.mpi a1_feni_cij__${i}__${j}
  done
done

exit 0
```
Once this has finished and you have successfully generated all the `<seedname>_cij__i__j.castep` files (check for warnings!) you will then have data to obtain the elastic constants.

**Fifth:** Run the command
```
python3 /path/to/castep_elastic_constants/elastics.py <seedname> <options>
```
to return the elastic constants (and associated uncertainties) for your structure, if desired. You can find information about the relevant options for this script in the source code.


Supported Strain Patterns
-------------------------

Cubic: e1+e4
Hexagonal: e3 and e1+e4
Trigonal-High (32, 3m, -3m): e1 and e3+e4
Trigonal-Low (3, -3): e1 and e3+e4
Tetragonal: e1+e4 and e3+e6
Orthorhombic: e1+e4 and e2+e5 and e3+e6
Monoclinic: e1+e4 and e3+e6 and e2 and e5
Triclinic: e1 to e6 separately 

Dependencies
------------

The script is written in Python, so you must have that installed. I use Python 2.5, but earlier versions might work too. 

SciPy (http://www.scipy.org/)

NumPy (https://numpy.org)

(Optional) For CML input/output, you'll need Golem (http://www.lexical.org.uk/golem/), and lxml (http://codespeak.net/lxml/).

(Optional) For graphical output, you'll need Matplotlib (http://matplotlib.sourceforge.net/).


