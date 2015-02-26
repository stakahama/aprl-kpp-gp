kppaermod
===

Program to modify outputs of KPP to run gas+aerosol simulations. gas-only simulation program also modified to receive same inputs as gas+aerosol model.


## Directory structure (user inputs)

Compound-specific information and initial conditions (e.g., "compounds/"):

* apinene_1/
  * {ROOT}.def
  * {ROOT}.kpp
  * mcm\_{ROOT}\_mass.txt (downloaded as mcm\_subset\_mass.txt from MCM web)
* apinene_2/
  * {ROOT}.def
  * {ROOT}.kpp
  * mcm\_{ROOT}\_mass.txt (downloaded as mcm\_subset\_mass.txt from MCM web)


Simulations (e.g., "simulations/") (here, the name of subfolders should be "run\_" followed by a three-digit integer):

* apinene_1/
	* run\_001/
	  * input.txt
	  * a0.txt
	* run\_002/
	  * input.txt
	  * a0.txt
* apinene_2/
	* run\_001/
	  * input.txt
	  * a0.txt
	* run\_002/
	  * input.txt
	  * a0.txt

Using the same name for the top level directory (e.g., "apinene\_1") for the compounds and simulation directories may be helpful.

Currently, "apinene\_1", "apinene\_2", etc. should contain runs which vary according to initial gas-phase concentrations and method for deriving photolysis rate constants (technically, the latter can be changed by replacing "photolysis.txt" but currently it is not grouped with input files). After generating the gas and total (gas+aerosol) simulation models ({ROOT}.exe), these parameters can be changed for various simulations ("run\_{DDD}/"):

* The duration and timesteps can be varied in "input.txt".
* The initial aerosol concentration value ("C<sub>OA</sub>") and integrator (*always "dlsode"*) is also specified in "input.txt".
* The initial aerosol mole fractions can be varied according to "a0.txt".


## Instructions

Add kppaermod to the list of paths in which executable are searched:

```
$ export PATH=~/git/projects/kppaermod:$PATH
```

### Calculate group abundances for each compound


*Run in compounds directory (e.g., "compounds/apinene\_1")*.


This uses [aprl-structsearch](https://bitbucket.org/stakahama/aprl-structsearch). Note that according to the instructions for aprl-structsearch, you should add the path of this program to the `PATH` environmental variable also:
```
$ export PATH=~/git/projects/aprl-structsearch:$PATH
```

I have created a script, search_struct.py, in kppaermod to facilitate generation of SIMPOL and FTIR groups, and also the vapor pressures at 298.15K and 358.15K (60 degrees C) using the aprl-structsearch program.


Command:
```
$ search_struct.py {ROOT} {PROGPATH}
```

Arguments:

* `ROOT`: 
* (optional) `PROGPATH`: path to aprl-structsearch if not on executable path


Example usage:
```
$ search_struct.py apinene
```

Outputs (in working directory):

* {ROOT}\_FTIRGroups.csv: table of abundances, compounds x groups
* {ROOT}\_SIMPOLGroups.csv: table of abundances, compounds x groups
* {ROOT}\_props\_298.csv: vapor pressures and enthalpies of vaporization at 298.15K
* {ROOT}\_props\_358.csv: vapor pressures and enthalpies of vaporization at 358.15K


### Build executables for gas-phase only ("gas") and gas+aerosol ("total") simulations

*Run in simulation directory (e.g., "simulations/apinene\_1/")*.

Command:
```
$ build_dual.py {ROOT} {LIGHT} {CPATH} {--skipbuild}
```
Arguments:

* `ROOT`:
* `LIGHT`: one of "original", "light", or "dark" (without quotes)
* `CPATH`: path to compounds directory
* (optional) `--skipbuild`: will not run kpp again but only use output of kpp (in "kppbuild/") to generate "exec\_gas/" and "exec\_total/"

Example usage:
```
$ build_dual.py apinene dark ../compounds/apinene_1
```

Outputs (in working directory):

* contents of {CPATH} are copied here
* kppbuild/: directory of kpp output
* exec\_gas/: executable for gas-phase simultions
* exec\_total/: executable for total-phase simulations

### Execute gas-phase and gas+aerosol simulations

*Run in simulation directory (e.g., "simulations/apinene\_1/").*

Command:
```
$ exec_dual.py {ROOT} {NUMBERS}
```

Arguments:

* `ROOT`: 
* `NUMBERS`: a single number, or several numbers concatenated by commas. e.g., "1", or "1,2" (without quotes) and so on.


Example usage:
```
$ exec_dual.py apinene 1,2
```

Outputs (in run directories):

* run\_{DDD}/gas/{ROOT}.dat
* run\_{DDD}/gas/{ROOT}\_formatted.txt
* run\_{DDD}/total/{ROOT}.dat
* run\_{DDD}/total/{ROOT}\_formatted.txt
* run\_{DDD}/total/{ROOT}\_aer.dat
* run\_{DDD}/total/{ROOT}\_aer\_formatted.txt
