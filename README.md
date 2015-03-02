kppaermod
===

Program to modify outputs of KPP to run gas+aerosol simulations. Gas-only simulation program also modified to receive same inputs as gas+aerosol model.


## User inputs

### Directory structure

The user should provide compound-specific information and initial conditions (e.g., in "compounds/") and simulations parameters (e.g., in "simulations/"). "photolysisfiles/" are provided and the appropriate "photolysis.txt" should be copied into the run subdirectory. Using the same name for the top level directory (e.g., "apinene\_1") for the compounds and simulation directories may be helpful. 

compounds/

* apinene_1/
	* {ROOT}.kpp
	* mcm\_{ROOT}\_mass.txt
* apinene_2/
	* {ROOT}.kpp
	* mcm\_{ROOT}\_mass.txt

simulations/

* apinene\_1/
	* run\_001/
  		* (for gas,total) photolysis.txt
  		* (for gas,total) [optional] input\_time.txt
  		* (for gas,total) [optional] input\_temp.txt
  		* (for gas,total) [optional] cgas\_init.def
  		* (for total) input\_partitioning.txt
  		* (for total) [optional] molefrac\_init.txt
	* run\_002/
  		* (for gas,total) photolysis.txt
  		* (for gas,total) [optional] input\_time.txt
  		* (for gas,total) [optional] input\_temp.txt
  		* (for gas,total) [optional] cgas\_init.def
  		* (for total) input\_partitioning.txt
  		* (for total) [optional] molefrac\_init.txt
* apinene\_2/ (*same structure as above*)

photolysisfiles/

* original/
	* photolysis.txt
* dark/
	* photolysis.txt
* constantlight/
	* photolysis.txt

The main objective is to build a program for a fixed mechanism (set of chemical reactions, species) to simulate over a range of temperatures, concentrations, and timesteps. After generating the gas and total (gas+aerosol) simulation models ("exec\_gas/{ROOT}.exe" or "exec\_total/{ROOT}.exe" in each simulation subdirectory), parameters can be changed through input files for various simulations ("run\_{DDD}/"). Note that runs using "input\_temp.txt" and "cgas\_init.txt" are untested and should be against a reference simulation. 

### File descriptions

Mechanism information:

- {ROOT}.def: combines organic and inorganic kpp files; specifies initial concentrations, temperature, and time parameters
- {ROOT}.kpp: generated from MCM web
- mcm\_{ROOT}\_mass.txt: table of masses and SMILES strings (downloaded as mcm\_subset\_mass.txt)

Simulations:

- photolysis.txt: input for kpp_constants.f90
- input\_time.txt: time in units of seconds.

		{TSTART}
		{DURATION}
		{DT}

- input\_temp.txt: temperature and conversion factor (ppb to molec/cm^3) (untested)

		{TEMP}
		{CFACTOR}

- cgas\_init.def: initial gas-phase concentrations (in ppb) in equation form as you would write in the .def file

		{COMPOUND1} = {PPB1}
		{COMPOUND1} = {PPB2}
		...

- input\_partitioning.txt: `M0` is the initial aerosol concentration; `PARTITION\_ON` is 1 (or 0 for no partitioning); `INTEGRATORCHECK` determines whether additional diagnostics are run (0=off, 1=on); `MINCONC` is the value (in ppb) at which minimum concentrations in gas and aerosol phases are maintained.

		{M0}
		{PARTITION_ON}
		{INTEGRATORCHECK}
                {MINCONC}

- molefrac\_init.txt: `IND` is the organic compound index and `a0` is the initial mole fraction; the first line is a label (for "harvest_parms.py") preceded by `#`

		{#COMMENT}
		{IND1} {a0(1)}
		{IND2} {a0(2)}


## Instructions

There are four main executable python scripts. The first should be run in the compound folder, and the rest in the simulation folder.

* search\_struct.py: generates SIMPOL and FTIR group tables; also property tables
* build\_dual.py: construct MCM/KPP executables for gas and aerosol simulations
* execute\_dual.py: run executables for a given folder of input parameters ("runpath")
* harvest\_parms.py: harvest parameters from one or more runpath folders

Add kppaermod/ to the list of paths in which executable are searched:

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

* `ROOT`: label for KPP
* [optional] `PROGPATH`: path to aprl-structsearch if not on executable path


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
$ build_dual.py {ROOT} {CPATH} {--skipbuild} {--onlygas} {--onlytotal}
```
Arguments:

* `ROOT`: label for KPP
* `CPATH`: path to compounds directory
* [optional] `--skipbuild`: will not run kpp again but only use output of kpp (in "kppbuild/") to generate "exec\_gas/" and "exec\_total/". Default is to build.
* [optional] `--onlygas`: will not run kpp again but only use output of kpp (in "kppbuild/") to generate "exec\_gas/". Default is to generate "exec\_total".
* [optional] `--onlytotal`: will not run kpp again but only use output of kpp (in "kppbuild/") to generate  "exec\_total/". Default is to generate "exec\_total".

Example usage:
```
$ build_dual.py apinene ../../compounds/apinene_1
```

Outputs (in working directory):

* contents of {CPATH} are copied here
* kppbuild/: kpp output
* exec\_gas/: executable for gas-phase simultions
* exec\_total/: executable for total-phase simulations

### Execute gas-phase and gas+aerosol simulations

*Run in simulation directory (e.g., "simulations/apinene\_1/").*

Command:
```
$ exec_dual.py {ROOT} {RUNPATH} {MODE}
```

Arguments:

* `ROOT`: label for KPP
* `RUNPATH`: name of input folder
* [optional] `MODE`: one of "gas,total", "gas", or "total" (without quotes). Default is "gas,total"

Example usage:
```
$ exec_dual.py apinene run_varyvoc_001
```

Outputs (in run directories):

* runpath/gas/{ROOT}.dat
* runpath/gas/{ROOT}\_formatted.txt
* runpath/total/{ROOT}.dat
* runpath/total/{ROOT}\_formatted.txt
* runpath/total/{ROOT}\_aer.dat
* runpath/total/{ROOT}\_aer\_formatted.txt

### Build executables for gas-phase only ("gas") and gas+aerosol ("total") simulations

*Run in simulation directory (e.g., "simulations/apinene\_1/")*.

Command:
```
$ harvest_parms.py {ROOT} {RUNPATH}
```

Arguments:

* `ROOT`: label for KPP
* [optional] `RUNPATH`: zero or more paths. if omitted, all paths beginning with "run\_" will be harvested

Example usage:
```
$ harvest_parms.py apinene
```

Outputs (in working directory):

* parameter_table.csv

