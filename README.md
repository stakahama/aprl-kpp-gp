APRL KPP G/P module [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.160810.svg)](https://doi.org/10.5281/zenodo.160810)
===

This program generates a gas-phase chemical kinetic model using the Kinetic Pre-Processor (KPP) [1][2] with the Master Chemical Mechanism (MCM) [3], and adds dynamic gas/particle (G/P) partitioning with vapor pressure estimation from the SIMPOL.1 group contribution model [4]. Further details are provided by Ruggeri *et al.* [5].

1. http://dx.doi.org/10.5194/acp-6-187-2006
2. https://github.com/barronh/kpp
3. http://mcm.leeds.ac.uk/MCM
4. http://dx.doi.org/10.5194/acp-8-2773-2008
5. http://dx.doi.org/10.5194/acp-16-8729-2016

This program is released under the GNU Public License v3.0 (LICENSE_GPLv3.txt). If used, please include a citation to our manuscript:

> Ruggeri, G., Bernhard, F. A., Henderson, B. H., and Takahama, S.: Model–measurement comparison of functional group abundance in α-pinene and 1,3,5-trimethylbenzene secondary organic aerosol formation, *Atmos. Chem. Phys.*, 16, 8729-8747, doi:10.5194/acp-16-8729-2016, 2016.

Main features:

* Vapor pressure estimates from MCM species and functional groups defined in APRL-SSP (https://github.com/stakahama/aprl-ssp). [Python]
* Gas-phase simulation code generated by KPP is modified to implement dynamic partitioning via operator splitting. [Python/Fortran]
* Generates 1) original gas-phase only and 2) gas-phase with G/P partitioning module to run with same inputs.

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

Note that the user will provide "cgas\_init.def"; a corresponding "cgas\_init.txt" file to be read by the Fortran program will be generated by "exec\_dual.py" described below.

The main objective is to build a program for a fixed mechanism (set of chemical reactions, species) to simulate over a range of temperatures, concentrations, and timesteps. After generating the gas and total (gas+aerosol) simulation models ("exec\_gas/{ROOT}.exe" or "exec\_total/{ROOT}.exe" in each simulation subdirectory), parameters can be changed through input files for various simulations ("run\_{DDD}/"). Note that runs using "input\_temp.txt" and "cgas\_init.txt" are untested and should be against a reference simulation.

### File descriptions

Mechanism information:

- {ROOT}.def: combines organic and inorganic kpp files; specifies initial concentrations, temperature, and time parameters
- {ROOT}.kpp: generated from MCM web
- mcm\_{ROOT}\_mass.txt: table of masses and SMILES strings (downloaded as mcm\_subset\_mass.txt)

Simulations:

- photolysis.txt: input for kpp_constants.f90
- input\_time.txt: time in units of seconds. Note that when partitioning is turned on, the operators are coupled as S1(DT)oS2(DT) so 2*DT is a full timestep for gas-phase chemistry + partitioning.

		{TSTART}
		{DURATION}
		{DT}

- input\_temp.txt: temperature and conversion factor (ppb to molec/cm^3) (untested)

		{TEMP}
		{CFACTOR}

- cgas\_init.def: initial gas-phase concentrations (in ppb) in equation form as you would write in the .def file

		{COMPOUND1} = {PPB1}
		{COMPOUND2} = {PPB2}
		...
where `COMPOUND1`, `COMPOUND2` are names of species. There will be a corresponding file generated by the program called cgas\_init.txt, and this will be the file read in by the Fortran program (uses species indices rather than species names).

- input\_partitioning.txt: `M0` is the initial aerosol concentration in micrograms per cubic meter; `PARTITIONING_MODE` is 0 for no partitioning, 1 for instantaneous (equilibrium) partitioning, and 2 for dynamic partitioning using LSODE; `ABSORPTIVE_MODE` is whether/when to turn on absorptive partitioning (see below); `INTEGRATORCHECK` determines whether additional diagnostics are run (0=off, 1=on); `MINCONC` is the value (in ppb) at which minimum concentrations in gas and aerosol phases are maintained; `MF` is a DLSODE option which controls the integration (10=Nonstiff, no Jacobian required; 21=User-supplies Jacobian-generating function (default); 22=Jacobian is internally generated). `MINCONC`=0 and `MF`=22 is recommended. `NCONC` is the fixed number concentration [m^-3] of particles and `DIAM_SEED` is the seed diameter [m] (enter 0.E0 if no seed).

		{M0}
		{PARTITIONING_MODE}
		{ABSORPTIVE_MODE}
		{INTEGRATORCHECK}
		{MINCONC}
		{MF}
		{NCONC}
		{DIAM_SEED}

`ABSORPTIVE_MODE` options:

* `0`: begin absorptive partitioning immediately.
* `1`: begin absorptive partitioning when C<sub>OA</sub> > 0 (does not use information about C<sub>OA,init</sub>).
* `2`: begin absorptive partitioning when C<sub>OA</sub> > C<sub>OA,init</sub>

Note that for the "extra solvent" simulation, set `ABSORPTIVE_MODE` to `0` and do not provide a molefrac\_init.txt file.

- molefrac\_init.txt: `IND` is the organic compound index and `a0` is the initial mole fraction; the first line is a label (for "harvest_parms.py") preceded by `#`

		{#COMMENT}
		{IND1} {a0(1)}
		{IND2} {a0(2)}
This input file is best generated by a a script (~/git/projects/aprl-kpp-gp/postprocess/a0\_initialize.R). If molefrac\_init.txt is not provded and `ABSORPTIVE_MODE` is `0` in input\_partitioning.txt -> "extra solvent" mode; `ABSORPTIVE_MODE` is `1` in input\_partitioning.txt -> "infinite sink" assumption (only until C<sub>OA</sub> > C<sub>OA,init</sub>, which generally occurs in the first time step).

## Instructions

There are four main executable python scripts. The first should be run in the compound/ folder, and the rest in the simulations/ folder.

* search\_struct.py: generates SIMPOL and FTIR group tables; also property tables
* build\_dual.py: construct MCM/KPP executables for gas and aerosol simulations
* execute\_dual.py: run executables for a given folder of input parameters ("runpath")
* harvest\_parms.py: harvest parameters from one or more runpath folders

Add aprl-kpp-gp/ to the list of paths in which executable are searched:

```
$ export PATH=~/git/projects/aprl-kpp-gp:$PATH
```

### Calculate group abundances for each compound


*Run in compounds directory (e.g., "compounds/apinene\_1")*.


This uses [aprl-structsearch](https://bitbucket.org/stakahama/aprl-structsearch). Note that according to the instructions for aprl-structsearch, you should add the path of this program to the `PATH` environmental variable also:
```
$ export PATH=~/git/projects/aprl-structsearch:$PATH
```

I have created a script, search_struct.py, in aprl-kpp-gp to facilitate generation of SIMPOL and FTIR groups, and also the vapor pressures at 298.15K and 358.15K (60 degrees C) using the aprl-structsearch program.

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
* runpath/gas/{ROOT}\_formatted.csv
* runpath/total/{ROOT}.dat
* runpath/total/{ROOT}\_formatted.csv
* runpath/total/{ROOT}\_aer.dat
* runpath/total/{ROOT}\_aer\_formatted.csv

Note that "exec\_dual.py" will also create a "cgas\_init.txt" file in the run directory for the Fortran program to read; each time the script is executed, "cgas\_init.txt" will be overwritten by the translated contents of the user-provided "cgas\_init.def".

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

### [optional] Generate molefrac\_init.txt file.

Command:
```
$ ~/git/projects/aprl-kpp-gp/postprocess/a0_initialize.R {ROOT} {TYPE} {RUNPATH}
```

Arguments:

* `ROOT`: label for KPP
* `TYPE` can be one of:
  * purecomponent
  * equalcomponent
  * initialequilibrium
  * gasphasecomp
  * recycledseed
* `RUNPATH`: name of input folder

Note that `ABSORPTIVE_MODE` in input\_partitioning.txt, in addition to the `TYPE` argument, will also affect the partitioning. For "infinitesink" or "extrasolvent", a0\_initialize.R does not need to be invoked (see explanation for molefrac\_init.txt).

Example usage:
```
$ ~/git/projects/aprl-kpp-gp/postprocess/a0_initialize.R apinene purecomponent run_001
```
### Full example

```
$ export MYPATH=/path/to/MCM/new # change to desired path

$ cd $MYPATH
$ cd compounds/apinene_1
$ search_struct.py apinene

$ cd $MYPATH
$ cd simulations/apinene_1
$ build_dual.py apinene ../../compounds/apinene_1
-> note that simulations/apinene_1/ has been populated (e.g., with kppbuild/, exec_gas/, exec_total/)

$ cd $MYPATH
$ mkdir simulations/apinene_1/run_001
$ cp -p photolysisfiles/dark/photolysis.txt simulations/apinene_1/run_001/
-> also add other input files (e.g., run a0\_initialize.R)

$ exec_dual.py apinene run_001
-> check in run_001/gas/ and run_001/total/ for outputs files

$ harvest_parms.py apinene
-> look at parameter_table.csv
```
