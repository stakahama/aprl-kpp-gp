#!/bin/sh

### example:
# $ /path/to/build.sh apinene kpp_constants_dark.f90 /path/to/compounds 
# $ /path/to/build.sh {ROOT} {KPPC} {CPATH}

### define paths on aprlpc1.epfl.ch:
# export KPP_HOME=/opt/kpp       #KPP
# export PATH=$PATH:/opt/kpp/bin #Executable PATHs
export PYTHONPATH=$PYTHONPATH:/opt/anaconda/lib/python2.7/site-packages/openbabel-1.8.1-py2.7-linux-x86_64.egg #pybel installation

# INPUT PARAMETERS #################################################################
ROOT=$1		#apinene
KPPC=$2		#../constantsfile/kpp_constants_dark.f90
CPATH=$3	#'../compounds'
####################################################################################

### get top directory
MAINPATH=`dirname $0`

### STEP 1 generate standard MCM code with kpp:
cp -pv $CPATH/${ROOT}/* .
cp -pv $MAINPATH/modules/* .
python $MAINPATH/scripts/editkpp_v3.py $ROOT
kpp ${ROOT}.def

### STEP 2 aerosol module:
cp -pv $KPPC kpp_constants.f90
python $MAINPATH/scripts/extract_org_indices_and_molmasses.py ${ROOT} mcm_${ROOT}_mass.txt
python $MAINPATH/scripts/kpp_edit_initialize.py ${ROOT} 
python $MAINPATH/scripts/kpp_edit_main.py ${ROOT}
python $MAINPATH/scripts/kpp_edit_makefile.py ${ROOT}
python $MAINPATH/scripts/kpp_generate_AERmodules.py ${ROOT}
python $MAINPATH/scripts/kpp_generate_SIMPOLGroups_extended.py ${ROOT} mcm_${ROOT}_mass.txt
make -f Makefile_${ROOT} 

## STEP 3 (I doubt you need all of these files):
cp -pv $MAINPATH/extra/* .
