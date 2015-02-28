#!/bin/sh
# EXPLANATIONS:
# start in a working directory to be able to acces 'kppAERmodFabian10' with a relative path as ../kppAERmodFabian10
# contents of this working directory: 
# - ROOT_template/before_kpp/mcm_ROOT_mass.txt
# - ROOT_template/before_kpp/ROOT.kpp
# - ROOT_template/before_kpp/ROOT.def
# - ROOT_template/inputs/input00.txt <- with a series of input files
#
# where root is your MCM compound. Define the same name in this script below:

# INPUT PARAMETERS #################################################################
#define root here:
export ROOT=tm135b

echo 'root:' $ROOT

# note: we need to define the path to the pybel installation. on aprlpc1.epfl.ch this is:
export PYTHONPATH=$PYTHONPATH:/opt/anaconda/lib/python2.7/site-packages/openbabel-1.8.1-py2.7-linux-x86_64.egg
####################################################################################




### STEP 1 generate standard MCM code with kpp: 

cp -r ${ROOT}_template/before_kpp/ ${ROOT}_template/after_kpp/
cd ${ROOT}_template/after_kpp/
python ../../../kppAERmodFabian10/kpp_copy_modules.py
printf "%s" "${ROOT}" | python editkpp.py
kpp ${ROOT}.def




### STEP 2 apply kppAERmodFabian10 on a copy of the MCM code:

cd ../..
cp -r ${ROOT}_template/after_kpp/ ${ROOT}
cd ${ROOT}
cp ../../kppAERmodFabian10/opkd* .
python ../../kppAERmodFabian10/extract_org_indices_and_molmasses.py ${ROOT} mcm_${ROOT}_mass.txt
python ../../kppAERmodFabian10/kpp_edit_initialize.py ${ROOT} 
python ../../kppAERmodFabian10/kpp_edit_main.py ${ROOT}
python ../../kppAERmodFabian10/kpp_edit_makefile.py ${ROOT}
python ../../kppAERmodFabian10/kpp_generate_AERmodules.py ${ROOT}
python ../../kppAERmodFabian10/kpp_generate_SIMPOLGroups_extended.py ${ROOT} mcm_${ROOT}_mass.txt
make -f Makefile_${ROOT} 

