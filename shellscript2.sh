#!/bin/bash
# EXPLANATIONS:
# start in a working directory ROOT/, created with the shell_script1_...
#
# where root is your MCM compound. Define the same name in this script below:
# Define also the series of input file numbers, for which you have input files in ../ROOT_template/in_outpus/
# e.g. INPUTS=('00' '01'); to run simulations with input00.txt and input01.txt
#
# CAREFUL THE DATA IN THE INPUT FILES WILL OVERRIDE THE THE ROOT.DEF-FILE. MAKE SURE YOU DO NOT CHANGE IT UNINTENTIONALLY.
#
# INPUT PARAMETERS #################################################################
#define root here:
ROOT=$1
NUMBERS=$2
IFS=',' read -d '' -ra TMPARRAY < <(printf '%s\0' "$NUMBERS")
declare -a INPUTS
for i in "${TMPARRAY[@]}"; do
    INPUTS[i]=`printf %03d $i`
done
echo 'root:' $ROOT '      input_files:' ${INPUTS[*]}
# note: we need to define the path to the pybel installation. on aprlpc1.epfl.ch this is:
export PYTHONPATH=$PYTHONPATH:/opt/anaconda/lib/python2.7/site-packages/openbabel-1.8.1-py2.7-linux-x86_64.egg
####################################################################################

# outputanalysisFabian.py will be called with the arguments: root, 'aer' or 'gas', file_number, 'all', 'mean', 'subset', number_of_values_to_plot

# N_SUBSET=8 			# here, you can define how many datapoints will be extracted/computed for mean and subset arguments 
N_SUBSET=${3:-8} # default value is 8

MAINPATH=`dirname $0`/extra

for i in "${INPUTS[@]}"; do
    runpath='run_'${i}
    ln -svf $PWD/mcm_${ROOT}_mass.txt $runpath/total
    cd $runpath/total
    mkdir mean
    mkdir subset
    mkdir all
    python ${MAINPATH}/outputanalysisFabian_.py ${ROOT} gas ${i} mean ${N_SUBSET}
    python ${MAINPATH}/outputanalysisFabian_.py ${ROOT} aer ${i} mean ${N_SUBSET}
    mv ${ROOT}_*_mean_n${N_SUBSET}* mean/
    python ${MAINPATH}/outputanalysisFabian_.py ${ROOT} gas ${i} subset ${N_SUBSET}
    python ${MAINPATH}/outputanalysisFabian_.py ${ROOT} aer ${i} subset ${N_SUBSET}
    mv ${ROOT}_*_subset_n${N_SUBSET}* subset/
    python ${MAINPATH}/outputanalysisFabian_.py ${ROOT} gas ${i} all 0
    python ${MAINPATH}/outputanalysisFabian_.py ${ROOT} aer ${i} all 0
    mv ${ROOT}_*_all_n0* all/
done
