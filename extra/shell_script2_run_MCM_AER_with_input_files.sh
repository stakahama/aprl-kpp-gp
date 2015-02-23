#!/bin/sh
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
export ROOT=tm135b
declare -a INPUTS=('00' '10' '01' '11' '02' '12' '03' '13');

echo 'root:' $ROOT '      input_files:' ${INPUTS[*]}
# note: we need to define the path to the pybel installation. on aprlpc1.epfl.ch this is:
export PYTHONPATH=$PYTHONPATH:/opt/anaconda/lib/python2.7/site-packages/openbabel-1.8.1-py2.7-linux-x86_64.egg
####################################################################################




### STEP 3 Run the codes and save outputs in ../results_${ROOT}:

cp -r ../${ROOT}_template/inputs/ in_outputs/
mkdir ../results_${ROOT}

for i in "${INPUTS[@]}"
do
./${ROOT}.exe $i
cp ${ROOT}.dat ${ROOT}_$i.dat
cp ${ROOT}_aer.dat ${ROOT}_aer_$i.dat
python process.py ${ROOT}_$i.dat
python process.py ${ROOT}_aer_$i.dat
cp ${ROOT}_${i}_formatted.txt ../results_${ROOT}/${ROOT}_${i}_formatted.txt 
cp ${ROOT}_aer_${i}_formatted.txt ../results_${ROOT}/${ROOT}_aer_${i}_formatted.txt 
done

echo 'Finished Simulations for root:' $ROOT '      input_files:' ${INPUTS[*]} '          ' 




### STEP 4 Analyze output:

# prepare for outputanalysisFabian.py, which will be executed from the results_${ROOT}/folder
cd ..
cp ${ROOT}/mcm_${ROOT}_mass.txt results_${ROOT}/
cp ${ROOT}/outputanalysisFabian.py results_${ROOT}/
#cp ${ROOT}/${ROOT}_SMILES.txt results_${ROOT}/
cp ${ROOT}/matmul.py results_${ROOT}/
cp ${ROOT}/simpolmod.py results_${ROOT}/

# outputanalysisFabian.py will be called with the arguments: root, 'aer' or 'gas', file_number, 'all', 'mean', 'subset', number_of_values_to_plot

export N_SUBSET=8 			# here, you can define how many datapoints will be extracted/computed for mean and subset arguments 

cd results_${ROOT}/
mkdir mean
mkdir subset
mkdir all
for i in "${INPUTS[@]}"
do
python outputanalysisFabian.py ${ROOT} gas ${i} mean ${N_SUBSET}
python outputanalysisFabian.py ${ROOT} aer ${i} mean ${N_SUBSET}
mv ${ROOT}_*_${i}_mean_n${N_SUBSET}* mean/
python outputanalysisFabian.py ${ROOT} gas ${i} subset ${N_SUBSET}
python outputanalysisFabian.py ${ROOT} aer ${i} subset ${N_SUBSET}
mv ${ROOT}_*_${i}_subset_n${N_SUBSET}* subset/
python outputanalysisFabian.py ${ROOT} gas ${i} all 0
python outputanalysisFabian.py ${ROOT} aer ${i} all 0
mv ${ROOT}_*_${i}_all_n0* all/
done
