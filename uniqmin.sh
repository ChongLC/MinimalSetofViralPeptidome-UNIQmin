#!/bin/bash
#$ -V

dir=/backup/user/ext/perdana/lichuin/cdhitObj2/testingStitchScript/
cd $dir
python3 U1_KmerGenerator.py
wait 
python3 U2.1_Singletons.py
wait
python3 U2.2_Multitons.py
wait
python3 U3.1_PreQualifiedMinSet.py
wait 
python3 U3.2_UnmatchedSingletons.py
wait 
python3 U4.1_Non-SingletonsDedup.py
wait 
python3 U4.2_Multi-OccurringPreMinSet.py
wait 
python3 U4.3_UnmatchedMulti-Occurring.py 
wait
mkdir minimalSet
cp seqfileZ.txt minimalSet/fileZ.txt
mkdir match
wait 
python3 U5.1_RemainingMinSet.py
wait
python3 U5.2_MinSet.py
