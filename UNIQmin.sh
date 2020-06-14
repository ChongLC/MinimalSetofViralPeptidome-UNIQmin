#!/bin/bash
#$ -V

dir=/backup/user/ext/perdana/lichuin/cdhitObj2/testingStitchScript/
cd $dir
python p1.py
wait 
python p2.py
wait
python p2_2.py
wait
python p3.py
wait 
python p3_2.py
wait 
python p4_1.py
wait 
python p4_2.py
wait 
python p4_3.py
wait
cp seqfileZ.txt fileZ.txt
mkdir match
wait 
python p5.py