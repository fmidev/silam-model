#!/bin/bash
rm virtualspecies.ini virtualspecies.spc
touch virtualspecies.ini virtualspecies.spc
for i in `seq 1 200`;
do
    echo $i
    istring=$( printf '%03d' $i )
    cp LISTsubstanceTemplate.txt tmp.txt
    sed -i -e "s/xxx/$istring/g" tmp.txt
    cat tmp.txt >> virtualspecies.ini
    printf "    VirS%03d     = IGNORE;       {Virtual species}\n" $i >> virtualspecies.spc
    #echo $spcstring
done
rm tmp.txt
