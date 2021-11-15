#!/bin/bash

# input_Wm_e_50-100_nlo.cfg  input_Wp_e_50-100_nlo.cfg

procId=${1}

if [ "$procId" == "" ]; then
        echo "usage: $0 processId  // processId = W+(e+ve: 101, m+vm: 102)  W-( e-ve~:-101, m-vm~: -102)";
      break;
fi

NameString=""
if [ "$procId" == "101" ]; then
        NameString="Wp_e"
elif [ "$procId" == "-101" ]; then
        NameString="Wm_e"
elif [ "$procId" == "102" ]; then
        NameString="Wp_mu"
elif [ "$procId" == "-102" ]; then
        NameString="Wm_mu"
elif [ "$procId" == "103" ]; then
        NameString="Wp_tau"
elif [ "$procId" == "-103" ]; then
        NameString="Wm_tau"
fi

#Mtmins="50 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000"
Mtmins="50"

for Mtmin in $Mtmins
do

   if [ $Mtmin -lt 80 ]; then 
	Mtmax=`expr $Mtmin + 7950` ;
        echo "START `date` ;"
	echo "./makeinputNLO.sh ${procId} $Mtmin $Mtmax ;"
	echo "../src/mcsanc input_${NameString}_${Mtmin}-${Mtmax}_nlo.cfg >& 0422_${NameString}_${Mtmin}.log ;"
	cp /your/working/directory/ewparam.cfg .
	source makeinputNLO.sh ${procId} $Mtmin $Mtmax ; ../src/mcsanc input_${NameString}_${Mtmin}-${Mtmax}_nlo.cfg >& 0422_${NameString}_${Mtmin}.log ;
	echo "END   `date` ; ${NameString}_${Mtmin}-${Mtmax} NLO Finish =="
   fi

done

