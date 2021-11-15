#!/bin/bash

# input_Wm_e_50-100_lo.cfg  input_Wp_e_50-100_lo.cfg

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

#input_Wm_e_50-100_lo.cfg

Mtmins="50"

for Mtmin in $Mtmins
do

   if [ $Mtmin -lt 80 ]; then 
	Mtmax=8000;
        echo "START `date` ;"
	echo "./makeinputLO.sh ${procId} $Mtmin $Mtmax ;"
	echo "../src/mcsanc input_${NameString}_${Mtmin}-${Mtmax}_lo.cfg >& 1126_${NameString}_${Mtmin}.log ;"
	cp -r /your/Working/Directory/ewparam.cfg .
	source makeinputLO.sh ${procId} $Mtmin $Mtmax ; ../src/mcsanc input_${NameString}_${Mtmin}-${Mtmax}_lo.cfg >& 1126_${NameString}_${Mtmin}.log ;
	echo "END   `date` ; ---> ${NameString}_${Mtmin}-${Mtmax} LO Finish !!"
   fi

done
