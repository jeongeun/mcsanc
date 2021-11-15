#!/bin/bash

TopDir=`pwd`
export SCRAM_ARCH=slc6_amd64_gcc630
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
cd /hcp/data/data02/jelee/FEWZ/CMSSW_10_1_9
eval `scramv1 runtime -sh`
cd -
export dirLHAPDF=/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/lhapdf/6.2.1
export PATH=$dirLHAPDF/bin:$PATH
export LD_LIBRARY_PATH=$dirLHAPDF/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$dirLHAPDF/lib/python2.7/site-packages:$PYTHONPATH

#exmaple : source makeAscii.sh Wm 0 1 or Wp 1 

##mcsanc-DY_CC_LO_EW_Wp_e_mt_800-900_13TeV-output.txt
Masses="50-8000"
#100-200
#200-300
#300-400
#400-500
#500-600
#600-700
#700-800
#800-900
#900-1000
#1000-2000
#2000-3000
#3000-4000
#4000-5000
#5000-6000
#6000-7000
#7000-8000"
# PDFs="NNPDF31_nnlo_as_0118_nf_4"
# PDFs="NNPDF31_lo_as_0130"

Mode=${1}  # Wm or Wp
Order=${2}
OrderString=""
if [ "$Order" == "0" ]; then
        OrderString="LO"
elif [ "$Order" == "1" ]; then
        OrderString="NLO"
fi
for mass in $Masses
do
MainString="mcsanc-DY_CC_${OrderString}_EW_${Mode}_e_mt_${mass}_13TeV-output"
Output="${MainString}.txt"

IntLineNum=`grep -n "histogram"  ${Output} | cut -f1 -d: | head -1`
IntLineNum=`expr ${IntLineNum} + 2`

Num=`grep -n "histograms:"  ${Output} | cut -f1 -d: | tail -1`
HistLineNum=`expr ${Num} + 2`

echo " cat ${Output} | head -${IntLineNum} | tail -3   >& ${OrderString}_EW_${Mode}_e_m34_${mass}_int.dat  "
cat ${Output} | head -${IntLineNum} | tail -2   >> ${OrderString}_EW_${Mode}_e_m34_int.dat ; 

echo " sed -n "${HistLineNum},\$p" ${Output}  >& ${OrderString}_EW_${Mode}_e_m34_${mass}_hist.dat  ;"
sed -n "${HistLineNum},\$p" ${Output}  >> ${OrderString}_EW_${Mode}_e_m34_hist.dat  ;

echo "${mass}" >> ${OrderString}_EW_${Mode}_e_m34_int.dat ;
chmod +x ${OrderString}_EW_${Mode}_e_m34_int.dat  ;
chmod +x ${OrderString}_EW_${Mode}_e_m34_hist.dat ;
done

##cat ${Output} | head -${IntLineNum} | tail -2   >& NLO_EW_${Mode}_e_m34_${Mi}-${Mf}_int.dat  
##cat ${Output} | sed -n '${HistLineNum},\$p'  >& NLO_EW_${Mode}_e_m34_${Mi}-${Mf}_hist.dat  
## chmod +x NLO_EW_${Mode}_e_m34_${Mi}-${Mf}_int.dat  
## chmod +x NLO_EW_${Mode}_e_m34_${Mi}-${Mf}_hist.dat 

#NNLO_W_M80-100_Wm_xsec.dat
# mv  ${OrderString}_W_M${3}-${4}_${1}_\*.dat ./data

### ./finish.sh run_W_M1-20_Wm_2 NNLO.run_W_M1-20_Wm_2.dat
### mv NNLO.run_W_M1-20_Wm_2.dat NNLO_W_M1-20_Wm.dat
### cat NNLO_W_M1-20_Wm.dat | head -68 | tail -2 >& NNLO_W_M1-20_Wm.xsec.dat
### cat NNLO_W_M1-20_Wm.dat | head -141 | tail -21 >&  NNLO_W_M1-20_Wm.hist.dat
