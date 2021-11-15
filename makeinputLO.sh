#!/bin/bash

TopDir=`pwd`
export SCRAM_ARCH=slc6_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_1_0
eval `scramv1 runtime -sh` # -- cmsenv

export dirLHAPDF=/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/lhapdf/6.2.1-gnimlf2    #/share/LHAPDF/NNPDF31_nlo_as_0118_luxqed/NNPDF31_nlo_as_0118_luxqed_0000.dat
export PATH=$dirLHAPDF/bin:$PATH
export LD_LIBRARY_PATH=$dirLHAPDF/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$dirLHAPDF/lib/python2.7/site-packages:$PYTHONPATH

cd $TopDir

procId=${1}
Mtmin=${2}
Mtmax=${3}
Bin=100
if [ $Mtmin -gt 999 ]; then
	Bin=1;
fi

if [ "$procId" == "" ]; then
        echo "usage: $0 processId Mtmin Mtmax  // processId = W+(e+ve: 101, m+vm: 102)  W-( e-ve~:-101, m-vm~: -102)";
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

Mtbin=`printf ${Mtmin}"-"${Mtmax}`
mainString="DY_CC_LO_EW_${NameString}_mt_${Mtbin}_13TeV"

INPUT="input_${NameString}_${Mtbin}_lo.cfg"

cat << EOF > ${INPUT}
&Process
  processId	= ${procId}
  !run_tag	= 'DY_CC_LO_EW_Wp_mt_50_13TeV_gf2'
  run_tag	= '${mainString}'  
  sqs0     	= 13000d0
  beams	        = 1,1
  !PDFSet	= 'MSTW2008nlo90cl.LHgrid'
  PDFSet	= 'NNPDF31_nnlo_as_0118_nf_4'
  PDFMember	= 0

  ! electroweak flags
  !               iqed, iew, iborn, ifgg, irun, iph, iho 
  iflew         =   1,   1,    0,     1,    1,   0,   0
  ! qcd flags
  iflqcd        =    0
  ! subtraction scheme 0/1/2 : none/MSbar/DIS
  imsb		= 2
  ! collinear particle recombination in calorimeter
  irecomb       = 1
/

&VegasPar
  relAcc	= 1d-50
  absAcc	= 1d-1
  nStart	= 10000000
  nIncrease	= 10000000
  nExplore	= 300000000
  maxEval	= 2100000000000
  seed	        = 41
  flags 	= 26
/

&KinCuts
  cutName	= 'm34',       'mtr',	'pt3',	'pt4',	'eta3',	'eta4',	'dR'
  cutFlag       = 1             0      0       0       0       0       0
  cutLow	= ${Mtmin}d0,  200d0, 25d0,	25d0,  -2.5d0, -2.5d0, 	0d0
  cutUp		= ${Mtmax}d0,  8d3,	8d3,	8d3,	2.5d0,	2.5d0, 	1d-1
/

! particle numbering 1+2 -> 3+4+5+... ( FIXME )
&FixedBinHist
  fbh_name	= 'm34','mtr','pt34',	'pt3',	'pt4', 'et34', 'et3', 'et4',  'phis','incl', 'for',  'bac'
  fbh_flag	= 1,	0,     0,       0,      0,     0,      0,      0,      0,     0,      0,      0
  fbh_low	= ${Mtmin}d0, 20d0,  0d0,	20d0,	20d0, -2.5d0,  0.0d0,  0.0d0,  0.0d0, 50d0,   50d0,   50d0
  fbh_up	= ${Mtmax}d0,70d0,  200d0,	60d0,	60d0,  2.5d0,  2.5d0,  2.5d0,  0.4d0, 14d3,   200d0,  200d0
  fbh_step	= ${Bin}d0,	1d0,   0.25d0,	1d0,	1d0,   0.1d0,  0.125d0,0.125d0,0.01d0,13.55d3,1d0,    1d0
/

&VarBinHist
  nvbh         	= 12,

  vbh_name(1)	= 'm34',
  vbh_flag(1)	= 0,
  vbh_nbins(1)	= 13,
  !vbh_bins(2,1:14) = 20d0 30d0 50d0 70d0 90d0 110d0 130d0 150d0 200d0 300d0 400d0 500d0 1000d0 1500d0,
  vbh_bins(1,1:14) = 6035d0 7285d0 8793d0 10614d0;

  vbh_name(2)	= 'mtr',
  vbh_flag(2)	= 0,
  vbh_nbins(2)	= 13,
  vbh_bins(2,1:14) = 20d0 30d0 50d0 70d0 90d0 110d0 130d0 150d0 200d0 300d0 400d0 500d0 1000d0 1500d0,

  vbh_name(3)	= 'pt34',
  vbh_flag(3)	= 0,
  vbh_nbins(3)	= 1,
  vbh_bins(3,1:2) = 0.00d0 0.25d0,

  vbh_name(4)	= 'pt3',
  vbh_flag(4)	= 0,
  vbh_nbins(4)	= 11,
  vbh_bins(4,1:12) = 0d0 0.2d0 0.4d0 0.6d0 0.8d0 1d0 1.2d0 1.4d0 1.6d0 1.8d0 2.0d0 2.5d0,

  vbh_name(5)	= 'pt4',
  vbh_flag(5)	= 0,
  vbh_nbins(5)	= 12,
  vbh_bins(5,1:13) = 0d0 0.2d0 0.4d0 0.6d0 0.8d0 1d0 1.2d0 1.4d0 1.6d0 1.8d0 2.0d0 2.5d0 2.6d0,

  vbh_name(6)	= 'et34',
  vbh_flag(6)	= 0,
  vbh_nbins(6)	= 1,
  vbh_bins(6,1:2) = 0.00d0 0.25d0,

  vbh_name(7)	= 'et3',
  vbh_flag(7)	= 0,
  vbh_nbins(7)	= 11,
  vbh_bins(7,1:12) = 0d0 0.21d0 0.42d0 0.63d0 0.84d0 1.05d0 1.37d0 1.52d0 1.74d0 1.95d0 2.18d0 2.50d0,

  vbh_name(8)	= 'et4',
  vbh_flag(8)	= 0,
  vbh_nbins(8)	= 12,
  vbh_bins(8,1:13) = 0d0 0.21d0 0.42d0 0.63d0 0.84d0 1.05d0 1.37d0 1.52d0 1.74d0 1.95d0 2.18d0 2.50d0 2.60d0,

  vbh_name(9)	= 'phis',
  vbh_flag(9)	= 0,
  vbh_nbins(9)	= 1,
  vbh_bins(9,1:2) = 0.00d0 0.25d0,

  vbh_name(10)	= 'incl',
  vbh_flag(10)	= 0,
  vbh_nbins(10)	= 1,
  vbh_bins(10,1:2) = 50d0 8d3

  vbh_name(11)	= 'for',
  vbh_flag(11)	= 0,
  vbh_nbins(11)	= 13,
  vbh_bins(11,1:14) = 20d0 30d0 50d0 70d0 90d0 110d0 130d0 150d0 200d0 300d0 400d0 500d0 1000d0 1500d0,

  vbh_name(12)	= 'bac',
  vbh_flag(12)	= 0,
  vbh_nbins(12)	= 13,
  vbh_bins(12,1:14) = 20d0 30d0 50d0 70d0 90d0 110d0 130d0 150d0 200d0 300d0 400d0 500d0 1000d0 1500d0

/
EOF
