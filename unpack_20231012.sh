#!/bin/bash
set -e

RunName=Run09_2_0095_Bi207_fadc_tpt_sync

InputPath=/home/DATA1/20240403
InputFile=${RunName} ##_NNNN.lmd

OutputPath=/home/shared/20240403
OutputFile=${RunName} ##_NNNN.root



#### For TAM calib
#rm -f Go4AutoSave.root
#ln -sf set_TamexControl.C~new set_TamexControl.C 
#go4analysis -file ${InputPath}/${InputFile}_0001.lmd -outevt-class TTamex_FullEvent -enable-store -store new.root -rate 
#ln -sf set_TamexControl.C~useold set_TamexControl.C 
#### end of TAM calib



ln -sf set_TamexControl.C~useold set_TamexControl.C 
mkdir -p ${OutputPath}

for i in {0007..0066}
do
	go4analysis -file ${InputPath}/${InputFile}_${i}.lmd -outevt-class TTamex_FullEvent -enable-store -store ${OutputPath}/${OutputFile}_${i}.root -rate
done
