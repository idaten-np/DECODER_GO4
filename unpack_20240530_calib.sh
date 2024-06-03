#!/bin/bash
set -e

RunName=Run30_0000_Jenny_calib_Eu3

InputPath=/home/DATA1/20240530
InputFile=${RunName} ##_NNNN.lmd

OutputPath=/home/DATA2/mbsrootfiles/20240530
OutputFile=${RunName} ##_NNNN.root



#### For TAM calib
#rm -f Go4AutoSave.root
#ln -sf set_TamexControl.C~new set_TamexControl.C 
#go4analysis -file ${InputPath}/${InputFile}_00*.lmd -outevt-class TTamex_FullEvent -rate 
#ln -sf set_TamexControl.C~useold set_TamexControl.C 
#### end of TAM calib



ln -sf set_TamexControl.C~useold set_TamexControl.C 
mkdir -p ${OutputPath}

for i in {0001..0008}
do
	go4analysis -file ${InputPath}/${InputFile}_${i}.lmd -outevt-class TTamex_FullEvent -enable-store -store ${OutputPath}/${OutputFile}_${i}.root -rate
done
