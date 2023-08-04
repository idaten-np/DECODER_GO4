#!/bin/bash
set -e


rm Go4AutoSave.root && ln -sf set_TamexControl.C~new set_TamexControl.C 
go4analysis -file /home/DATA1/test/WR_PEX_run10_75_0001.lmd -outevt-class TTamex_FullEvent -enable-store -store new.root -rate 
ln -sf set_TamexControl.C~useold set_TamexControl.C 
#go4analysis -file /home/DATA1/test/WR_PEX_run10_75_0*.lmd -outevt-class TTamex_FullEvent -enable-store -store new.root -rate
go4analysis -file /home/DATA1/test/WR_PEX_run10_75_0001.lmd -outevt-class TTamex_FullEvent -enable-store -store mbs0075_1.root -rate
go4analysis -file /home/DATA1/test/WR_PEX_run10_75_0002.lmd -outevt-class TTamex_FullEvent -enable-store -store mbs0075_2.root -rate
go4analysis -file /home/DATA1/test/WR_PEX_run10_75_0003.lmd -outevt-class TTamex_FullEvent -enable-store -store mbs0075_3.root -rate
go4analysis -file /home/DATA1/test/WR_PEX_run10_75_0004.lmd -outevt-class TTamex_FullEvent -enable-store -store mbs0075_4.root -rate
go4analysis -file /home/DATA1/test/WR_PEX_run10_75_0005.lmd -outevt-class TTamex_FullEvent -enable-store -store mbs0075_5.root -rate
go4analysis -file /home/DATA1/test/WR_PEX_run10_75_0006.lmd -outevt-class TTamex_FullEvent -enable-store -store mbs0075_6.root -rate

