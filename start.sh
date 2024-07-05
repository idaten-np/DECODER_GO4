#!/bin/bash
set -x
mbsrun=$1
vmerun=$2

go4analysis -lib ./libGo4UserAnalysis.so -stream 10.32.6.213 -http 5000 -rate -outevt-class TTamex_FullEvent -enable-store -store run${mbsrun}_${vmerun}_online.root
