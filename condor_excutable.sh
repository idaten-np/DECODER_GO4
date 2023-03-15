dir=`pwd`
idir=$1
logdir=$dir/$idir/log
mkdir -p $idir
mkdir -p $logdir
outfile=$2_`date +"%Y%m%d"`_$3
subfile=$outfile.sub

filename=$3
inputpath=/home/gamma/KHALA_TPTAMEX_atKU/lmd
inputfile=$filename.lmd
outputfile=$filename.root

cat > $subfile << EOF
Executable   = /home/eaw0224/Go4/go4-6.2.0_install/bin/go4analysis
Arguments    = -file $inputfile -outevt-class TTamex_FullEvent -enable-store -store $outputfile
Initialdir   = $dir/$idir
Log          = $logdir/$outfile.log
Error        = $logdir/$outfile.err
Output       = $logdir/$outfile.out
Input        = 
should_transfer_files   = YES
transfer_input_files = $inputpath/$inputfile, TTamex_FullProc.h, G__Tamex_Full.cxx, TTamex_FullParam.o, TTamex_FullEvent.o, TTamex_FullProc.o, G__Tamex_Full.o, libGo4UserAnalysis.so, libGo4UserAnalysis_rdict.pcm, libGo4UserAnalysis.rootmap, set_TamexControl.C, Go4AutoSave.root
transfer_output_files=$outputfile
Universe     = vanilla
GetEnv       = True
Request_cpus = 1
Notification = Never
#Requirements = machine=="kunpl-node10"
#Notify_user  = CHANGE_THIS_TO_YOUR_EMAIL
Queue
EOF


condor_submit $subfile
rm $subfile

