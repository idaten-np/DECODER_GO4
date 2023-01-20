
##################################################
##        convert LMD file to root tree         ##
##################################################

# set environment (just like /path/to/thisroot.sh)
. /cvmfs/eel.gsi.de/bin/go4login


# if you changed Go4 settings you have to
make clean && make
# before go4analysis


# file
go4analysis -file 220602_pdc_bjt_twin_peaks_th_fst_15000_beam_0001.lmd -outevt-class TTamex_FullEvent -enable-store -store new.root -rate

# from stream (if you want to slurp data from server)
go4analysis -stream x86l-127 -outevt-class TTamex_FullEvent -enable-store -store new.root -rate





##################################################
##       configure your ROOT tree leaves        ##
##################################################

edit TTamex_FullProc.h


each pair of entries forms a difference ( t(n+1) - t(n))
and gets stored into a leaf of the root tree.
The root tree leaves start with 0

                ( 0 )  ( 1 )   ( 2 )  ... etc 
#define CHA_ID { 1,34,  2,35,  21,54,  22,55, 3,36, 4,37, 7,40, 8,41, 9,42,10,43,11,44,12,45,13,46,14,47,15,48,16,49,17,50,18,51,19,52,20,53,21,54,22,55,23,56,24,57,25,58,26,59,27,60,28,61,29,62,30,63,31,64,32,65 }


for leaves to represent a TOT its always  n, n+33   (n is LE, n+33 is TE))

for leaves to represent a Leading Edge, create a difference vs ch 0 (trigger)

if you then want to correlate different detectors, ROOT-Draw a difference of their respective LE leaves.

Always take the difference to ch0 on the first Tamex (TAM_ID 0)


##################################################
##           ROOT tree draw examples            ##
##################################################

root -l new.root

// fast slow correlation example 2D hist
root [6] AnalysisxTree->Draw("fTimeDiff[4]:fTimeDiff[5]>>banana(1000,0,2000000,1000,0,500000)","","colz")
(long long) 6205179

// slow channel TOT spectra 1D hist
AnalysisxTree->Draw("fTimeDiff[3]>>abc(1000,0,2000000)")
root [9] AnalysisxTree->Draw("fTimeDiff[3]>>SLW_TOT(1000,0,2000000)")
