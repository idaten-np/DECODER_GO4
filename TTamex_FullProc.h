
//-------------------------------------------------------------
//        Go4 Release Package v3.03-05 (build 30305)
//                      05-June-2008
//---------------------------------------------------------------
//   The GSI Online Offline Object Oriented (Go4) Project
//   Experiment Data Processing at EE department, GSI
//---------------------------------------------------------------
//
//Copyright (C) 2000- Gesellschaft f. Schwerionenforschung, GSI
//                    Planckstr. 1, 64291 Darmstadt, Germany
//Contact:            http://go4.gsi.de
//----------------------------------------------------------------
//This software can be used under the license agreements as stated
//in Go4License.txt file which is part of the distribution.
//----------------------------------------------------------------
#ifndef TUNPACKPROCESSOR_H
#define TUNPACKPROCESSOR_H

#define WR_TIME_STAMP     1   // white rabbit time stamp is head of data
#ifdef WR_TIME_STAMP
#endif // WR_TIME_STAMP

#define MAX_SPEZIAL 1000000000

#define VETO_EVT 1

#define IDATEN_MONITOR 1
#ifdef IDATEN_MONITOR
#define TREND_INTV 60e9 /// ns
#define TREND_N 30
#endif // IDATEN_MONITOR

#define ONLINE_CALIB 1
#ifdef ONLINE_CALIB
#define MAX_ENERGY_OI 2000 // keV
#define Npar 4
#endif // ONLINE_CALIB

#ifdef WR_TIME_STAMP
#define SUB_SYSTEM_ID      0x0100
#define TS__ID_L16         0x03e1
#define TS__ID_M16         0x04e1
#define TS__ID_H16         0x05e1
#define TS__ID_X16         0x06e1
#endif // WR_TIME_STAMP

#define STATISTIC 200000

#define DUMP_BAD_EVENT 1

#define COARSE_CT_RANGE  0x800  // 11 bits

#define MAX_SSY        1                // maximum number of sub-systems (readout pcs in nxm system)
#define MAX_SFP        2
#define MAX_TAM        7                // maximum febex/tamex per sfp

#define MAX_CHA_old_INPUT 33                // A) maximum physical input channels per module. must be modulo 4
#define MAX_CHA_old       MAX_CHA_old_INPUT * 2 // B) leading egdes + trailing edges + qtc trailing edges
#define MAX_CHA_phy 16
#define MAX_CHA_tam MAX_CHA_phy *2 +1 // last for TTS

//#define N_CAL_EVT               (ULong64_t) 200000
#define N_CAL_EVT               (ULong64_t) 36*50000

#define N_PHY_TREND_PRINT       (ULong64_t) 1000000

//#define N_DELTA_T   400000
#define N_DELTA_T   1000*1000
#define N_BIN       100000
#define N_TIM       10000

#define N_TR_BINS   100000  
#define N_COARSE    30

#define CYCLE_TIME    (Double_t) 5000

//#define TRIG_WIN_SIZE      200     // in clock cycles 
#define HITPAT_CT_RANGE    20

#define N_BIN_T    600
#define RESET_VAL -100000

#include "TGo4EventProcessor.h"

//#include "TTamex_FullEvent.h"

class TTamex_FullParam;
class TTamex_FullEvent;
class TGo4Fitter;

class TTamex_FullProc : public TGo4EventProcessor {
	public:
		TTamex_FullProc() ;
		TTamex_FullProc(const char* name);
		virtual ~TTamex_FullProc() ;

		Bool_t BuildEvent(TGo4EventElement* target); // event processing function

	private:
		TGo4MbsEvent  *fInput; //!
		TTamex_FullEvent* fOutput; //!

		TTamex_FullParam* fPar;

		Bool_t fCalibrationDone;// flag if calibration is ready

		TH1   *h_box[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_old];  // box histogram in SFP id / TAMEX id / CHANNEL nr coordinates

		TH1   *h_err_box[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_old];  // box histogram in SFP id / TAMEX id / CHANNEL nr coordinates

		TH1   *h_tim_2[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_tam];
		TH1   *h_cct[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_tam];
		TH1   *h_sum_2[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_tam];  // sum histogram in SFP id / TAMEX id / CHANNEL nr coordinates

		TH1   *h_p_sum_ab;
		TH2   *h_p_tota_vs_a;
		TH2   *h_p_totb_vs_b;
		TH2   *h_p_diff_ba_sum_ab;                

#ifdef IDATEN_MONITOR
		TH2 *h2_PCHA_STOT[MAX_SSY][MAX_SFP][MAX_TAM];
		TH2 *h2_PCHA_FTOT[MAX_SSY][MAX_SFP][MAX_TAM];

		TH1 *h1_STOT[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_phy];
		TH1 *h1_FTOT[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_phy];
		TH2 *h2_STOT_FTOT[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_phy];
		TH2 *h2_STOT_FTle[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_phy];
		TH1 *h1_FTle[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_phy];
		TH2 *h2_trend_STOT[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_phy];

		TH1 *h1_Multiplicity[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_tam][3]; // 0 for trailing, 1 for leading, 2 for TOT
		TH1 *h1_Multiplicity_LaBr3;
#endif // IDATEN_MONITOR

#ifdef ONLINE_CALIB
		Double_t par_f1_STOT_Energy[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_phy][Npar] = {0};
#ifdef IDATEN_MONITOR
		TH2 *h2_PCHA_Energy[MAX_SSY][MAX_SFP][MAX_TAM];
		TH1 *h1_Energy[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_phy];
		TH2 *h2_Energy_FTle[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_phy];

		TH1 *h1_Energy_LaBr3;
		TH2 *h2_Energy_FTle_LaBr3;


#endif // IDATEN_MONITOR
#endif // ONLINE_CALIB

		TGo4Picture      *fPicture;

		ClassDef(TTamex_FullProc,1)
};

static  UInt_t l_err_catch = 0;
static  UInt_t l_prev_err_catch = 0;
static  UInt_t l_err_ssy [MAX_CHA_old];
static  UInt_t l_err_sfp [MAX_CHA_old];
static  UInt_t l_err_tam [MAX_CHA_old];
static  UInt_t l_err_cha [MAX_CHA_old];
static  UInt_t l_prev_err_ssy [MAX_CHA_old];
static  UInt_t l_prev_err_sfp [MAX_CHA_old];
static  UInt_t l_prev_err_tam [MAX_CHA_old];
static  UInt_t l_prev_err_cha [MAX_CHA_old];
static  UInt_t l_num_err;
static  UInt_t l_prev_num_err;
#ifdef WR_TIME_STAMP
static	ULong64_t l_wr_ts00;
static ULong64_t test_good=0;
static ULong64_t test_bad1=0;
static ULong64_t test_bad2=0;
static ULong64_t test_bad3=0;
#endif // WR_TIME_STAMP

#endif //TUNPACKPROCESSOR_H


//----------------------------END OF GO4 SOURCE FILE ---------------------


#ifndef __TPHIT
#define __TPHIT
struct TPHit
{
      UInt_t 	SSY;
      UInt_t 	SFP;
      UInt_t 	TAM;
      Int_t  	PCHA;
      Double_t  STOT;
      Double_t  STle;
      Double_t  FTOT;
      Double_t  FTle;
      Double_t  CalE;
      Double_t  TTS;
};
#endif // __TPHIT


