// N.Kurz, EE, GSI, 18-Jun-2012
//
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

#include <stdlib.h>
#include <sys/time.h>
#include <stdint.h>


#include "TTamex_FullProc.h"
#include "TTamex_FullEvent.h"

#include "Riostream.h"

using namespace std;

#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"
#include "snprintf.h"

#include "TGo4Analysis.h"
#include "TGo4MbsEvent.h"
#include "TGo4WinCond.h"
#include "TGo4PolyCond.h"
#include "TGo4CondArray.h"
#include "TGo4Picture.h"
#include "TTamex_FullParam.h"
#include "TGo4Fitter.h"
#include "TLatex.h"

static  FILE *fd_out;
static  struct timeval s_time;
static  time_t l_time;


static  UInt_t l_err_flg = 0;
static  UInt_t l_err_ct [MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_old];
static  UInt_t l_hit_ct [MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_old];
static  Double_t d_err_rate;

//***********************************************************
TTamex_FullProc::TTamex_FullProc() : TGo4EventProcessor("Proc")
{
	cout << "**** TTamex_FullProc: Create instance " << endl;
}
//***********************************************************
TTamex_FullProc::~TTamex_FullProc()
{
	cout << "**** TTamex_FullProc: Delete instance " << endl;
}
//***********************************************************
// this one is used in standard factory
TTamex_FullProc::TTamex_FullProc(const char* name) : TGo4EventProcessor(name)
{
	cout << "**** TTamex_FullProc: Create instance " << name << endl;

	fPar=dynamic_cast<TTamex_FullParam*> (MakeParameter("TamexControl", "TTamex_FullParam","set_TamexControl.C"));
	fCalibrationDone=kFALSE;

	Text_t chis [256];
	Text_t chead[256];
	Int_t iSSY, iSFP, iTAM, iCHA, iTCHA, iPCHA;

	// Creation of histograms (check if restored from auto save file):
	if(GetHistogram("didi")==0)
	{
		for (iSSY=0; iSSY<MAX_SSY; iSSY++)
		{
			for (iSFP=0; iSFP<MAX_SFP; iSFP++)
			{
				for (iTAM=0; iTAM<MAX_TAM; iTAM++)
				{
					for (iCHA=0; iCHA<MAX_CHA_old; iCHA++)
					{
						sprintf (chis,"Error_Box_Diagrams/SUB%d/SFP%d/TAMEX%02d/Error_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iSSY, iSFP, iTAM, iCHA);
						sprintf (chead,"ERR_BOX diagram");
						h_err_box[iSSY][iSFP][iTAM][iCHA] = MakeTH1 ('I', chis, chead, N_BIN_T, 0, N_BIN_T);
					}
				}
			}
		}

		// box histograms in / sub-system / sfp / tamex / channel  coordinates

		for (iSSY=0; iSSY<MAX_SSY; iSSY++)
		{
			for (iSFP=0; iSFP<MAX_SFP; iSFP++)
			{
				for (iTAM=0; iTAM<MAX_TAM; iTAM++)
				{
					for (iCHA=0; iCHA<MAX_CHA_old; iCHA++)
					{
						sprintf (chis,"Box Diagrams/SUB%d/SFP%d/TAMEX %02d/SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iSSY, iSFP, iTAM, iCHA);
						sprintf (chead,"BOX diagram");
						h_box[iSSY][iSFP][iTAM][iCHA] = MakeTH1 ('I', chis, chead, N_BIN_T, 0, N_BIN_T);
					}
				}
			}
		}      

		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_tam; iCHA++)
		{
			sprintf (chis,"COARSE_COUNTER/SUB%d/SFP%d/TAMEX%02d/COARSE_COUNTER_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"COARSE COUNTER");
			h_cct[iSSY][iSFP][iTAM][iCHA] = MakeTH1 ('I', chis, chead, COARSE_CT_RANGE, 0, COARSE_CT_RANGE);
		}      
		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_tam; iCHA++)
		{
			sprintf (chis,"CALIB/SUB%d/SFP%d/TAMEX%02d/CHA%02d/CALIB TIME_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iCHA, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"CALIB TIME");
			h_tim_2[iSSY][iSFP][iTAM][iCHA] = MakeTH1 ('I', chis, chead, N_BIN_T, 0, N_BIN_T);
		}      
		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_tam; iCHA++)
		{
			sprintf (chis,"CALIB/SUB%d/SFP%d/TAMEX%02d/CHA%02d/CALIB SUM_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iCHA, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"CALIB SUM");
			h_sum_2[iSSY][iSFP][iTAM][iCHA] = MakeTH1 ('I', chis, chead, N_BIN_T, 0, N_BIN_T);
		}      
#ifdef IDATEN_MONITOR
		sprintf (chis,"MULTIPLICITY_LaBr3_hit");
		sprintf (chead,"MULTIPLICITY_LaBr3_hit");
		h1_Multiplicity_LaBr3 = MakeTH1 ('I', chis, chead, 20, 0, 20);

		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_tam; iCHA++)
		{
			sprintf (chis,"MULTIPLICITY0/SUB%d/SFP%d/TAMEX%02d/MULTIPLICITY0_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"MULTIPLICITY0");
			h1_Multiplicity[iSSY][iSFP][iTAM][iCHA][0] = MakeTH1 ('I', chis, chead, 10, 0, 10);
			sprintf (chis,"MULTIPLICITY1/SUB%d/SFP%d/TAMEX%02d/MULTIPLICITY1_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"MULTIPLICITY1");
			h1_Multiplicity[iSSY][iSFP][iTAM][iCHA][1] = MakeTH1 ('I', chis, chead, 10, 0, 10);
			sprintf (chis,"MULTIPLICITY2/SUB%d/SFP%d/TAMEX%02d/MULTIPLICITY2_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"MULTIPLICITY2");
			h1_Multiplicity[iSSY][iSFP][iTAM][iCHA][2] = MakeTH1 ('I', chis, chead, 10, 0, 10);
		}


		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_phy; iCHA++)
		{
			sprintf (chis,"By_PCha/SUB%d/SFP%d/TAMEX%02d/CHA%02d/SlowTOT_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iCHA, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"SlowTOT");
			h1_STOT[iSSY][iSFP][iTAM][iCHA] = MakeTH1 ('I', chis, chead,
					COARSE_CT_RANGE/2, 0, COARSE_CT_RANGE*CYCLE_TIME/4
					);
		}      
		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_phy; iCHA++)
		{
			sprintf (chis,"By_PCha/SUB%d/SFP%d/TAMEX%02d/CHA%02d/FastTOT_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iCHA, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"FastTOT");
			h1_FTOT[iSSY][iSFP][iTAM][iCHA] = MakeTH1 ('I', chis, chead,
					COARSE_CT_RANGE/4, 0, COARSE_CT_RANGE*CYCLE_TIME/4/8
					);
		}      
		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_phy; iCHA++)
		{
			sprintf (chis,"By_PCha/SUB%d/SFP%d/TAMEX%02d/CHA%02d/STOT_FTOT_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iCHA, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"h2_STOT_FTOT");
			h2_STOT_FTOT[iSSY][iSFP][iTAM][iCHA] = MakeTH2 ('I', chis, chead,
					COARSE_CT_RANGE/4/2, 0, COARSE_CT_RANGE*CYCLE_TIME/4,
					COARSE_CT_RANGE/2, 0, COARSE_CT_RANGE*CYCLE_TIME/4/8
					);
		}      
		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_phy; iCHA++)
		{
			sprintf (chis,"By_PCha/SUB%d/SFP%d/TAMEX%02d/CHA%02d/FTle-TTS_SUB%d_SFP%d_TAM%02d_CHA%02d", iSSY, iSFP, iTAM, iCHA, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"Fle-TTS");
			h1_FTle_TTS[iSSY][iSFP][iTAM][iCHA] = MakeTH1 ('I', chis, chead,
					COARSE_CT_RANGE, -COARSE_CT_RANGE*CYCLE_TIME/2, COARSE_CT_RANGE*CYCLE_TIME/2
					);
		}      
		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) 
		{
			sprintf (chis,"By_PCha/SUB%d/SFP%d/TAMEX%02d/PCHA_SlowTOT_SUB%d_SFP%d_TAM%02d", iSSY, iSFP, iTAM, iSSY, iSFP, iTAM);
			sprintf (chead,"PCHA_SlowTOT");
			h2_PCHA_STOT[iSSY][iSFP][iTAM] = MakeTH2 ('I', chis, chead,
					MAX_CHA_phy, 0, MAX_CHA_phy, 
					COARSE_CT_RANGE/2, 0, COARSE_CT_RANGE*CYCLE_TIME/4
					);
		}      
		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) 
		{
			sprintf (chis,"By_PCha/SUB%d/SFP%d/TAMEX%02d/PCHA_FastTOT_SUB%d_SFP%d_TAM%02d", iSSY, iSFP, iTAM, iSSY, iSFP, iTAM);
			sprintf (chead,"PCHA_FastTOT");
			h2_PCHA_FTOT[iSSY][iSFP][iTAM] = MakeTH2 ('I', chis, chead,
					MAX_CHA_phy, 0, MAX_CHA_phy, 
					COARSE_CT_RANGE/4, 0, COARSE_CT_RANGE*CYCLE_TIME/4/8
					);
		}
		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_phy; iCHA++)
		{
			sprintf (chis,"TREND/SUB%d/SFP%d/TAMEX%2d/CHA%2d/trendSTOT SUB %d SFP %d TAM %2d CHA %2d", iSSY, iSFP, iTAM, iCHA, iSSY, iSFP, iTAM, iCHA);
			sprintf (chead,"trend_STOT");
			h2_trend_STOT[iSSY][iSFP][iTAM][iCHA] = MakeTH2 ('I', chis, chead,
					TREND_N, 0, TREND_N,
					COARSE_CT_RANGE/4, 0, COARSE_CT_RANGE*CYCLE_TIME/4
					);
		}
#endif // IDATEN_MONITOR
	

		// box histograms for test channels only

/*		l_nc = MAX_CHA_old_AN >>2; // nr. of columns in picture
		l_nr = 4;              // nr. of rows    in picture
		fPicture = new TGo4Picture("Box test channels","Box");
		fPicture->SetDivision(l_nr, l_nc);
		fPicture->SetDrawHeader(kTRUE);
		l_c = 0;
		l_r = 0;
		for (l_i=0; l_i<MAX_CHA_old_AN; l_i++)
		{
			printf ("box l_r: %d, l_c: %d \n", l_r, l_c);
			fPicture->Pic(l_r,l_c)->AddObject(h_tim[l_i]);
			sprintf (c_tmp, "%s", &c_m_c[l_i][0]);
			printf ("box string: %s \n", c_tmp);
			TLatex* l = new TLatex(0.82, 0.13, c_tmp);
			l->SetNDC(kTRUE);
			//l->SetTextSize (30); // geht nicht??
			fPicture->Pic(l_r,l_c)->AddSpecialObject(l);
			if (l_c < (l_nc-1))
			{
				l_c++;
			}
			else
			{
				l_c = 0;
				l_r++;
			} 
		}  
		AddPicture(fPicture);
*/

		cout << "**** TTamex_FullProc: Created histograms and pictures" << endl;

		// open output file 
		if ((fd_out = fopen ("./tamex1_tdc_time_diff.txt", "a")) == NULL)
		{
			perror ("fopen");
			exit (0);
		}
		else
		{
			printf ("opened file: %s \n", "time_diff.txt"); 
		}
		fprintf (fd_out," %s", "Date");
		fprintf (fd_out, "%s", "\n\n");
		gettimeofday (&s_time, NULL);
		time (&l_time);

		for (iSSY=0; iSSY<MAX_SSY; iSSY++)
		{
			for (iSFP=0; iSFP<MAX_SFP; iSFP++)
			{
				for (iTAM=0; iTAM<MAX_TAM; iTAM++)
				{
					for (iCHA=0; iCHA<MAX_CHA_old; iCHA++)
					{
						l_err_ct [iSSY][iSFP][iTAM][iCHA] = 0;
						l_hit_ct [iSSY][iSFP][iTAM][iCHA] = 0;
					}
				}
			}
		}
	}
	else // got them from autosave file, restore pointers
	{
	}
}
//-----------------------------------------------------------
// event function
Bool_t TTamex_FullProc::BuildEvent(TGo4EventElement* target)
{  // called by framework. We dont fill any output event here at all
	Int_t      l_h, l_i, l_j, l_k;
	Int_t     size;
	Int_t     iSSY, iSFP, iTAM, iCHA, iTCHA, iPCHA;
	Bool_t    SlowFast; // 0 for fast, 1 for slow
						//UInt_t    *pl_se_dat;
						//UInt_t    *pl_tmp;
	uint32_t    *pl_se_dat;
	uint32_t    *pl_tmp;

#ifdef WR_TIME_STAMP
	uint32_t	l_wr_id;
	ULong64_t	l_wr_value;
	ULong64_t	l_wr_ts;
#endif // WR_TIME_STAMP
	UInt_t     l_padd;
	UInt_t     l_trig_type;
	UInt_t     l_ssy_idx;
	UInt_t     l_sfp_idx;
	UInt_t     l_tam_idx;
	UInt_t     l_cha_idx;
	UInt_t     l_tdc_dat; // 4bytes=32bits

	UInt_t     l_cha_head;  
	UInt_t     l_dat_size;
	UInt_t     l_tdc_head=0;
	UInt_t     l_tdc_size;
	UInt_t     l_tdc_trail=0;

	UInt_t     l_dat_len;  
	UInt_t     l_dat_len_byte;  

	UInt_t     l_dat;
	//  UInt_t     l_ix;

	Int_t     l_ch_tim;
	Int_t     l_ch_ix;
	Int_t     l_edge_type;
	Int_t     l_coarse_ct;

	Int_t	l_hct2[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_tam][3];
	Int_t   l_hct2_LaBr3;

	std::vector<UInt_t> 	v_SSY;
	std::vector<UInt_t> 	v_SFP;
	std::vector<UInt_t> 	v_TAM;
	std::vector<Int_t> 	v_TCHA;
	std::vector<Int_t> 	v_Edge_type;
	std::vector<Int_t> 	v_Coarse_ct;
	std::vector<Int_t> 	v_tdl;
	Double_t     d_tts;

	std::vector<UInt_t> 	v_SSY2;
	std::vector<UInt_t> 	v_SFP2;
	std::vector<UInt_t> 	v_TAM2;
	std::vector<Int_t> 	  v_PCHA;
	std::vector<Double_t> v_TOT;
	std::vector<Int_t> 	  v_Tle_cct;
	std::vector<Double_t> v_Fine_time;
	std::vector<Bool_t>   v_branch;
	std::vector<Double_t> v_TTS;

#ifdef DUMP_BAD_EVENT
	UInt_t     *pl_tdc_data=0;
	UInt_t     l_tdc_data_cnt=0;
#endif // DUMP_BAD_EVENT

	static Int_t l_bad_evt_found=0;

	static Double_t  d_diff;
	static Int_t     l_coarse_diff;

	static Double_t  d_finetimecal[MAX_SSY][MAX_SFP][MAX_TAM][MAX_CHA_tam][N_BIN_T];

	static ULong64_t  l_evt_ct=0;
	static ULong64_t  l_phy_evt_ct=0;
	static UInt_t     l_phy_tr_i=0;

	static UInt_t     l_trg_wind;
	static UInt_t     l_pre;
	static UInt_t     l_post;
	static UInt_t     l_trg_wind_len;
	static UInt_t     l_first=1;
	static UInt_t     l_first_tr=1;

	TGo4MbsSubEvent* psubevt;

	fInput = (TGo4MbsEvent* ) GetInputEvent();
	if(fInput == 0)
	{
		cout << "AnlProc: no input event !"<< endl;
		return kFALSE;
	}
	if(fInput->GetTrigger() > 11)
	{
		cout << "**** TTamex_FullProc: Skip trigger event"<<endl;
		return kFALSE;
	} 
	// Note that one has to loop over all subevents and select them by
	// crate number:   psubevt->GetSubcrate(),
	// procid:         psubevt->GetProcid(),
	// and/or control: psubevt->GetControl()
	// here we use only crate number

	fOutput= dynamic_cast<TTamex_FullEvent*> (target);
	if(fOutput == 0)
	{
		cout << "AnlProc: no ouput event !"<< endl;
		return kFALSE;
	}


	// start new calibration
	if(fPar->resetCalibration)
	{
		printf ("found request for new calibration \n");
		fPar->resetCalibration=kFALSE;
		fCalibrationDone=kFALSE;
		l_phy_tr_i = 0; //  reset trending bin, begin from 0

		for (iSSY=0; iSSY<MAX_SSY; iSSY++)
		{
			for (iSFP=0; iSFP<MAX_SFP; iSFP++)
			{
				for (iTAM=0; iTAM<MAX_TAM; iTAM++)
				{
					for (iCHA=0; iCHA<MAX_CHA_old; iCHA++)
					{
						l_err_ct [iSSY][iSFP][iTAM][iCHA] = 0;
						l_hit_ct [iSSY][iSFP][iTAM][iCHA] = 0;
					}
				}
			}
		}

		l_phy_evt_ct = 0;
		printf ("clear all histograms and event counter \n");

		//clear all histograms
		TGo4Analysis* an = TGo4Analysis::Instance();
		if (an) {
			an->ClearObjects("Histograms");
			// an->ClearObjects("Conditions"); // to clear conditions statistics
		}

		printf ("start collecting new calibration data \n");
		fflush (stdout);
	}

	l_evt_ct++;
	l_phy_evt_ct++; // no difference up to now

	//printf ("next     event \n"); sleep (1);

	for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_tam; iCHA++)
	{
		l_hct2[iSSY][iSFP][iTAM][iCHA][0]=0;
		l_hct2[iSSY][iSFP][iTAM][iCHA][1]=0;
		l_hct2[iSSY][iSFP][iTAM][iCHA][2]=0;
	}
	l_hct2_LaBr3 = 0;
	v_SSY       	.clear();
	v_SFP       	.clear();
	v_TAM       	.clear();
	v_TCHA       	.clear();
	v_Edge_type 	.clear();
	v_Coarse_ct 	.clear();
	v_tdl        	.clear();

	v_SSY2      .clear();
	v_SFP2      .clear();
	v_TAM2      .clear();
	v_PCHA      .clear();
	v_TOT       .clear();
	v_Tle_cct   .clear();
	v_Fine_time .clear();
	v_branch    .clear();
	v_TTS       .clear();

	l_err_catch = 0;
	l_num_err   = 0;

	l_bad_evt_found = 0; 
	fInput->ResetIterator();
	while((psubevt = fInput->NextSubEvent()) != 0) // loop over sub-events
	{
		//printf ("next sub-event \n"); sleep (1);

		l_ssy_idx = psubevt->GetProcid();
		if (l_ssy_idx == 100) {l_ssy_idx = 0;}
		else if (l_ssy_idx == 200) {l_ssy_idx = 1;}
		else {fprintf(stderr,"wrong Procid !!\n");}
		//l_ssy_idx = 0;

		//printf ("sub-event procid: %d, %d \n",  psubevt->GetProcid(), l_ssy_idx); fflush (stdout); 

		pl_se_dat = (uint32_t *)psubevt->GetDataField();
		pl_tmp = pl_se_dat;

#ifdef WR_TIME_STAMP
		// 5 first 32 bits must be white rabbit time stamp
		l_dat = *pl_tmp++;
		if (l_dat != SUB_SYSTEM_ID)
		{
			printf ("ERROR>> 1. data word is not sub-system id: \n");
			printf ("should be: 0x%x, but is: 0x%x\n", SUB_SYSTEM_ID, l_dat);
			goto bad_event;
		}
		l_wr_ts=0;

		l_dat = *pl_tmp++;
		l_wr_id = l_dat>>16;
		if (l_wr_id !=TS__ID_L16)
		{
			printf ("ERROR>> 2. data word does not contain 0-15 16bit identifier: \n");
			printf ("should be: 0x%x, but is: 0x%x\n", TS__ID_L16, l_wr_id);
		}
		l_wr_value = (l_dat & 0x0000ffff);
		l_wr_ts += l_wr_value << 0;

		l_dat = *pl_tmp++;
		l_wr_id = l_dat>>16;
		if (l_wr_id != TS__ID_M16)
		{
			printf ("ERROR>> 3. data word does not contain 16-31 16bit identifier: \n");
			printf ("should be: 0x%x, but is: 0x%x\n", TS__ID_M16, l_wr_id);
		}
		l_wr_value = (l_dat & 0x0000ffff);
		l_wr_ts += l_wr_value << 16;

		l_dat = *pl_tmp++;
		l_wr_id = l_dat>>16;
		if (l_wr_id != TS__ID_H16)
		{
			printf ("ERROR>> 4. data word does not contain 32-47 16bit identifier: \n");
			printf ("should be: 0x%x, but is: 0x%x\n", TS__ID_H16, l_wr_id);
		}
		l_wr_value = (l_dat & 0x0000ffff);
		l_wr_ts += l_wr_value << 32;

		l_dat = *pl_tmp++;
		l_wr_id = l_dat>>16;
		if (l_wr_id != TS__ID_X16)
		{
			printf ("ERROR>> 4. data word does not contain 48-63 16bit identifier: \n");
			printf ("should be: 0x%x, but is: 0x%x\n", TS__ID_H16, l_wr_id);
		}
		l_wr_value = (l_dat & 0x0000ffff);
		l_wr_ts += l_wr_value << 48;

		//fprintf(stdout,"l_wr_ts = 0x%llx = %llu\n", l_wr_ts, l_wr_ts);
		fOutput->Set_WR_TS(l_wr_ts);

#endif // WR_TIME_STAMP

		l_dat_len = psubevt->GetDlen();
		l_dat_len_byte = (l_dat_len - 2) * 2; 
		//fprintf(stdout,"l_dat_len %d l_dat_len_byte %d \n", l_dat_len, l_dat_len_byte);
		//printf("%x l_trg_wind\n",*pl_tmp);
		l_trg_wind = *pl_tmp++;

		if (l_first == 1)
		{
			l_pre  = ((l_trg_wind & 0xffff0000) >> 16);  // in clock cycles 
			l_post =  (l_trg_wind & 0x0000ffff);         // in clock cycles
			l_trg_wind_len = l_pre + l_post;

			printf ("trigger window: before trigger: %4d ns \n", l_pre  * 5);
			printf ("                after  trigger: %4d ns \n", l_post * 5);
			printf ("                total:        : %4d ns \n", (l_pre + l_post) * 5);

			l_first = 0;
		}

		//printf ("l_dat_len_byte: %d \n", l_dat_len_byte);
		for (l_i=0; l_i<100; l_i++)
		{
			//printf("%x l_padd\n",*pl_tmp);
			l_padd = *pl_tmp++;
			//printf ("l_padd: 0x%x \n", l_padd); 
			if ( (l_padd & 0xfff00000) != 0xadd00000 )
			{
				//printf ("%d padding words \n", l_i);
				pl_tmp--; 
				break;
			}
		}

		//sleep (1);

		// loop over all payload data 32 bit words
		while ( (pl_tmp - pl_se_dat) < (l_dat_len_byte/4) )
		{
			//printf("%x l_dat\n",*pl_tmp);
			l_dat = *pl_tmp++;   // must be padding word or channel header
								 //printf ("l_dat 0x%x \n", l_dat);

			if ( (l_dat & 0xff) == 0x34) //channel header
			{

				l_cha_head = l_dat;
				//printf ("l_cha_head: 0x%x \n", l_cha_head);

				//	0x34	= (l_cha_head & 0x000000ff) >>  0;
				l_trig_type 	= (l_cha_head & 0x00000f00) >>  8;
				l_sfp_idx   	= (l_cha_head & 0x0000f000) >> 12; // JAM indices for tree output event
				l_tam_idx   	= (l_cha_head & 0x00ff0000) >> 16;
				l_cha_idx   	= (l_cha_head & 0xff000000) >> 24; // shall always be 0 of tamex tdc
																   // tdc channels are coded in tdc data
																   //l_tdc_id   = l_tam_id;
																   //l_sfp_id    = 0;         // analysis is limited to one afp chain currently

				if (l_ssy_idx >= (MAX_SSY))
				{
					printf ("ERROR>> l_ssy_idx: %d \n", l_ssy_idx);  fflush (stdout);
					l_bad_evt_found = 1; 
					goto bad_event; 
				}
				if (l_sfp_idx >= (MAX_SFP))
				{
					printf ("ERROR>> l_spf_idx: %d \n", l_sfp_idx);  fflush (stdout);
					l_bad_evt_found = 1; 
					goto bad_event; 
				}
				if (l_tam_idx >= (MAX_TAM))
				{
					printf ("ERROR>> l_tam_idx: %d \n", l_tam_idx); fflush (stdout);
					l_bad_evt_found = 1; 
					goto bad_event; 
				}

				//printf("%x l_dat_size\n",*pl_tmp);
				l_dat_size = *pl_tmp++;

				//printf("%x l_tdc_head\n",*pl_tmp);
				l_tdc_head = *pl_tmp++;    // must be 0xaa header
				if ( ((l_tdc_head & 0xff000000) >> 24) != 0xaa)
				{
					printf ("ERROR>> header id is not 0xaa \n");
					l_bad_evt_found = 1; 
					goto bad_event; 
				}

				// now tdc data between 0xaa und 0xbb
				l_tdc_size = (l_dat_size/4) - 2;     // in longs/32bit
													 //printf("l_tdc_size %d l_dat_size %d\n", l_tdc_size, l_dat_size); // 9, 44

#ifdef DUMP_BAD_EVENT
				l_tdc_data_cnt = l_tdc_size;
				pl_tdc_data = pl_tmp;
#endif // DUMP_BAD_EVENT  	

				for (l_i=0; l_i< (Int_t) l_tdc_size; l_i++)   // loop over tdc data
				{
					l_err_flg = 0;
					//printf("%x l_tdc_dat\n",*pl_tmp);
					l_tdc_dat = *pl_tmp++;
					//printf ("raw tdc data: 0x%x\n", l_tdc_dat);
					//printf ("check %d \n", (l_tdc_dat & 0xe0000000) >> 29);

					if (((l_tdc_dat & 0xe0000000) >> 29) == 4)      // tdc channel data
					{
						//data_type	= (l_tdc_dat & 0xe0000000) >> 29;	// 4 for tdc 
						l_ch_ix		= (l_tdc_dat & 0x1fc00000) >> 22;	// 0x1fc: 7 bits: 0-127
						l_ch_tim	= (l_tdc_dat & 0x003ff000) >> 12;	// fine time
						l_edge_type	= (l_tdc_dat & 0x00000800) >> 11;	// 1 for leading edge, 0 for trailing edge
						l_coarse_ct	= (l_tdc_dat & 0x000007ff) >>  0;
						if (l_ch_ix > (MAX_CHA_old_INPUT-1))
						{
							goto bad_event; 
						}

						/*printf("l_sfp_idx %d l_tam_idx %d l_cha_idx %d l_ch_ix %d edge %d l_coarse_ct %d l_ch_tim %d\n", // TODO: save it to ttree!
						  l_sfp_idx,
						  l_tam_idx,
						  l_cha_idx,
						  (l_tdc_dat & 0x1fc00000) >> 22, 
						  (l_tdc_dat & 0x00000800) >> 11, 
						  (l_tdc_dat & 0x000007ff) >>  0,
						  (l_tdc_dat & 0x003ff000) >> 12
						  );*/

						//fOutput->AddFlipTime(l_ssy_idx, l_sfp_idx, l_tam_idx, l_ch_ix-1, l_edge_type, l_coarse_ct, l_ch_tim);
						v_SSY       	.push_back(l_ssy_idx);
						v_SFP       	.push_back(l_sfp_idx);
						v_TAM       	.push_back(l_tam_idx);
						v_TCHA       	.push_back(l_ch_ix-1);
						v_Edge_type 	.push_back(l_edge_type);
						v_Coarse_ct 	.push_back(l_coarse_ct);
						v_tdl        	.push_back(l_ch_tim);


						/*if(fCalibrationDone)
						  {
						  fprintf(stdout,"v_SSY        %u ", v_SSY       	.back());
						  fprintf(stdout,"v_SFP        %u ", v_SFP       	.back());
						  fprintf(stdout,"v_TAM        %u ", v_TAM       	.back());
						  fprintf(stdout,"v_TCHA       %d ", v_TCHA       .back());
						  fprintf(stdout,"v_Edge_type  %d ", v_Edge_type 	.back());
						  fprintf(stdout,"v_Coarse_ct  %d ", v_Coarse_ct 	.back());
						  fprintf(stdout,"v_tdl        %d ", v_tdl        .back());
						  fprintf(stdout,"\n");
						  fflush(stdout);
						  }*/

						// h_box is filled regardless if error or not
						h_box    [l_ssy_idx][l_sfp_idx][l_tam_idx][l_ch_ix]->Fill (l_ch_tim);
						l_hit_ct [l_ssy_idx][l_sfp_idx][l_tam_idx][l_ch_ix]++;
						if(l_ch_ix!=0 && l_edge_type==1)
						{
							h_cct [l_ssy_idx][l_sfp_idx][l_tam_idx][l_ch_ix-1]->Fill(l_coarse_ct);
							h_tim_2 [l_ssy_idx][l_sfp_idx][l_tam_idx][l_ch_ix-1]->Fill(l_ch_tim);
							//fprintf(stdout,"h_tim_2 [l_ssy_idx%u][l_sfp_idx%u][l_tam_idx%u][l_ch_ix-1 %d]->Fill(l_ch_tim%d);\n", l_ssy_idx,l_sfp_idx,l_tam_idx,l_ch_ix-1,l_ch_tim);
						}
						if(l_ch_ix!=0) l_hct2[l_ssy_idx][l_sfp_idx][l_tam_idx][l_ch_ix-1][l_edge_type]++;
						else if(l_ch_ix==0)
						{
							h_cct [l_ssy_idx][l_sfp_idx][l_tam_idx][MAX_CHA_tam-1]->Fill(l_coarse_ct);
							h_tim_2 [l_ssy_idx][l_sfp_idx][l_tam_idx][MAX_CHA_tam-1]->Fill(l_ch_tim);
							l_hct2[l_ssy_idx][l_sfp_idx][l_tam_idx][MAX_CHA_tam-1][l_edge_type]++;
						}

						//printf ("l_prev_num_err: %d \n", l_prev_num_err);
						if (l_prev_err_catch == 1)
						{
							for (l_j=0; l_j< (Int_t) l_prev_num_err; l_j++)
							{
								//printf ("err cha: %d %d %d %d \n", l_prev_err_ssy[l_j], l_prev_err_sfp[l_j],  l_prev_err_tam[l_j], l_prev_err_cha[l_j]);   

								if (    (l_ssy_idx == l_prev_err_ssy[l_j])
										&& (l_sfp_idx == l_prev_err_sfp[l_j]) 
										&& (l_tam_idx == l_prev_err_tam[l_j]) 
										&& (l_ch_ix   == (Int_t) l_prev_err_cha[l_j]) )
								{
									//printf ("bloedi \n"); fflush (stdout);

									h_err_box [l_ssy_idx][l_sfp_idx][l_tam_idx][l_ch_ix]->Fill (l_ch_tim);  
								}
							}
						}

						if (l_ch_tim == 0x3ff)                           // error in fpga fine time measurement 
						{
							l_err_flg = 1;
							l_err_ct[l_ssy_idx][l_sfp_idx][l_tam_idx][l_ch_ix]++;

							l_err_catch = 1;
							l_err_ssy[l_num_err] = l_ssy_idx;
							l_err_sfp[l_num_err] = l_sfp_idx;
							l_err_tam[l_num_err] = l_tam_idx;
							l_err_cha[l_num_err] = l_ch_ix;
							l_num_err++;
							if (l_num_err > MAX_CHA_old)
							{
								printf ("ERROR>> more 0x3ff errors found than storage prepared, exiting... \n");
								exit (0);
							}
						}
					} 
					else if (((l_tdc_dat & 0xe0000000) >> 29) == 0x3)   // channel epoch counter
					{
						//printf ("epoch found: 0x%x \n", (l_tdc_dat & 0xfffffff)); fflush (stdout);
					}
					else if (((l_tdc_dat & 0xff000000) >> 24) == 0xee)  // error word   
					{
					}
					else
					{
						printf ("ERROR>> unknown data type\t0x%x \n", l_tdc_dat);
					}
				}

				// tdc trailer
				//printf("%x l_tdc_trail\n",*pl_tmp);
				l_tdc_trail = *pl_tmp++;
				if ( ((l_tdc_trail & 0xff000000) >> 24) != 0xbb)
				{
					printf ("ERROR>> trailer id is not 0xbb, ");
					printf ("SUB: %d, SFP: %d, TAM: %d \n", l_ssy_idx, l_sfp_idx, l_tam_idx);
					l_bad_evt_found = 1; 
					goto bad_event; 
				}
			}
			else
			{
				printf ("ERROR>> data word neither channel header nor padding word \n");
			}       
		}
	}

	for (l_i=0; l_i< (Int_t) l_num_err; l_i++)
	{
		l_prev_err_ssy[l_i] = l_err_ssy[l_i]; 
		l_prev_err_sfp[l_i] = l_err_sfp[l_i]; 
		l_prev_err_tam[l_i] = l_err_tam[l_i]; 
		l_prev_err_cha[l_i] = l_err_cha[l_i]; 
	}
	l_prev_err_catch = l_err_catch;
	l_prev_num_err   = l_num_err;            

	if ( (l_evt_ct % STATISTIC) == 0)
	{
		printf ("\nfine time errors \n"); 
		for (l_h=0; l_h<MAX_SSY; l_h++)
		{
			for (l_i=0; l_i<MAX_SFP; l_i++)
			{
				for (l_j=0; l_j<MAX_TAM; l_j++)
				{
					for (l_k=0; l_k<MAX_CHA_old; l_k++)
					{
						if (l_err_ct [l_h][l_i][l_j][l_k] != 0)
						{
							d_err_rate = (Double_t)l_err_ct [l_h][l_i][l_j][l_k] / (Double_t)l_hit_ct [l_h][l_i][l_j][l_k];
							printf ("SUB: %d, SFP: %d, TAM: %2d, CHA: %2d: \t#errors: %5d \t#hits: %10d \t(%f) ",
									l_h, l_i, l_j, l_k, l_err_ct [l_h][l_i][l_j][l_k], l_hit_ct [l_h][l_i][l_j][l_k], d_err_rate);
							if ( (d_err_rate  > 0.01) && (d_err_rate  <= 0.02) )
							{
								printf ("  >>>> 1              error rate > 1 %%, <= 2 %% \n");
							}
							else if ( (d_err_rate  > 0.02) && (d_err_rate  <= 0.03) )
							{
								printf ("  >>>>>>>> 2          error rate > 2 %%, <= 3 %% \n");
							}
							else if ( (d_err_rate  > 0.03) && (d_err_rate  <= 0.04) )
							{
								printf ("  >>>>>>>>>>>> 3      error rate > 3 %%, <= 4 %% \n");
							}
							else if ( d_err_rate  > 0.04 )
							{
								printf ("  >>>>>>>>>>>>>>>> 4  error rate > 4 %% \n");
							}
							else
							{
								printf ("\n");
							}  
						}
					}
				}
			}
		}
	}

	//--------------------- do calibration ---------------------------------

	if (!fCalibrationDone) if( (l_phy_evt_ct >= N_CAL_EVT) || fPar->useOldCalibration)
	{
		printf ("charly! start calibaration \n");
		if(fPar->useOldCalibration)
		{
			fprintf(stdout,"fPar->useOldCalibration\n"); fflush(stdout);
			fprintf(fd_out,"fPar->useOldCalibration\n"); fflush(fd_out);
			/*for(iSSY=0; iSSY<MAX_SSY; iSSY++) for(iSFP=0; iSFP<MAX_SFP; iSFP++) for(iTAM=0; iTAM<MAX_TAM; iTAM++) for(iCHA=0; iCHA<MAX_CHA_tam; iCHA++)
			  {
			  for(l_i=0; l_i<N_BIN_T; l_i++)
			  d_finetimecal[iSSY][iSFP][iTAM][iCHA][l_i] = ((double) h_sum_2[iSSY][iSFP][iTAM][iCHA]->GetBinContent(l_i+1) / (double) h_tim_2[iSSY][iSFP][iTAM][iCHA]->GetEntries()) * CYCLE_TIME;
			  }*/
		}
		//else
		{
			for(iSSY=0; iSSY<MAX_SSY; iSSY++) for(iSFP=0; iSFP<MAX_SFP; iSFP++) for(iTAM=0; iTAM<MAX_TAM; iTAM++) for(iCHA=0; iCHA<MAX_CHA_tam; iCHA++)
			{
				Int_t sum=0;
				//fprintf(stdout,"h_tim_2[iSSY%u][iSFP%u][iTAM%u][iCHA%d]->GetEntries()==%.0f\n",iSSY,iSFP,iTAM,iCHA, h_tim_2[iSSY][iSFP][iTAM][iCHA]->GetEntries());
				if(h_tim_2[iSSY][iSFP][iTAM][iCHA]->GetEntries()<N_BIN_T*4) // 600*4
				{
					fprintf(stdout, "not enough entries to calibrate SSY%u SFP%u TAM%u TCHA%d %0.f ... set as 0.\n", iSSY, iSFP, iTAM, iCHA, h_tim_2[iSSY][iSFP][iTAM][iCHA]->GetEntries());
					for(l_i=0; l_i<N_BIN_T; l_i++) d_finetimecal[iSSY][iSFP][iTAM][iCHA][l_i] = 0;
					//for(l_i=0; l_i<N_BIN_T; l_i++) d_finetimecal[iSSY][iSFP][iTAM][iCHA][l_i] = (Double_t)(l_i * CYCLE_TIME)/N_BIN_T;
					//d_finetimecal[iSSY][iSFP][iTAM][iCHA]={0};
					continue;
				}
				else 
				{
					for(l_i=0; l_i<N_BIN_T; l_i++)
					{
						sum += (int)h_tim_2[iSSY][iSFP][iTAM][iCHA]->GetBinContent(l_i+1);
						//fprintf(stdout,"sum += h_tim_2[iSSY%u][iSFP%u][iTAM%u][iCHA%d]->GetBinContent(l_i+1 %d)%f;->%d\n", iSSY,iSFP,iTAM,iCHA,l_i,h_tim_2[iSSY][iSFP][iTAM][iCHA]->GetBinContent(l_i+1),sum);
						h_sum_2[iSSY][iSFP][iTAM][iCHA]->SetBinContent(l_i+1, sum);
						//d_finetimecal[iSSY][iSFP][iTAM][iCHA][l_i] = (double)sum / (double) h_tim_2[iSSY][iSFP][iTAM][iCHA]->GetEntries() * CYCLE_TIME;
					}
					if(sum!=h_tim_2[iSSY][iSFP][iTAM][iCHA]->GetEntries())
					{
						fprintf(stderr,"sum!=h_tim_2[iSSY][iSFP][iTAM][iCHA]->GetEntries()\n");
					}
					for(l_i=0; l_i<N_BIN_T; l_i++)
					{
						d_finetimecal[iSSY][iSFP][iTAM][iCHA][l_i] = (Double_t)h_sum_2[iSSY][iSFP][iTAM][iCHA]->GetBinContent(l_i+1) / (Double_t) sum * CYCLE_TIME;
					}
				}
			}
		}
#ifdef WR_TIME_STAMP
		l_wr_ts00 = l_wr_ts + TREND_INTV;
#endif // WR_TIME_STAMP
		fCalibrationDone=kTRUE;
		printf ("calibration finished \n");  
		fflush (stdout);

	}


	//------------------------ end of calibration --------------------------

	//------------- calculate times, fill histogramms ---------------------- 

	if(fCalibrationDone)
	{
#ifdef IDATEN_MONITOR
#ifdef WR_TIME_STAMP
		if (l_wr_ts > l_wr_ts00)
		{
			fprintf(stdout, "l_wr_ts %llu l_wr_ts00 %llu diff %llu d %llu\n", l_wr_ts, l_wr_ts00, l_wr_ts - l_wr_ts00, TREND_INTV); fflush(stdout);
			l_wr_ts00 += TREND_INTV;
			for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_phy; iCHA++) 
			{
				for (int ibiny=0; ibiny<COARSE_CT_RANGE/4; ibiny++) 
				{
					for (int ibinx=0; ibinx<TREND_N-1; ibinx++) 	
					{
						h2_trend_STOT[iSSY][iSFP][iTAM][iCHA]->SetBinContent(1+ibinx, 1+ibiny, h2_trend_STOT[iSSY][iSFP][iTAM][iCHA]->GetBinContent(1+ibinx+1, 1+ibiny));
					}
					h2_trend_STOT[iSSY][iSFP][iTAM][iCHA]->SetBinContent(TREND_N, 1+ibiny, 0);
				}
			}
		}

#endif // WR_TIME_STAMP
#endif // IDATEN_MONITOR


		size=v_SSY.size();
		//fprintf(stdout, "v_SSY.size() %d\n", size);
		for (int i=0; i<size; i++)
		{
			if(v_TCHA[i]==-1)
			{ 
				d_tts = (Double_t)(v_Coarse_ct[i] * CYCLE_TIME) - d_finetimecal[v_SSY[i]][v_SFP[i]][v_TAM[i]][MAX_CHA_tam-1][v_tdl[i]];
				//fprintf(stdout,"d_tts = (Double_t)(v_Coarse_ct[i%d]%d * CYCLE_TIME) - d_finetimecal[v_SSY[i]][v_SFP[i]][v_TAM[i]][MAX_CHA_tam-1][v_tdl[i]]=%.0f;\n",i, v_Coarse_ct[i],d_tts);
				continue;
			}
			if(v_tdl[i]==RESET_VAL) continue;
			for (int j=i; j<size; j++)
			{
				if(v_TCHA[j]==-1) continue;
				if(v_tdl[j]==RESET_VAL) continue;

				if(v_SSY[i]!=v_SSY[j]) continue;
				if(v_SFP[i]!=v_SFP[j]) continue;
				if(v_TAM[i]!=v_TAM[j]) continue;
				if(v_TCHA[i]!=v_TCHA[j]) continue;
				if(v_Edge_type[i]==1 && v_Edge_type[j]==0)
				{
					iSSY = v_SSY[i];
					iSFP = v_SFP[i];
					iTAM = v_TAM[i];
					iTCHA = v_TCHA[i];
					iPCHA = iTCHA/2;
					SlowFast = iTCHA%2; // 0 for fast, 1 for slow
					l_hct2[iSSY][iSFP][iTAM][iTCHA][2]++;

					l_coarse_diff = -v_Coarse_ct[i] + v_Coarse_ct[j];
					if(l_coarse_diff<0) l_coarse_diff += COARSE_CT_RANGE;
					d_diff = (Double_t)(l_coarse_diff * CYCLE_TIME) + d_finetimecal[iSSY][iSFP][iTAM][iTCHA][v_tdl[i]] - d_finetimecal[iSSY][iSFP][iTAM][iTCHA][v_tdl[j]];

					v_SSY2      .push_back(iSSY);
					v_SFP2      .push_back(iSFP);
					v_TAM2      .push_back(iTAM);
					v_PCHA      .push_back(iPCHA);
					v_TOT       .push_back(d_diff);
					v_Tle_cct   .push_back(v_Coarse_ct[i]);
					v_Fine_time .push_back(d_finetimecal[iSSY][iSFP][iTAM][iTCHA][v_tdl[i]]);
					v_branch    .push_back(SlowFast);
					v_TTS       .push_back(d_tts);

					/*fprintf(stdout,"v_SSY2      %u ", v_SSY2      .back());
					  fprintf(stdout,"v_SFP2      %u ", v_SFP2      .back());
					  fprintf(stdout,"v_TAM2      %u ", v_TAM2      .back());
					  fprintf(stdout,"v_TCHA      %d ", v_TAM[i]);
					  fprintf(stdout,"v_PCHA      %d ", v_PCHA      .back());
					  fprintf(stdout,"v_TOT       %f ", v_TOT       .back());
					  fprintf(stdout,"v_Tle_cct   %d ", v_Tle_cct   .back());
					  fprintf(stdout,"v_Fine_time %f ", v_Fine_time .back());
					  fprintf(stdout,"v_branch    %d ",(int) v_branch    .back());
					  fprintf(stdout,"\n");*/

#ifdef IDATEN_MONITOR
					if(SlowFast)
					{
						h1_STOT[iSSY][iSFP][iTAM][iPCHA]->Fill(d_diff);
						h2_trend_STOT[iSSY][iSFP][iTAM][iPCHA]->Fill(TREND_N-1,d_diff);
						h2_PCHA_STOT[iSSY][iSFP][iTAM]->Fill(iPCHA,d_diff);
					}
					else
					{
						h1_FTOT[iSSY][iSFP][iTAM][iPCHA]->Fill(d_diff);
						h2_PCHA_FTOT[iSSY][iSFP][iTAM]->Fill(iPCHA,d_diff);
					}
#endif // IDATEN_MONITOR
					break;
				}
			}
		}


		size=v_SSY2.size();
		//fprintf(stdout, "v_SSY2.size() %d\n", size);
		for (int itot=0; itot<size; itot++)
		{
			for (int jtot=0; jtot<size; jtot++)
			{
				if(v_SSY2[itot]!=v_SSY2[jtot]) continue;
				if(v_SFP2[itot]!=v_SFP2[jtot]) continue;
				if(v_TAM2[itot]!=v_TAM2[jtot]) continue;
				if(v_PCHA[itot]!=v_PCHA[jtot]) continue;

				if(v_branch[itot]==1 && v_branch[jtot]==0)
				{
					l_coarse_diff = v_Tle_cct[itot]-v_Tle_cct[jtot];
					if(l_coarse_diff<0) l_coarse_diff += COARSE_CT_RANGE;
					if(0<l_coarse_diff && l_coarse_diff<12)
					{
						//fprintf(stdout,"l_coarse_diff %d, h2_STOT_FTOT[v_SSY2[itot] %u][v_SFP[itot] %u][v_TAM2[itot] %u][v_PCHA[itot] %d]->Fill(v_TOT[itot] %f, v_TOT[jtot] %f);\n",
						//        l_coarse_diff, v_SSY2[itot], v_SFP[itot], v_TAM2[itot], v_PCHA[itot], v_TOT[itot], v_TOT[jtot] );
						h2_STOT_FTOT[v_SSY2[itot]][v_SFP2[itot]][v_TAM2[itot]][v_PCHA[itot]]->Fill(v_TOT[itot], v_TOT[jtot]);

						fOutput->AddHit(v_SSY2[itot], v_SFP2[itot], v_TAM2[itot], v_PCHA[itot],
								v_TOT[itot], (Double_t)(v_Tle_cct[itot]*CYCLE_TIME)-v_Fine_time[itot],
								v_TOT[jtot], (Double_t)(v_Tle_cct[jtot]*CYCLE_TIME)-v_Fine_time[jtot],
								v_TTS[itot]
								);
						if ( v_TAM2[itot]==0 || v_TAM2[itot]==1 || (v_TAM2[itot]==2 && v_PCHA[itot]<4) )
							l_hct2_LaBr3++;
						Double_t ftle_tts = (Double_t)(v_Tle_cct[jtot]*CYCLE_TIME)-v_Fine_time[jtot] - v_TTS[itot];
						if(ftle_tts < -CYCLE_TIME*COARSE_CT_RANGE/2) ftle_tts += CYCLE_TIME*COARSE_CT_RANGE;
						if(ftle_tts >  CYCLE_TIME*COARSE_CT_RANGE/2) ftle_tts -= CYCLE_TIME*COARSE_CT_RANGE;
						h1_FTle_TTS[iSSY][iSFP][iTAM][iPCHA]->Fill(ftle_tts);
						//h1_FTle_TTS[iSSY][iSFP][iTAM][iPCHA]->Fill((Double_t)(v_Tle_cct[jtot]*CYCLE_TIME)-v_Fine_time[jtot] - v_TTS[itot]);
					}
				}
			}
		}


#ifdef IDATEN_MONITOR
		for (iSSY=0; iSSY<MAX_SSY; iSSY++) for (iSFP=0; iSFP<MAX_SFP; iSFP++) for (iTAM=0; iTAM<MAX_TAM; iTAM++) for (iCHA=0; iCHA<MAX_CHA_tam; iCHA++)
		{
			h1_Multiplicity[iSSY][iSFP][iTAM][iCHA][0]->Fill(l_hct2[iSSY][iSFP][iTAM][iCHA][0]);
			h1_Multiplicity[iSSY][iSFP][iTAM][iCHA][1]->Fill(l_hct2[iSSY][iSFP][iTAM][iCHA][1]);
			h1_Multiplicity[iSSY][iSFP][iTAM][iCHA][2]->Fill(l_hct2[iSSY][iSFP][iTAM][iCHA][2]);
		}
		h1_Multiplicity_LaBr3->Fill(l_hct2_LaBr3);

#endif // IDATEN_MONITOR
	





	}


	// trending
	if((fCalibrationDone)
			&& ((l_phy_evt_ct % N_PHY_TREND_PRINT) == 0) )
	{
		gettimeofday (&s_time, NULL);
		//printf ("gettime: %d \n", (int)s_time.tv_sec); 

		fprintf (fd_out, "%d", (int)s_time.tv_sec); 



	}
	fflush (fd_out);
	fflush (stdout);

bad_event:

#ifdef DUMP_BAD_EVENT
	if (l_bad_evt_found == 1)
	{
		printf ("ERROR>> found bad event: tdc header:  0x%x \n", l_tdc_head);
		for (l_i=0; l_i<(Int_t) l_tdc_data_cnt; l_i++)
		{
			printf (" 0x%x ", *pl_tdc_data++);
			if (((l_i+1) % 8) == 0) { printf ("\n");}
		}
		if (((l_i+1) % 8) != 0) { printf ("\n");}
		printf ("                         tdc trailer: 0x%x \n\n", l_tdc_trail);
		fflush (stdout);
		//sleep (3);
	}
#endif // DUMP_BAD_EVENT

	return kTRUE;
}

//----------------------------END OF GO4 SOURCE FILE ---------------------
