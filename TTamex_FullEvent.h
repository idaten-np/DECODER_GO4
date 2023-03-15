// $Id: TTamex_FullEvent.h 2627 2019-10-01 08:02:45Z linev $
//-----------------------------------------------------------------------
//       The GSI Online Offline Object Oriented (Go4) Project
//         Experiment Data Processing at EE department, GSI
//-----------------------------------------------------------------------
// Copyright (C) 2000- GSI Helmholtzzentrum fuer Schwerionenforschung GmbH
//                     Planckstr. 1, 64291 Darmstadt, Germany
// Contact:            http://go4.gsi.de
//-----------------------------------------------------------------------
// This software can be used under the license agreements as stated
// in Go4License.txt file which is part of the distribution.
//-----------------------------------------------------------------------

#ifndef TTAMEXXEVENT_H
#define TTAMEXXEVENT_H

#include "TGo4EventElement.h"
#include <vector>

#include "TTamex_FullProc.h"

#define MAX_CHA_old_AN_DIFF    MAX_CHA_old_AN / 2
// JAM 3-jun-2022: we only use the preselected "analysis" channels from event processor here. This is redefinition of TTamex_FullProc.h
// proper way would be to put all defines here, but we don't do this not to confuse legacy users...
#define MAX_HITS 1



class TTamex_FullEvent : public TGo4EventElement {
   public:
      TTamex_FullEvent();
      TTamex_FullEvent(const char* name);
      virtual ~TTamex_FullEvent();

      /** Method called by the framework to clear the event element. */
      void Clear(Option_t *t="");
      
      /* add new timstamp to buffer*/
     /* void SetTimeDiff(UChar_t channel, Double_t value)
      {
        if(channel>=MAX_CHA_old_AN_DIFF) return;
        fTimeDiff[channel]=value;
        fprintf(stdout,"fTimeDiff[%d]=%.0f;\n", (int) channel, value); fflush(stdout);
      }
      Double_t GetTimeDiff(UInt_t channel)
      {
	      if(channel>=MAX_CHA_old_AN_DIFF) return -1; // TODO: proper error handling maybe..
	      return fTimeDiff[channel];
      }*/

      void AddFlipTime(UInt_t ssy, UInt_t sfp, UInt_t tam, Int_t tcha, Int_t edge_type, Int_t coarse_ct, Int_t tdl)
      {
	      flip_SSY       	.push_back(ssy);
	      flip_SFP       	.push_back(sfp);
	      flip_TAM       	.push_back(tam);
	      flip_TCHA     	.push_back(tcha);
	      flip_Edge_type 	.push_back(edge_type);
	      flip_Coarse_ct 	.push_back(coarse_ct);
	      flip_TDL 	      .push_back(tdl);
        //printf("AddFlipTime size %lu capacity %lu\n", flip_SSY.size(), flip_SSY.capacity());
      }

      void AddHit(UInt_t ssy, UInt_t sfp, UInt_t tam, Int_t pcha, Double_t stot, Double_t stle, Double_t ftot, Double_t ftle, Double_t tts )
      {
          hit_SSY   .push_back( ssy );
          hit_SFP   .push_back( sfp );
          hit_TAM   .push_back( tam );
          hit_PCHA  .push_back( pcha );
          hit_STOT  .push_back( stot );
          hit_STle  .push_back( stle );
          hit_FTOT  .push_back( ftot );
          hit_FTle  .push_back( ftle );
          hit_TTS   .push_back( tts );

          /*fprintf(stdout,"size %zu\n", hit_SSY.size());
          fprintf(stdout,"hit_SSY   %u\n"  , hit_SSY   .back());
          fprintf(stdout,"hit_SFP   %u\n"  , hit_SFP   .back());
          fprintf(stdout,"hit_TAM   %u\n"  , hit_TAM   .back());
          fprintf(stdout,"hit_PCHA  %d\n"  , hit_PCHA  .back());
          fprintf(stdout,"hit_STOT  %.0f\n", hit_STOT  .back());
          fprintf(stdout,"hit_STle  %.0f\n", hit_STle  .back());
          fprintf(stdout,"hit_FTOT  %.0f\n", hit_FTOT  .back());
          fprintf(stdout,"hit_FTle  %.0f\n", hit_FTle  .back());
          fprintf(stdout,"hit_TTS   %.0f\n", hit_TTS   .back());
          fprintf(stdout,"\n"); fflush(stdout);*/

      }


      /* Total number of timestamps in buffer for specific id. to be used in readout loops of second analysis step*/
//      UInt_t NumTimestamps(UChar_t channel)
//        {
//          if(channel>=MAX_CHA_old_AN) return 0;
//          return fTimeStamp[channel].size();
//        }

   private:
      //Double_t fTimeDiff[MAX_CHA_old_AN_DIFF];
      std::vector<UInt_t> 	flip_SSY;
      std::vector<UInt_t> 	flip_SFP;
      std::vector<UInt_t> 	flip_TAM;
      std::vector<Int_t>  	flip_TCHA;
      std::vector<Int_t> 	  flip_Edge_type;
      std::vector<Int_t> 	  flip_Coarse_ct;
      std::vector<Int_t> 	  flip_TDL;

      std::vector<UInt_t> 	hit_SSY;
      std::vector<UInt_t> 	hit_SFP;
      std::vector<UInt_t> 	hit_TAM;
      std::vector<Int_t>  	hit_PCHA;
      std::vector<Double_t> hit_STOT;
      std::vector<Double_t> hit_STle;
      std::vector<Double_t> hit_FTOT;
      std::vector<Double_t> hit_FTle;
      std::vector<Double_t> hit_TTS;

      

   ClassDef(TTamex_FullEvent,1)
};
#endif //TEVENT_H



