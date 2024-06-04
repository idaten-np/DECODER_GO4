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


class TTamex_FullEvent : public TGo4EventElement {
   public:
      TTamex_FullEvent();
      TTamex_FullEvent(const char* name);
      virtual ~TTamex_FullEvent();

      /** Method called by the framework to clear the event element. */
      void Clear(Option_t *t="");

#ifdef WR_TIME_STAMP
	  void Set_WR_TS(ULong64_t ts)
	  {
		  WR_TS=ts;
	  }
#endif // WR_TIME_STAMP
      
      void AddHit(UInt_t ssy, UInt_t sfp, UInt_t tam, Int_t pcha, Double_t stot, Double_t stle, Double_t ftot, Double_t ftle, Double_t energy, Double_t tts )
      {
          hit_SSY   .push_back( ssy );
          hit_SFP   .push_back( sfp );
          hit_TAM   .push_back( tam );
          hit_PCHA  .push_back( pcha );
          hit_STOT  .push_back( stot );
          hit_STle  .push_back( stle );
          hit_FTOT  .push_back( ftot );
          hit_FTle  .push_back( ftle );
          hit_CalE  .push_back( energy );
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
	  size_t GetN() {return hit_SSY.size();}


   private:
#ifdef WR_TIME_STAMP
   		ULong64_t WR_TS;
#endif // WR_TIME_STAMP
      std::vector<UInt_t> 	hit_SSY;
      std::vector<UInt_t> 	hit_SFP;
      std::vector<UInt_t> 	hit_TAM;
      std::vector<Int_t>  	hit_PCHA;
      std::vector<Double_t> hit_STOT;
      std::vector<Double_t> hit_STle;
      std::vector<Double_t> hit_FTOT;
      std::vector<Double_t> hit_FTle;
      std::vector<Double_t> hit_CalE;
      std::vector<Double_t> hit_TTS;

      

   ClassDef(TTamex_FullEvent,1)
};
#endif //TEVENT_H



