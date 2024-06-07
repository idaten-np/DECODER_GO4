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
      
      void AddHit(UInt_t ssy, UInt_t sfp, UInt_t tam, Int_t pcha, Double_t stot, Double_t stle, Double_t ftot, Double_t ftle, Double_t energy, Double_t tts);
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



