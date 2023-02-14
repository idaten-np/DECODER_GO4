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

#define MAX_CHA_AN_DIFF    MAX_CHA_AN / 2
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
      void SetTimeDiff(UChar_t channel, Double_t value)
      {
        if(channel>=MAX_CHA_AN_DIFF) return;
        fTimeDiff[channel]=value;
      }
      Double_t GetTimeDiff(UInt_t channel)
      {
	      if(channel>=MAX_CHA_AN_DIFF) return -1; // TODO: proper error handling maybe..
	      return fTimeDiff[channel];
      }

      void AddFlipTime(UInt_t ssy, UInt_t sfp, UInt_t tam, Int_t cha, Int_t edge_type, Int_t coarse_ct, Int_t fine_time)
      {
	      SSY       	.push_back(ssy);
	      SFP       	.push_back(sfp);
	      TAM       	.push_back(tam);
	      CHA       	.push_back(cha);
	      Edge_type 	.push_back(edge_type);
	      Coarse_ct 	.push_back(coarse_ct);
	      Fine_time 	.push_back(fine_time);
		//printf("AddFlipTime size %lu capacity %lu\n", SSY.size(), SSY.capacity());
      }



      /* Total number of timestamps in buffer for specific id. to be used in readout loops of second analysis step*/
//      UInt_t NumTimestamps(UChar_t channel)
//        {
//          if(channel>=MAX_CHA_AN) return 0;
//          return fTimeStamp[channel].size();
//        }

   private:

      Double_t fTimeDiff[MAX_CHA_AN_DIFF];
      std::vector<UInt_t> 	SSY;
      std::vector<UInt_t> 	SFP;
      std::vector<UInt_t> 	TAM;
      std::vector<Int_t> 	CHA;
      std::vector<Int_t> 	Edge_type;
      std::vector<Int_t> 	Coarse_ct;
      std::vector<Int_t> 	Fine_time;



   ClassDef(TTamex_FullEvent,1)
};
#endif //TEVENT_H



