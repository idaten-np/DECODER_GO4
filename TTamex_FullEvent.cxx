// $Id: TTamex_FullEvent.cxx 2627 2019-10-01 08:02:45Z linev $
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

#include "TTamex_FullEvent.h"

#include "TGo4Log.h"


//***********************************************************
TTamex_FullEvent::TTamex_FullEvent() :
   TGo4EventElement()
{
   TGo4Log::Info("TTamex_FullEvent: Create instance");
}
//***********************************************************
TTamex_FullEvent::TTamex_FullEvent(const char* name) :
   TGo4EventElement(name)
{
   TGo4Log::Info("TTamex_FullEvent: Create instance %s", name);
}
//***********************************************************
TTamex_FullEvent::~TTamex_FullEvent()
{
   TGo4Log::Info("TTamex_FullEvent: Delete instance");
}

//-----------------------------------------------------------

void TTamex_FullEvent::Clear(Option_t *t)
{
#ifdef WR_TIME_STAMP
		WR_TS=0;
#endif // WR_TIME_STAMP
      /*fprintf(stdout, "hit_SSY.size() %zu\n",hit_SSY.size()); 
      fprintf(stdout, "hit_SFP.size() %zu\n",hit_SFP.size()); 
      fprintf(stdout, "hit_TAM.size() %zu\n",hit_TAM.size()); 
      fprintf(stdout, "hit_PCHA.size() %zu\n",hit_PCHA.size()); 
      fprintf(stdout, "hit_STOT.size() %zu\n",hit_STOT.size()); 
      fprintf(stdout, "hit_STle.size() %zu\n",hit_STle.size()); 
      fprintf(stdout, "hit_FTOT.size() %zu\n",hit_FTOT.size()); 
      fprintf(stdout, "hit_FTle.size() %zu\n",hit_FTle.size()); 
      fprintf(stdout, "hit_TTS.size() %zu\n",hit_TTS.size()); fflush(stdout);*/
      hit_SSY               .clear();
      hit_SFP               .clear();
      hit_TAM               .clear();
      hit_PCHA              .clear();
      hit_STOT              .clear();
      hit_STle              .clear();
      hit_FTOT              .clear();
      hit_FTle              .clear();
      hit_CalE              .clear();
      hit_TTS               .clear();



}

