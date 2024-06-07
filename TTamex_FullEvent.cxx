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

void TTamex_FullEvent::AddHit(UInt_t ssy, UInt_t sfp, UInt_t tam, Int_t pcha, Double_t stot, Double_t stle, Double_t ftot, Double_t ftle, Double_t energy, Double_t tts )
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
}

void TTamex_FullEvent::Clear(Option_t *t)
{
#ifdef WR_TIME_STAMP
		WR_TS=0;
#endif // WR_TIME_STAMP
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

