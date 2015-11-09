//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Odata.h"

ClassImp(Odata)

////////////////////////////////////////////////////////////////////////////////////
Odata::Odata(const int aChunkDuration, const int aOverlapDuration, const int aVerbosity){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // save parameters
  ChunkDuration      = (int)fabs(aChunkDuration);
  OverlapDuration    = (int)fabs(aOverlapDuration);
  fVerbosity         = aVerbosity;
   
  // print durations
  if(fVerbosity>1){
    cout<<"Odata::Odata: the nominal durations are:"<<endl;
    cout<<"              Chunk duration   = "<<ChunkDuration<<"s"<<endl;
    cout<<"              Overlap duration = "<<OverlapDuration<<"s"<<endl;
  }

  // no timing yet
  fSegments              = new Segments();
  seg                    = 0;
  ChunkStart             = -1;
  OverlapDurationCurrent = OverlapDuration;
}

////////////////////////////////////////////////////////////////////////////////////
Odata::~Odata(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Odata::~Odata"<<endl;
  delete fSegments;
}

////////////////////////////////////////////////////////////////////////////////////
bool Odata::SetSegments(Segments *aSegments){
////////////////////////////////////////////////////////////////////////////////////
  if(aSegments==NULL||!aSegments->GetStatus()){
    cerr<<"Odata::SetSegments: invalid input segments"<<endl;
    return false;
  }

  // copy input segments
  delete fSegments;
  fSegments = new Segments(aSegments->GetStarts(),aSegments->GetEnds());

  // reset timing
  seg                    = 0; // sets on first segment
  ChunkStart             = -1;// means beginning segment
  OverlapDurationCurrent = OverlapDuration;

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Odata::NewChunk(void){
////////////////////////////////////////////////////////////////////////////////////
  if(seg>=fSegments->GetNsegments()){
    cerr<<"Odata::NewChunk: end of segments"<<endl;
    return false;
  }

  // current segment is too short --> move to next segment
  if((int)fSegments->GetLiveTime(seg)<ChunkDuration){
    seg++;
    ChunkStart=-1;
    return NewChunk();
  }

  // first chunk for this segment = initialization
  if(ChunkStart<0) ChunkStart=(int)fSegments->GetStart(seg)-ChunkDuration+OverlapDuration;

  // new test chunk
  OverlapDurationCurrent = OverlapDuration;// reset current overlap
  int start_test  = ChunkStart+ChunkDuration-OverlapDuration;
  int stop_test   = start_test + ChunkDuration;
  
  // end of current segment --> move to next segment
  if(start_test==(int)fSegments->GetEnd(seg)-OverlapDuration){
    seg++;
    ChunkStart=-1;
    return NewChunk();
  }
    
  // chunk ends after current segment end --> adjust overlap
  if(stop_test>(int)fSegments->GetEnd(seg)){
    ChunkStart=(int)fSegments->GetEnd(seg)-ChunkDuration;
    OverlapDurationCurrent=start_test+OverlapDuration-ChunkStart;
    seg++;
    return true;
  }

  ChunkStart=start_test;
  return true;
}



