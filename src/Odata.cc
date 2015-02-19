//////////////////////////////////////////////////////////////////////////////
//  Author : florent robinet (LAL - Orsay): robinet@lal.in2p3.fr
//////////////////////////////////////////////////////////////////////////////
#include "Odata.h"

ClassImp(Odata)

////////////////////////////////////////////////////////////////////////////////////
Odata::Odata(const int aChunkDuration, const int aSegmentDuration,
	     const int aOverlapDuration, const int aVerbosity){ 
////////////////////////////////////////////////////////////////////////////////////
 
  // save parameters
  ChunkDuration      = (int)fabs(aChunkDuration);
  SegmentDuration    = (int)fabs(aSegmentDuration);
  OverlapDuration    = (int)fabs(aOverlapDuration);
  fVerbosity         = aVerbosity;

  // adjust durations
  if(!OverlapDuration){
    cerr<<"Odata::Odata: the nominal overlap duration is incorrect. Change it to: ---> 2s"<<endl;
    OverlapDuration=2;
  }
  if(OverlapDuration%2){
    cerr<<"Odata::Odata: the nominal overlap duration is not even. Change it to: ---> "<<OverlapDuration+1<<endl;
    OverlapDuration++;
  }

  if(SegmentDuration<2*OverlapDuration){
    SegmentDuration=NextPowerOfTwo(2*OverlapDuration+1);
    cerr<<"Odata::Odata: the nominal segment duration is too short. Change it to: ---> "<<SegmentDuration<<endl;
  }
  
  if(ChunkDuration<SegmentDuration){
    ChunkDuration=SegmentDuration;
    cerr<<"Odata::Odata: the nominal chunk duration is too short. Change it to: ---> "<<ChunkDuration<<endl;    
  }

  if((ChunkDuration-OverlapDuration)%(SegmentDuration-OverlapDuration)){
    ChunkDuration=((ChunkDuration-OverlapDuration)/(SegmentDuration-OverlapDuration)+1)*(SegmentDuration-OverlapDuration)+OverlapDuration;
    cerr<<"Odata::Odata: the nominal chunk duration does not macth the segment/overlap structure. Change it to: ---> "<<ChunkDuration<<endl;
  }
 
  // print durations
  if(fVerbosity>1){
    cout<<"Odata::Odata: the nominal durations are:"<<endl;
    cout<<"              Chunk duration   = "<<ChunkDuration<<"s"<<endl;
    cout<<"              Segment duration = "<<SegmentDuration<<"s"<<endl;
    cout<<"              Overlap duration = "<<OverlapDuration<<"s"<<endl;
  }

  // no timing yet
  fSegments  = NULL;
  seg        = -1;
  ChunkStart = -1;
  ChunkStop  = -1;
  NSegments  = (ChunkDuration-OverlapDuration)/(SegmentDuration-OverlapDuration);
  OverlapDurationCurrent=OverlapDuration;
}

////////////////////////////////////////////////////////////////////////////////////
Odata::~Odata(void){
////////////////////////////////////////////////////////////////////////////////////
  if(fVerbosity>1) cout<<"Odata::~Odata"<<endl;
}

////////////////////////////////////////////////////////////////////////////////////
bool Odata::SetSegments(Segments *aSegments){
////////////////////////////////////////////////////////////////////////////////////
  if(aSegments==NULL||!aSegments->GetLiveTime()){
    cerr<<"Odata::SetSegments: no live time in input segments"<<endl;
    return false;
  }

  fSegments = aSegments;

  // reset timing
  seg        = 0; // sets on first segment
  ChunkStart = (int)(fSegments->GetStart(seg));
  ChunkStop  = -1;// flags the first chunk
  NSegments  = (ChunkDuration-OverlapDuration)/(SegmentDuration-OverlapDuration);
  OverlapDurationCurrent=OverlapDuration;

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Odata::NewChunk(void){
////////////////////////////////////////////////////////////////////////////////////
  if(ChunkStart<0){
    cerr<<"Odata::NewChunk: no segment to read"<<endl;
    return false;
  }
  
  
  // First chunk
  if(ChunkStop<0){
    ChunkStop=ChunkStart+ChunkDuration;
    OverlapDurationCurrent=ChunkStart+OverlapDuration;
  }
  // New chunk
  else{
    OverlapDurationCurrent=ChunkStop;
    ChunkStart=ChunkStop-OverlapDuration;
    ChunkStop=ChunkStart+ChunkDuration;
  }
  
  // test chunk against segments and update chunk if necessary
  while(!TestChunk());
  
  // overlap with previous chunk
  OverlapDurationCurrent-=ChunkStart;

  // correct if we start a new segment
  OverlapDurationCurrent=TMath::Max(OverlapDuration,OverlapDurationCurrent);

  // end of data segments --> stop
  if(seg>=fSegments->GetNsegments()){
    if(fVerbosity) cout<<"Odata::NewChunk: end of data segments"<<endl;
    ChunkStart=-1;
    return false;
  }
  
  // new values
  NSegments=(ChunkStop-ChunkStart-OverlapDuration)/(SegmentDuration-OverlapDuration);

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
bool Odata::TestChunk(void){
  
  // TEST1: the segment left-over is too small to be processed --> move to next segment
  if(ChunkStart>=fSegments->GetEnd(seg)-OverlapDuration){
    seg++;// move to next segment
    if(seg>=fSegments->GetNsegments()) return true;// end of data segments: stop the loop
    ChunkStart=(int)(fSegments->GetStart(seg));
    ChunkStop=ChunkStart+ChunkDuration;
    return false;// the validity of this case needs to be tested again
  }

  // TEST2: simple case: inside the current segment --> OK
  if(ChunkStart>=fSegments->GetStart(seg)&&ChunkStop<=fSegments->GetEnd(seg))
    return true;// stop the loop

  // TEST3: the chunk stops after the segment end (edge effect)
  if(ChunkStart<fSegments->GetEnd(seg)&&ChunkStop>fSegments->GetEnd(seg)){
    if(ChunkStop-ChunkStart>SegmentDuration)//  --> shorten chunk
      ChunkStop-=(SegmentDuration-OverlapDuration);
    else{//  --> left-over to be processed
      ChunkStop=(int)(fSegments->GetEnd(seg));
      ChunkStart=ChunkStop-SegmentDuration;
    }
    return false;// the validity of this case needs to be tested again
  }

  // TEST4: the chunk start before current segment (due to TEST3) --> move to next segment
  if(ChunkStart<fSegments->GetStart(seg)){
    seg++;// move to next segment
    if(seg>=fSegments->GetNsegments()) return true;// end of data segments
    ChunkStart=(int)(fSegments->GetStart(seg));
    ChunkStop=ChunkStart+ChunkDuration;
    return false;// the validity of this case needs to be tested again
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////////
int Odata::GetSegmentTimeStart(const int aNseg){
////////////////////////////////////////////////////////////////////////////////////
  if(aNseg<0||aNseg>=NSegments){
    cerr<<"Odata::GetSegmentTimeStart: segment "<<aNseg<<" cannot be found in the data chunk"<<endl;
    return -1;
  }

  return ChunkStart+aNseg*(SegmentDuration-OverlapDuration);
}

////////////////////////////////////////////////////////////////////////////////////
int Odata::GetSegmentTimeEnd(const int aNseg){
////////////////////////////////////////////////////////////////////////////////////
  if(aNseg<0||aNseg>=NSegments){
    cerr<<"Odata::GetSegmentTimeEnd: segment "<<aNseg<<" cannot be found in the data chunk"<<endl;
    return -1;
  }

  return ChunkStart+SegmentDuration+aNseg*(SegmentDuration-OverlapDuration);
}


